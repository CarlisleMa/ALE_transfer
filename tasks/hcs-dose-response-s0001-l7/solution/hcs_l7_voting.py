#!/usr/bin/env python3
"""Level-7 reference solution: the level-4 analysis (classical image analysis without Cellpose) with one
step in front of the plate statistics: the plate controls are looked at before they are used.

  per image   exactly as at level 4 (hcs_l4_voting.py): radial-symmetry voting for the nuclei, SOP rules for
              fragments and the counting frame, voting parameters chosen on the annotated fields, histogram-
              mode background, ring 1-4 px outside each counted nucleus, power-spectrum focus metric.
              Nothing of this needs the plate layout: an image, its annotation file and its row in the two
              per-well tables carry the same name.
  controls    the layout is laid over the images in the two ways a plate fits into a stage holder (as the
              files are named; turned by 180 degrees, i.e. the image named for row r, column c shows plate
              well (9 - r, 13 - c)). For each: plate model on the wells the layout calls vehicle, then
                * the marker response of the wells the layout calls reference (the reference inhibitor
                  must be inhibited), and
                * the share of the vehicle wells that lie more than a factor 1.3 off the plate model on
                  either readout (vehicle wells must look alike).
              The orientation in which the controls behave is used; if neither or both do, the analysis
              stops instead of guessing.
  statistics  as at level 4 (edge ratio + Theil-Sen plate model on the vehicle wells, grid-search fit of the
              SOP model), with the layout in the orientation found.

The second solution (hcs_l7_cellpose.py) shares no code with this file and finds the orientation without the
controls (direction of the dose series).

The switches reproduce what a careless analyst would deliver; they exist so that the task's verifier can be
attacked with real outputs built on a correct base (see the family's battery).
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--orientation", choices=["auto", "name", "turned"], default="auto")   # name = layout joined by file name
ap.add_argument("--controls", choices=["layout", "swapped"], default="layout")         # swapped = vehicle <-> reference
ap.add_argument("--rows", choices=["image", "plate"], default="image")   # naming of the rows of the per-well tables
ap.add_argument("--measurements")    # per-image table (well, n_cells, marker_per_cell, qc_pass) instead of the images
ap.add_argument("--counting", choices=["sop", "fragments", "noframe"], default="sop")   # fragments = fragments counted
ap.add_argument("--qc", choices=["focus", "none"], default="focus")
ap.add_argument("--readout", choices=["percell", "wholewell"], default="percell")
ap.add_argument("--background", choices=["local", "global", "none"], default="local")
ap.add_argument("--plate", choices=["model", "none"], default="model")
ap.add_argument("--fit", choices=["sop", "free"], default="sop")
ap.add_argument("--dump", action="store_true")    # also write the responses of every well (responses.tsv)
ap.add_argument("--input", default="/task/input")
ap.add_argument("--output", default="/task/output")
args = ap.parse_args()

inp, out = Path(args.input), Path(args.output)
out.mkdir(parents=True, exist_ok=True)
lay = pd.read_csv(inp / "plate_layout.tsv", sep="\t")
lay["r"] = lay["well"].str[0].map(lambda s: ord(s) - 65)
lay["c"] = lay["well"].str[1:].astype(int) - 1
NR, NC = int(lay["r"].max()) + 1, int(lay["c"].max()) + 1
lay["edge"] = lay["r"].isin([0, NR - 1]) | lay["c"].isin([0, NC - 1])
names = list(lay["well"])                            # the images are named like the wells of the layout
manual = {p.stem: sum(1 for ln in p.read_text().splitlines()[1:] if ln.strip())
          for p in sorted((inp / "annotations").glob("*.tsv"))}


def turned(well):
    """The well at the same place of the holder when the plate lies the other way round."""
    return f"{chr(65 + NR - 1 - (ord(well[0]) - 65))}{NC + 1 - int(well[1:]):02d}"


# ---- per image: counts, marker per cell, focus (as at level 4) -------------------------------------
PX_UM = 1.3
FRAGMENT_AREA = np.pi / 4 * (9.5 / PX_UM) ** 2      # px^2 at half maximum; SOP: fragments <= 6.5 um, nuclei >= 11 um across
TOL_ABS, TOL_REL = 3, 0.10                          # SOP section 2


def measure_images():
    import tifffile
    from scipy import ndimage as ndi
    from skimage.feature import peak_local_max
    from skimage.segmentation import expand_labels, watershed

    dark = tifffile.imread(inp / "flatfield" / "dark.tif").astype(np.float64)
    gainmap = {}
    for ch in (1, 2):
        f = tifffile.imread(inp / "flatfield" / f"flat_ch{ch}.tif").astype(np.float64) - dark
        f = ndi.uniform_filter(f, 9)
        gainmap[ch] = f / np.median(f)

    def read(well, ch):
        return (tifffile.imread(inp / "images" / f"{well}_ch{ch}.tif").astype(np.float64) - dark) / gainmap[ch]

    def nuclear_level(im):
        sm = ndi.gaussian_filter(im, 1.0)
        b = np.median(sm)
        return float(np.median(sm[sm > b + 0.5 * (np.percentile(sm, 99) - b)]) - b)

    nuc_img = {w: read(w, 1) for w in names}
    # typical nuclear signal of the plate. Taken over all fields, not over the wells the layout calls vehicle:
    # nothing in the image analysis may depend on a layout that has not been checked yet.
    level = float(np.median([nuclear_level(nuc_img[w]) for w in names]))
    prep = {}

    def prepared(w):
        """Smoothed image, background, foreground, gradient and the size of the nuclei of this field relative
        to the usual 5.5 px half-width (treated nuclei can be larger; all radii scale with it)."""
        if w not in prep:
            sm = ndi.gaussian_filter(nuc_img[w], 1.0)
            bg = float(np.median(sm))
            fg = sm > bg + 0.25 * level
            comp, n = ndi.label(fg)
            scale = 1.0
            if n:
                idx = np.arange(1, n + 1)
                big = ndi.sum(fg, comp, idx) >= 60
                if big.sum() >= 3:
                    half_width = np.median(ndi.maximum(ndi.distance_transform_edt(fg), comp, idx)[big])
                    scale = float(np.clip((half_width - 0.8) / 5.5, 0.85, 2.0))
            gy, gx = np.gradient(sm)
            prep[w] = (sm, bg, fg, scale, gy, gx)
        return prep[w]

    def segment(w, edge_thr, vote_sigma):
        """Returns (labels, ids of counted nuclei, ids of all nuclei incl. those outside the counting frame)."""
        sm, bg, fg, scale, gy, gx = prepared(w)
        H, W = sm.shape
        mag = np.hypot(gy, gx)
        local = np.maximum(ndi.maximum_filter(sm, 7) - bg, 0.25 * level)
        edge = (mag / local > edge_thr) & ndi.binary_dilation(fg, iterations=2)
        ys, xs = np.nonzero(edge)
        uy, ux = gy[ys, xs] / mag[ys, xs], gx[ys, xs] / mag[ys, xs]
        votes = np.zeros((H, W))
        for r in np.arange(4.5 * scale, 7.0 * scale + 0.01, 0.5):
            vy, vx = np.round(ys + r * uy).astype(int), np.round(xs + r * ux).astype(int)
            ok = (vy >= 0) & (vy < H) & (vx >= 0) & (vx < W)
            np.add.at(votes, (vy[ok], vx[ok]), 1.0)
        votes = ndi.gaussian_filter(votes, vote_sigma * scale)
        pk = peak_local_max(votes, min_distance=max(2, int(round(4 * scale))), threshold_abs=3.0,
                            labels=fg.astype(int), exclude_border=False)
        lab = np.zeros((H, W), dtype=int)
        if not len(pk):
            return lab, [], []
        seeds = np.zeros((H, W), dtype=int)
        seeds[tuple(pk.T)] = np.arange(1, len(pk) + 1)
        lab = watershed(-votes, seeds, mask=fg)
        high = set(np.unique(np.concatenate([lab[-1, :], lab[:, -1]]))) - {0}     # bottom / right edge
        low = set(np.unique(np.concatenate([lab[0, :], lab[:, 0]]))) - {0}        # top / left edge
        peak = ndi.maximum(sm, lab, np.arange(1, len(pk) + 1))
        R = int(round(10 * scale))
        counted, nuclei = [], []
        for i, (y, x) in enumerate(pk, start=1):
            fragment = False
            if i not in low and i not in high and args.counting != "fragments":
                y0, y1, x0, x1 = max(0, y - R), min(H, y + R + 1), max(0, x - R), min(W, x + R + 1)
                reg, _ = ndi.label(sm[y0:y1, x0:x1] > bg + 0.5 * (peak[i - 1] - bg))
                hit = np.unique(reg[lab[y0:y1, x0:x1] == i])
                hit = hit[hit > 0]
                fragment = (not len(hit)) or np.isin(reg, hit).sum() < FRAGMENT_AREA
            if fragment:
                continue
            nuclei.append(i)
            if i not in high or args.counting == "noframe":
                counted.append(i)
        return lab, counted, nuclei

    # tuning on the annotated fields: the setting with the smallest worst-case error relative to the SOP limit
    best = None
    for edge_thr in (0.12, 0.15):
        for vote_sigma in (1.2, 1.5):
            worst = max(abs(len(segment(w, edge_thr, vote_sigma)[1]) - m) / max(TOL_ABS, TOL_REL * m)
                        for w, m in manual.items())
            if best is None or worst < best[0]:
                best = (worst, edge_thr, vote_sigma)
    print(f"voting parameters: edge {best[1]}, sigma {best[2]}; worst annotated field at {best[0]:.2f} of the SOP limit")

    def hist_mode(a):
        lo, hi = np.percentile(a, [0.5, 85])
        width = max((hi - lo) / 120.0, 0.5)
        h, e = np.histogram(a, bins=np.arange(lo, hi + width, width))
        h = ndi.uniform_filter1d(h.astype(float), 7)
        k = int(np.argmax(h))
        return 0.5 * (e[k] + e[k + 1])

    fy = np.fft.fftfreq(dark.shape[0])[:, None]
    fx = np.fft.fftfreq(dark.shape[1])[None, :]
    fr = np.hypot(fy, fx)
    LOW, HIGH, NOISE = (fr > 0.015) & (fr < 0.07), (fr > 0.14) & (fr < 0.30), fr > 0.42

    def sharpness(img):
        p = np.abs(np.fft.fft2(img - img.mean())) ** 2
        floor = p[NOISE].mean()
        low = p[LOW].mean() - floor
        if low < 30 * floor:
            return np.nan                              # (almost) empty field: cannot be judged, keep
        return max(p[HIGH].mean() - floor, 0.0) / low

    rec = []
    for w in names:
        m2 = read(w, 2)
        lab, counted, nuclei = segment(w, best[1], best[2])
        ids = np.array(counted, dtype=int)
        cells = np.where(np.isin(lab, ids), lab, 0)
        ring = np.where(expand_labels(lab, 1) > 0, 0, expand_labels(cells, 4))
        bg = hist_mode(m2)
        sums = ndi.sum(m2, ring, ids) if len(ids) else np.array([])
        area = ndi.sum(np.ones_like(m2), ring, ids) if len(ids) else np.array([])
        has = area > 0                                   # a nucleus enclosed by its neighbours has no ring
        raw = float(np.mean(sums[has] / area[has])) if has.any() else np.nan
        rec.append((w, len(ids), raw, bg, float(m2.mean()), sharpness(nuc_img[w])))
    m = pd.DataFrame(rec, columns=["image", "n_cells", "raw_ring", "bg", "field_mean", "sharp"])
    m["qc_pass"] = np.where(m["sharp"] < 0.35 * np.nanmedian(m["sharp"]), 0, 1) if args.qc == "focus" else 1
    return m


if args.measurements:        # per-image values handed in (used to attack the verifier with planted measurements)
    given = pd.read_csv(args.measurements, sep="\t", na_values=["NA"], keep_default_na=True)
    meas = pd.DataFrame({"image": given["well"], "n_cells": given["n_cells"].astype(int),
                         "raw_ring": given["marker_per_cell"].astype(float), "bg": 0.0,
                         "field_mean": given["marker_per_cell"].astype(float), "sharp": np.nan,
                         "qc_pass": given["qc_pass"].astype(int)})
    meas = meas.set_index("image").loc[names].reset_index()
else:
    meas = measure_images()

# background of the marker channel: a property of the image, like the count and the focus state
img_good = meas["qc_pass"] == 1
if args.background == "local":
    meas["marker_per_cell"] = meas["raw_ring"] - meas["bg"]
    meas["field_signal"] = meas["field_mean"] - meas["bg"]
elif args.background == "none":
    meas["marker_per_cell"] = meas["raw_ring"]
    meas["field_signal"] = meas["field_mean"]


# ---- the layout laid over the images ------------------------------------------------------------------
def theil_sen(x, y):
    i, j = np.triu_indices(len(x), 1)
    dx = x[j] - x[i]
    keep = dx != 0
    return float(np.median((y[j] - y[i])[keep] / dx[keep]))


def joined(orientation):
    """Layout rows (plate wells) with the measurements of the image that shows each of them."""
    t = lay.copy()
    t["image"] = t["well"] if orientation == "name" else t["well"].map(turned)
    t = t.merge(meas, on="image", how="left", validate="one_to_one")
    if args.controls == "swapped":                   # the two control roles exchanged, nothing else
        t["role"] = t["role"].map({"vehicle": "reference", "reference": "vehicle"}).fillna(t["role"])
    if args.background == "global":                  # one constant for the plate, from the vehicle wells
        g = float(np.median(t.loc[(t["role"] == "vehicle") & (t["qc_pass"] == 1), "bg"]))
        t["marker_per_cell"] = t["raw_ring"] - g
        t["field_signal"] = t["field_mean"] - g
    t["marker_readout"] = t["marker_per_cell"] if args.readout == "percell" else t["field_signal"]
    return t


def vehicle_expectation(t, col, plate_model=True):
    v = t[(t["role"] == "vehicle") & (t["qc_pass"] == 1) & (t[col] > 0)]
    if not plate_model:
        return np.full(len(t), float(v[col].mean()))
    y = np.log(v[col].to_numpy(dtype=float))
    c, r, e = v["c"].to_numpy(float), v["r"].to_numpy(float), v["edge"].to_numpy(bool)
    bc = br = le = 0.0
    for _ in range(6):
        res = y - bc * c - br * r
        le = (np.median(res[e]) - np.median(res[~e])) if e.any() and (~e).any() else 0.0
        z = y - le * e
        bc = theil_sen(c, z - br * r)
        br = theil_sen(r, z - bc * c)
    a = np.median(y - le * e - bc * c - br * r)
    return np.exp(a + le * t["edge"].to_numpy(float) + bc * t["c"].to_numpy(float) + br * t["r"].to_numpy(float))


def control_check(orientation):
    """How the wells the layout calls vehicle and reference look when the layout is applied this way."""
    t = joined(orientation)
    good = t["qc_pass"] == 1
    em, en = vehicle_expectation(t, "marker_per_cell"), vehicle_expectation(t, "n_cells")
    rm, rn = t["marker_per_cell"] / em, t["n_cells"] / en
    veh, ref = (t["role"] == "vehicle") & good, (t["role"] == "reference") & good
    off = ((np.abs(np.log(np.clip(rm[veh], 1e-3, None))) > np.log(1.3))
           | (np.abs(np.log(np.clip(rn[veh], 1e-3, None))) > np.log(1.3)))
    return {"orientation": orientation, "reference_wells": int(ref.sum()), "vehicle_wells": int(veh.sum()),
            "reference_marker_response_median": float(np.nanmedian(rm[ref])),
            "reference_wells_below_half": int((rm[ref] < 0.5).sum()),
            "vehicle_wells_off_plate_model": int(off.sum()),
            "vehicle_marker_min": float(rm[veh].min()), "vehicle_marker_max": float(rm[veh].max()),
            "vehicle_count_min": float(rn[veh].min()), "vehicle_count_max": float(rn[veh].max())}


def behaves(c):
    return c["reference_marker_response_median"] < 0.5 and c["vehicle_wells_off_plate_model"] <= 0.15 * c["vehicle_wells"]


checks = [control_check(o) for o in ("name", "turned")]
pd.DataFrame(checks).to_csv(out / "control_check.tsv", sep="\t", index=False, float_format="%.4g")
print(pd.DataFrame(checks).T.to_string(header=False))
if args.orientation == "auto":
    ok = [c["orientation"] for c in checks if behaves(c)]
    if len(ok) != 1:
        raise SystemExit("the plate controls do not single out one orientation of the layout "
                         f"(orientations in which they behave: {ok}); not analysed")
    orientation = ok[0]
else:
    orientation = args.orientation
print("layout applied:", "as the files are named" if orientation == "name" else "turned by 180 degrees")

# ---- plate statistics (as at level 4) ------------------------------------------------------------------
tab = joined(orientation)
good = tab["qc_pass"] == 1
tab["resp_m"] = tab["marker_readout"] / vehicle_expectation(tab, "marker_readout", args.plate == "model")
tab["resp_n"] = tab["n_cells"] / vehicle_expectation(tab, "n_cells", args.plate == "model")

if args.dump:
    tab[["well", "image", "role", "compound", "dose_uM", "qc_pass", "n_cells", "marker_per_cell", "resp_m", "resp_n"]].to_csv(
        out / "responses.tsv", sep="\t", index=False, float_format="%.4g")

LE = np.linspace(-3.0, 3.5, 196)[:, None, None, None]
HH = np.linspace(0.5, 3.0, 26)[None, :, None, None]
BB = np.linspace(0.0, 0.9, 37)[None, None, :, None]
TT = np.linspace(0.6, 2.4, 37)[None, None, :, None]


def half_point(d, y):
    ok = np.isfinite(y)
    d, y = d[ok], y[ok]
    if len(d) < 4:
        return "inactive"
    if args.fit == "free":
        # textbook default: both plateaus free, potency = EC50 (midpoint between the fitted plateaus)
        best_fit = None
        for bottom in np.linspace(0.0, 2.4, 49):
            pred = bottom + (TT - bottom) / (1 + (d[None, None, None, :] / 10 ** LE) ** HH)
            sse = ((pred - y[None, None, None, :]) ** 2).sum(axis=3)
            i, j, k = np.unravel_index(int(np.argmin(sse)), sse.shape)
            if best_fit is None or sse[i, j, k] < best_fit[0]:
                best_fit = (sse[i, j, k], LE[i, 0, 0, 0], TT[0, 0, k, 0], bottom)
        _, le, top, bottom = best_fit
        x50 = 10 ** le
        return f"{x50:.4g}" if abs(top - bottom) >= 0.25 and d.min() <= x50 <= d.max() else "inactive"
    pred = BB + (1 - BB) / (1 + (d[None, None, None, :] / 10 ** LE) ** HH)
    sse = ((pred - y[None, None, None, :]) ** 2).sum(axis=3)
    i, j, k = np.unravel_index(int(np.argmin(sse)), sse.shape)
    le, h, b = LE[i, 0, 0, 0], HH[0, j, 0, 0], BB[0, 0, k, 0]
    if b >= 0.5:
        return "inactive"
    x50 = 10 ** le * ((1 - b) / (0.5 - b) - 1) ** (1 / h)
    return f"{x50:.4g}" if x50 <= d.max() else "inactive"


lines = []
for name, s in tab[(tab["role"] == "compound") & good].groupby("compound"):
    d = s["dose_uM"].to_numpy(float)
    lines.append((name, half_point(d, s["resp_m"].to_numpy(float)), half_point(d, s["resp_n"].to_numpy(float))))
pd.DataFrame(lines, columns=["compound", "ic50_marker_uM", "gi50_count_uM"]).to_csv(
    out / "dose_response.tsv", sep="\t", index=False)
print(pd.DataFrame(lines, columns=["compound", "ic50_marker_uM", "gi50_count_uM"]).to_string(index=False))

# ---- per-well tables: a row carries the name of the image it was measured on (--rows plate: the name of
# the plate well that image shows, which is the other consistent way to file them) ------------------------
key = (lambda w: w) if args.rows == "image" or orientation == "name" else turned
per_image = tab.set_index("image").loc[names].reset_index()
wm = pd.DataFrame({"well": per_image["image"].map(key), "n_cells": per_image["n_cells"].astype(int),
                   "marker_per_cell": per_image["marker_per_cell"].map(lambda v: f"{v:.2f}" if np.isfinite(v) else "NA"),
                   "qc_pass": per_image["qc_pass"].astype(int), "image": per_image["image"],
                   "plate_well": per_image["well"]})
wm.sort_values("well").to_csv(out / "well_measurements.tsv", sep="\t", index=False)
n_of = dict(zip(per_image["image"], per_image["n_cells"]))
pd.DataFrame([(key(w), m, int(n_of[w])) for w, m in sorted(manual.items())],
             columns=["field", "manual_count", "pipeline_count"]).to_csv(out / "validation.tsv", sep="\t", index=False)
