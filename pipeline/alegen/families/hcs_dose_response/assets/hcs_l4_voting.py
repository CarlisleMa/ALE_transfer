#!/usr/bin/env python3
"""Level-4 reference solution: classical image analysis without Cellpose.

  nuclei      radial-symmetry voting: outline pixels vote for the centre along their gradient, peaks of
              the vote map are nuclei. Overlapping nuclei keep most of their outline, so each still
              collects its votes; chromatin texture votes incoherently.
  fragments   SOP section 1: an object whose half-maximum region is smaller than a 9.5 um circle
              (between the largest fragment, 6.5 um, and the smallest nucleus, 11 um across) is a fragment
  frame       SOP section 1: objects on the bottom or right image edge are not counted, objects on the
              top or left edge are (and are not tested for size: they are cut nuclei)
  tuning      SOP section 2: the voting parameters are chosen on the annotated fields (smallest
              worst-case error relative to the criterion), the comparison is written to validation.tsv
  background  mode of the field histogram
  marker      mean over a ring 1-4 px outside each counted nucleus
  focus       band ratio of the power spectrum, noise-corrected
  plate model edge ratio from medians + Theil-Sen slopes on the vehicle wells
  curve fit   exhaustive grid search of the SOP model (upper asymptote fixed at 1)

The second solution (hcs_l4_cellpose.py) shares no code with this file and uses a different method at
every step (Cellpose nuclei model with settings chosen on the annotated fields, mask area, median
background far from nuclei, 0-4 px ring, edge-gradient focus metric, least-squares plate model,
scipy least_squares fit).

The switches reproduce what a careless analyst would deliver; they exist so that the task's verifier
can be attacked with real outputs built on a correct base (see the family's battery).
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi
from skimage.feature import peak_local_max
from skimage.segmentation import expand_labels, watershed

ap = argparse.ArgumentParser()
ap.add_argument("--counting", choices=["sop", "fragments", "noframe"], default="sop")   # fragments = fragments counted
ap.add_argument("--qc", choices=["focus", "none"], default="focus")
ap.add_argument("--readout", choices=["percell", "wholewell"], default="percell")
ap.add_argument("--background", choices=["local", "global", "none"], default="local")
ap.add_argument("--plate", choices=["model", "none"], default="model")
ap.add_argument("--fit", choices=["sop", "free"], default="sop")
args = ap.parse_args()

inp, out = Path("/task/input"), Path("/task/output")
out.mkdir(parents=True, exist_ok=True)
lay = pd.read_csv(inp / "plate_layout.tsv", sep="\t")
lay["r"] = lay["well"].str[0].map(lambda s: ord(s) - 65)
lay["c"] = lay["well"].str[1:].astype(int) - 1
lay["edge"] = lay["r"].isin([0, lay["r"].max()]) | lay["c"].isin([0, lay["c"].max()])
manual = {p.stem: sum(1 for ln in p.read_text().splitlines()[1:] if ln.strip())
          for p in sorted((inp / "annotations").glob("*.tsv"))}

dark = tifffile.imread(inp / "flatfield" / "dark.tif").astype(np.float64)
gainmap = {}
for ch in (1, 2):
    f = tifffile.imread(inp / "flatfield" / f"flat_ch{ch}.tif").astype(np.float64) - dark
    f = ndi.uniform_filter(f, 9)
    gainmap[ch] = f / np.median(f)


def read(well, ch):
    return (tifffile.imread(inp / "images" / f"{well}_ch{ch}.tif").astype(np.float64) - dark) / gainmap[ch]


# ---- nuclei --------------------------------------------------------------------------------------
PX_UM = 1.3
FRAGMENT_AREA = np.pi / 4 * (9.5 / PX_UM) ** 2      # px^2 at half maximum; SOP: fragments <= 6.5 um, nuclei >= 11 um across
TOL_ABS, TOL_REL = 3, 0.10                          # SOP section 2


def nuclear_level(im):
    sm = ndi.gaussian_filter(im, 1.0)
    b = np.median(sm)
    return float(np.median(sm[sm > b + 0.5 * (np.percentile(sm, 99) - b)]) - b)


nuc_img = {w: read(w, 1) for w in lay["well"]}
LEVEL = float(np.median([nuclear_level(nuc_img[w]) for w in lay.loc[lay["role"] == "vehicle", "well"]]))
_prep = {}


def prepared(w):
    """Smoothed image, background, foreground, gradient and the size of the nuclei of this field relative
    to the usual 5.5 px half-width (treated nuclei can be larger; all radii scale with it)."""
    if w not in _prep:
        sm = ndi.gaussian_filter(nuc_img[w], 1.0)
        bg = float(np.median(sm))
        fg = sm > bg + 0.25 * LEVEL
        comp, n = ndi.label(fg)
        scale = 1.0
        if n:
            idx = np.arange(1, n + 1)
            big = ndi.sum(fg, comp, idx) >= 60
            if big.sum() >= 3:
                half_width = np.median(ndi.maximum(ndi.distance_transform_edt(fg), comp, idx)[big])
                scale = float(np.clip((half_width - 0.8) / 5.5, 0.85, 2.0))
        gy, gx = np.gradient(sm)
        _prep[w] = (sm, bg, fg, scale, gy, gx)
    return _prep[w]


def segment(w, edge_thr, vote_sigma):
    """Returns (labels, ids of counted nuclei, ids of all nuclei incl. those outside the counting frame)."""
    sm, bg, fg, scale, gy, gx = prepared(w)
    H, W = sm.shape
    mag = np.hypot(gy, gx)
    local = np.maximum(ndi.maximum_filter(sm, 7) - bg, 0.25 * LEVEL)
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
for w in lay["well"]:
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
tab = lay.merge(pd.DataFrame(rec, columns=["well", "n_cells", "raw_ring", "bg", "field_mean", "sharp"]), on="well")
tab["qc_pass"] = np.where(tab["sharp"] < 0.35 * np.nanmedian(tab["sharp"]), 0, 1) if args.qc == "focus" else 1
good = tab["qc_pass"] == 1
if args.background == "local":
    bgv = tab["bg"]
elif args.background == "global":
    bgv = float(np.median(tab.loc[(tab["role"] == "vehicle") & good, "bg"]))
else:
    bgv = 0.0
tab["marker_per_cell"] = tab["raw_ring"] - bgv
tab["marker_readout"] = tab["marker_per_cell"] if args.readout == "percell" else tab["field_mean"] - bgv


def theil_sen(x, y):
    i, j = np.triu_indices(len(x), 1)
    dx = x[j] - x[i]
    keep = dx != 0
    return float(np.median((y[j] - y[i])[keep] / dx[keep]))


def vehicle_expectation(col):
    v = tab[(tab["role"] == "vehicle") & good]
    if args.plate == "none":
        return np.full(len(tab), float(v[col].mean()))
    y = np.log(v[col].to_numpy(dtype=float))
    c, r, e = v["c"].to_numpy(float), v["r"].to_numpy(float), v["edge"].to_numpy(bool)
    bc = br = le = 0.0
    for _ in range(6):
        res = y - bc * c - br * r
        le = np.median(res[e]) - np.median(res[~e])
        z = y - le * e
        bc = theil_sen(c, z - br * r)
        br = theil_sen(r, z - bc * c)
    a = np.median(y - le * e - bc * c - br * r)
    return np.exp(a + le * tab["edge"].to_numpy(float) + bc * tab["c"].to_numpy(float) + br * tab["r"].to_numpy(float))


tab["resp_m"] = tab["marker_readout"] / vehicle_expectation("marker_readout")
tab["resp_n"] = tab["n_cells"] / vehicle_expectation("n_cells")

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
wm = tab[["well", "n_cells", "marker_per_cell", "qc_pass"]].copy()
wm["marker_per_cell"] = wm["marker_per_cell"].map(lambda v: f"{v:.2f}" if np.isfinite(v) else "NA")
wm.to_csv(out / "well_measurements.tsv", sep="\t", index=False)
n_of = dict(zip(tab["well"], tab["n_cells"]))
pd.DataFrame([(w, m, int(n_of[w])) for w, m in sorted(manual.items())],
             columns=["field", "manual_count", "pipeline_count"]).to_csv(out / "validation.tsv", sep="\t", index=False)
