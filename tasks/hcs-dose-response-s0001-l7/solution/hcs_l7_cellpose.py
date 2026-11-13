#!/usr/bin/env python3
"""Level 7, independent second solution: Cellpose with settings validated on the annotated fields (the
level-4 second solution), plus a check of which way round the plate layout fits the images.
No code shared with hcs_l7_voting.py (the reference solution), a different method at every step; the
orientation of the plate is found from the dose series, without using the control wells.

It is run as a sequence of stages (see solve.sh); intermediate results live in --work.

  prepare                      illumination correction of channel 1, one plate-wide intensity scale
  try --model M --diameter D --flow F
                               count the annotated fields with one Cellpose setting under the SOP's
                               counting rules and append the comparison with the manual counts
  choose                       pick the setting that meets the SOP criterion on every annotated field
                               (smallest worst-case error); stop if none does
  segment                      segment all wells with the chosen setting
  measure                      marker per cell, field background, focus metric
  orient                       the layout is laid over the images as the files are named and turned by 180
                               degrees; for both, the rank correlation of each dose series with the
                               concentration is computed on readouts that were normalised without the
                               control wells (marker: robust trend along the image columns over all fields;
                               cell number: median of the outer ring / of the inner fields). A treatment can
                               lower a readout with rising concentration; series that rise instead are what a
                               turned plate produces. The orientation with the larger balance of falling over
                               rising series is used; if the balance is not clear the run stops.
  analyse                      focus QC, plate-position model, dose-response fits, deliverables, with the
                               layout in the orientation found; per-well rows carry the name of the image

Counting rules of the SOP as implemented here: an object smaller than a 9.5 um circle (between the
largest fragment, 6.5 um, and the smallest nucleus, 11 um) is a fragment unless it lies on the top or
left image edge (a cut nucleus); an object on the bottom or right edge is not counted.

--force (go on with a setting that fails the criterion) and --counting all (every mask >= 20 px) exist
so that the outputs of a careless analyst can be produced. They are not sensible options.
"""
import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi
from scipy.optimize import least_squares
from skimage import segmentation

ap = argparse.ArgumentParser()
ap.add_argument("stage", choices=["prepare", "try", "choose", "segment", "measure", "orient", "analyse"])
ap.add_argument("--model", default="nuclei")
ap.add_argument("--diameter", type=float, default=17.0)
ap.add_argument("--flow", type=float, default=0.4)
ap.add_argument("--counting", choices=["sop", "all"], default="sop")          # all = every mask >= 20 px
ap.add_argument("--fit", choices=["sop", "free"], default="sop")              # free = 4PL with free plateaus, EC50
ap.add_argument("--force", action="store_true")       # choose: go on with the best setting although it fails
ap.add_argument("--input", default="/task/input")
ap.add_argument("--output", default="/task/output")
ap.add_argument("--work", default="/tmp/hcs_l4_work")
ap.add_argument("--threads", type=int, default=int(os.environ.get("HCS_THREADS", "2")))
args = ap.parse_args()
inp, out, work = Path(args.input), Path(args.output), Path(args.work)
work.mkdir(parents=True, exist_ok=True)

PX_UM = 1.3                                            # assay notes: 1.3 um per pixel
MIN_AREA = np.pi / 4 * (9.5 / PX_UM) ** 2              # px^2, see module docstring
TOL_ABS, TOL_REL = 3, 0.10                             # SOP section 2

layout = pd.read_csv(inp / "plate_layout.tsv", sep="\t")
layout["row"] = [ord(w[0]) - ord("A") for w in layout["well"]]
layout["col"] = [int(w[1:]) - 1 for w in layout["well"]]
nrow, ncol = layout["row"].max() + 1, layout["col"].max() + 1
layout["edge"] = ((layout["row"] == 0) | (layout["row"] == nrow - 1) |
                  (layout["col"] == 0) | (layout["col"] == ncol - 1)).astype(float)
wells = list(layout["well"])                           # the image files carry the well names of the layout


def opposite(well):
    """Name of the well at the opposite position of the plate (the plate turned by 180 degrees)."""
    return f"{chr(ord('A') + nrow - 1 - (ord(well[0]) - ord('A')))}{ncol + 1 - int(well[1:]):02d}"


def image_table():
    """Per-image measurements with the focus QC (a property of the image, not of the layout)."""
    m = pd.read_csv(work / "meas.tsv", sep="\t").rename(columns={"well": "image"})
    m["qc_pass"] = (~(m["focus"] < 0.6 * np.nanmedian(m["focus"]))).astype(int)
    m["marker_per_cell"] = m["raw_ring"] - m["bg_local"]
    return m


def laid_over(m, orientation):
    """Layout rows (plate wells) with the measurements of the image that shows each of them."""
    t = layout.copy()
    t["image"] = t["well"] if orientation == "name" else t["well"].map(opposite)
    return t.merge(m, on="image", validate="one_to_one")
manual = {p.stem: len(pd.read_csv(p, sep="\t")) for p in sorted((inp / "annotations").glob("*.tsv"))}

dark = tifffile.imread(inp / "flatfield" / "dark.tif").astype(float)


def flat_of(ch):
    f = ndi.gaussian_filter(tifffile.imread(inp / "flatfield" / f"flat_ch{ch}.tif").astype(float) - dark, 4)
    return f / f.mean()


def load(well, ch, flat):
    return (tifffile.imread(inp / "images" / f"{well}_ch{ch}.tif").astype(float) - dark) / flat


def counted_ids(lab):
    """Labels the SOP counts in one field (section 1: fragments, counting frame)."""
    ids, area = np.unique(lab[lab > 0], return_counts=True)
    if args.counting == "all":
        return [int(i) for i, a in zip(ids, area) if a >= 20]
    low = set(np.unique(np.concatenate([lab[0, :], lab[:, 0]]))) - {0}       # top / left edge
    high = set(np.unique(np.concatenate([lab[-1, :], lab[:, -1]]))) - {0}    # bottom / right edge
    return [int(i) for i, a in zip(ids, area) if i not in high and (a >= MIN_AREA or i in low)]


def cellpose_labels(names, model_type, diameter, flow):
    import torch
    from cellpose import models
    torch.set_num_threads(args.threads)
    prep = np.load(work / "prep.npz")
    lo, hi = float(prep["lo"]), float(prep["hi"])
    imgs = [np.clip((prep[w].astype(float) - lo) / (hi - lo), 0, None) for w in names]
    model = models.CellposeModel(gpu=False, model_type=model_type)
    res = model.eval(imgs, diameter=diameter, channels=[0, 0], batch_size=8, normalize=False,
                     flow_threshold=flow, cellprob_threshold=0.0)[0]
    return {w: np.asarray(m).astype(np.int32) for w, m in zip(names, res)}


# ------------------------------------------------------------------------------------------------
if args.stage == "prepare":
    flat1 = flat_of(1)
    nuc = {w: load(w, 1, flat1) for w in wells}
    # one intensity scale for the plate, from all fields (the layout has not been checked at this point)
    lo = float(np.median([np.percentile(nuc[w], 1) for w in wells]))
    hi = float(np.median([np.percentile(nuc[w], 99) for w in wells]))
    np.savez_compressed(work / "prep.npz", lo=lo, hi=hi, **{w: nuc[w].astype(np.float32) for w in wells})
    (work / "calibration.tsv").unlink(missing_ok=True)          # a new run starts a new comparison table
    print(f"{len(wells)} wells, {len(manual)} annotated fields, intensity scale {lo:.0f}..{hi:.0f}")

elif args.stage == "try":
    labs = cellpose_labels(sorted(manual), args.model, args.diameter, args.flow)
    rows = []
    for w in sorted(manual):
        n = len(counted_ids(labs[w]))
        rows.append({"model": args.model, "diameter": args.diameter, "flow": args.flow, "field": w,
                     "manual_count": manual[w], "pipeline_count": n,
                     "within": int(abs(n - manual[w]) <= max(TOL_ABS, TOL_REL * manual[w]))})
    tab = pd.DataFrame(rows)
    f = work / "calibration.tsv"
    tab.to_csv(f, sep="\t", index=False, mode="a", header=not f.exists())
    print(tab.to_string(index=False))

elif args.stage == "choose":
    cal = pd.read_csv(work / "calibration.tsv", sep="\t")
    cal["excess"] = (cal["pipeline_count"] - cal["manual_count"]).abs() / np.maximum(TOL_ABS, TOL_REL * cal["manual_count"])
    summ = cal.groupby(["model", "diameter", "flow"], as_index=False).agg(worst=("excess", "max"), ok=("within", "min"))
    print(summ.to_string(index=False))
    good = summ[summ["ok"] == 1].sort_values("worst")
    if good.empty and not args.force:
        raise SystemExit("no setting meets the SOP criterion on the annotated fields")
    best = (good if not good.empty else summ.sort_values("worst")).iloc[0]
    (work / "choice.json").write_text(json.dumps({"model": best["model"], "diameter": float(best["diameter"]),
                                                  "flow": float(best["flow"])}))
    print("chosen:", best["model"], best["diameter"], best["flow"])

elif args.stage == "segment":
    c = json.loads((work / "choice.json").read_text())
    labs = cellpose_labels(wells, c["model"], c["diameter"], c["flow"])
    np.savez_compressed(work / "labels.npz", **labs)
    print("segmented", len(labs), "wells with", c)

elif args.stage == "measure":
    labs = np.load(work / "labels.npz")
    prep = np.load(work / "prep.npz")
    flat2 = flat_of(2)
    rows = []
    for w in wells:
        lab = labs[w]
        ids = counted_ids(lab)
        cells = np.where(np.isin(lab, ids), lab, 0)
        m2 = load(w, 2, flat2)
        ring = np.where(lab > 0, 0, segmentation.expand_labels(cells, 4))
        far = ndi.distance_transform_edt(lab == 0) > 12
        bg_local = float(np.median(m2[far])) if far.sum() > 400 else float(np.percentile(m2, 5))
        per_cell = np.array([])
        if ids:
            area = ndi.sum(np.ones_like(m2), ring, ids)
            per_cell = (ndi.sum(m2, ring, ids)[area > 0]) / area[area > 0]   # enclosed nuclei have no ring
        # focus: sharpest edges of channel 1 relative to the signal of the same field
        n1 = prep[w].astype(float)
        g = ndi.gaussian_gradient_magnitude(n1, 1.0)
        amp = np.percentile(n1, 99.9) - np.median(n1)
        noise = 1.4826 * np.median(np.abs(n1 - np.median(n1)))
        focus = float(np.percentile(g, 99.9) / amp) if amp > 12 * noise else np.nan
        rows.append({"well": w, "n_cells": len(ids), "raw_ring": float(np.mean(per_cell)) if len(per_cell) else np.nan,
                     "bg_local": bg_local, "focus": focus})
    pd.DataFrame(rows).to_csv(work / "meas.tsv", sep="\t", index=False)
    print(pd.DataFrame(rows).describe().to_string())

elif args.stage == "orient":
    from scipy.stats import spearmanr
    m = image_table()
    m["irow"] = [ord(w[0]) - ord("A") for w in m["image"]]
    m["icol"] = [int(w[1:]) - 1 for w in m["image"]]
    m["ring"] = (m["irow"] == 0) | (m["irow"] == nrow - 1) | (m["icol"] == 0) | (m["icol"] == ncol - 1)
    use = (m["qc_pass"] == 1) & (m["marker_per_cell"] > 0)
    # marker: the staining level changes along the columns. Most fields are not inhibited, so a line through
    # the dense band of log(marker) against the image column gives the untreated level without naming a control.
    x, y = m.loc[use, "icol"].to_numpy(float), np.log(m.loc[use, "marker_per_cell"].to_numpy(float))
    top = np.array([np.percentile(y[x == c], 62.5) for c in range(ncol)])
    slope, icpt = np.polyfit(np.arange(ncol), top, 1)
    for _ in range(12):
        keep = np.abs(y - (icpt + slope * x)) < 0.25
        slope, icpt = np.polyfit(x[keep], y[keep], 1)
    m["marker_rel"] = m["marker_per_cell"] / np.exp(icpt + slope * m["icol"])
    # cell number: outer ring and inner fields have their own level; both sets are the same either way round
    level = m[m["qc_pass"] == 1].groupby("ring")["n_cells"].median()
    m["count_rel"] = m["n_cells"] / m["ring"].map(level)
    rows = []
    for orientation in ("name", "turned"):
        t = laid_over(m, orientation)
        s = t[(t["role"] == "compound") & (t["qc_pass"] == 1)]
        for name, g in s.groupby("compound"):
            for readout in ("marker_rel", "count_rel"):
                gg = g[np.isfinite(g[readout])]
                rho = spearmanr(gg["dose_uM"], gg[readout]).correlation if len(gg) >= 6 else np.nan
                rows.append({"orientation": orientation, "compound": name, "readout": readout, "wells": len(gg), "rho": rho})
    tabr = pd.DataFrame(rows)
    tabr.to_csv(work / "series_direction.tsv", sep="\t", index=False, float_format="%.3f")
    summary = {}
    for orientation, g in tabr.groupby("orientation"):
        summary[orientation] = {"falling": int((g["rho"] <= -0.7).sum()), "rising": int((g["rho"] >= 0.7).sum())}
        summary[orientation]["balance"] = summary[orientation]["falling"] - summary[orientation]["rising"]
    print(json.dumps(summary))
    best, other = sorted(summary, key=lambda o: -summary[o]["balance"])
    if summary[best]["balance"] < 4 or summary[best]["balance"] - summary[other]["balance"] < 4:
        raise SystemExit("the dose series do not single out one orientation of the layout; not analysed")
    (work / "orientation.json").write_text(json.dumps({"orientation": best, "series": summary}))
    print("layout applied:", "as the files are named" if best == "name" else "turned by 180 degrees")

elif args.stage == "analyse":
    orientation = json.loads((work / "orientation.json").read_text())["orientation"]
    meas = laid_over(image_table(), orientation)        # one row per plate well, with the image that shows it
    out.mkdir(parents=True, exist_ok=True)
    ok = meas["qc_pass"] == 1
    vehicle = (meas["role"] == "vehicle") & ok

    def expected_vehicle(col):
        v = meas.loc[vehicle]
        X = np.column_stack([np.ones(len(v)), v["col"], v["row"], v["edge"]])
        beta, *_ = np.linalg.lstsq(X, np.log(v[col].astype(float)), rcond=None)
        Xa = np.column_stack([np.ones(len(meas)), meas["col"], meas["row"], meas["edge"]])
        return np.exp(Xa @ beta)

    meas["resp_marker"] = meas["marker_per_cell"] / expected_vehicle("marker_per_cell")
    meas["resp_count"] = meas["n_cells"] / expected_vehicle("n_cells")

    def potency(d, y):
        d, y = np.asarray(d, dtype=float), np.asarray(y, dtype=float)
        keep = np.isfinite(y)
        d, y = d[keep], y[keep]
        if len(d) < 4:
            return "inactive"
        if args.fit == "free":
            # the textbook default: four free parameters, potency = fitted EC50 (midpoint between the
            # fitted plateaus) when it lies inside the tested range and the curve spans >= 0.25
            best = None
            for le0 in np.linspace(-2.0, 2.0, 5):
                r = least_squares(lambda p: p[3] + (p[2] - p[3]) / (1 + (d / 10 ** p[0]) ** p[1]) - y,
                                  [le0, 1.0, max(y.max(), 0.1), y.min()],
                                  bounds=([-4.0, 0.4, 0.0, -0.5], [4.0, 4.0, 5.0, 5.0]))
                if best is None or r.cost < best.cost:
                    best = r
            le, h, top, bottom = best.x
            ec50 = 10 ** le
            return f"{ec50:.4g}" if abs(top - bottom) >= 0.25 and d.min() <= ec50 <= d.max() else "inactive"

        def model(p):
            return p[2] + (1 - p[2]) / (1 + (d / 10 ** p[0]) ** p[1])

        best = None
        for le0 in np.linspace(-2.0, 2.0, 5):
            for b0 in (0.05, 0.6):
                r = least_squares(lambda p: model(p) - y, [le0, 1.0, b0],
                                  bounds=([-4.0, 0.4, 0.0], [4.0, 4.0, 0.95]))
                if best is None or r.cost < best.cost:
                    best = r
        le, h, b = best.x
        if b >= 0.5:
            return "inactive"
        ic50 = 10 ** le * ((1 - b) / (0.5 - b) - 1) ** (1 / h)
        return f"{ic50:.4g}" if ic50 <= d.max() else "inactive"

    res = []
    cmp_wells = meas[(meas["role"] == "compound") & ok]
    for name in sorted(set(layout.loc[layout["role"] == "compound", "compound"])):
        s = cmp_wells[cmp_wells["compound"] == name].sort_values("dose_uM")
        res.append({"compound": name, "ic50_marker_uM": potency(s["dose_uM"], s["resp_marker"]),
                    "gi50_count_uM": potency(s["dose_uM"], s["resp_count"])})
    pd.DataFrame(res).to_csv(out / "dose_response.tsv", sep="\t", index=False)

    # per-well tables: a row carries the name of the image (and annotation file) it was measured on
    wm = meas[["image", "n_cells", "marker_per_cell", "qc_pass", "well"]].rename(
        columns={"image": "well", "well": "plate_well"}).sort_values("well")
    wm["marker_per_cell"] = [f"{v:.2f}" if np.isfinite(v) else "NA" for v in wm["marker_per_cell"]]
    wm.to_csv(out / "well_measurements.tsv", sep="\t", index=False)
    n_of = dict(zip(meas["image"], meas["n_cells"]))
    val = pd.DataFrame([{"field": w, "manual_count": manual[w], "pipeline_count": int(n_of[w])} for w in sorted(manual)])
    val.to_csv(out / "validation.tsv", sep="\t", index=False)
    meas.to_csv(out / "qc_table.tsv", sep="\t", index=False, float_format="%.5g")
    print(val.to_string(index=False))
    print(pd.DataFrame(res).to_string(index=False))
