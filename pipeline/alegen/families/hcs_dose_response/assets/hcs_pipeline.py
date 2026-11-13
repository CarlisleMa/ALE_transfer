#!/usr/bin/env python3
"""Reference pipeline: 2-channel plate images -> per-well measurements -> per-compound potencies.

Default = the oracle (Cellpose `nuclei` model for segmentation). The switches reproduce what a
careless analyst would deliver; they exist so that the task's verifier can be attacked with real
outputs (see the family's battery), not because they are sensible options.

  --segmenter cellpose|watershed|threshold   threshold = global Otsu + connected components
  --plate     model|none                     none = normalise to the plate mean of the vehicle wells
  --readout   percell|wholewell              wholewell = mean marker intensity of the field
  --background local|global|none             global = one plate-wide constant, none = camera offset only
  --qc        focus|none                     none = every well passes
  --cellpose-norm plate|image                image = Cellpose's default per-image percentile scaling
"""
import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi
from scipy.optimize import least_squares
from skimage import feature, filters, measure, morphology, segmentation

ap = argparse.ArgumentParser()
ap.add_argument("--segmenter", choices=["cellpose", "watershed", "threshold"], default="cellpose")
ap.add_argument("--plate", choices=["model", "none"], default="model")
ap.add_argument("--readout", choices=["percell", "wholewell"], default="percell")
ap.add_argument("--background", choices=["local", "global", "none"], default="local")
ap.add_argument("--qc", choices=["focus", "none"], default="focus")
ap.add_argument("--cellpose-norm", choices=["plate", "image"], default="plate")
ap.add_argument("--input", default="/task/input")
ap.add_argument("--output", default="/task/output")
ap.add_argument("--threads", type=int, default=int(os.environ.get("HCS_THREADS", "2")))
args = ap.parse_args()
inp, out = Path(args.input), Path(args.output)
out.mkdir(parents=True, exist_ok=True)

layout = pd.read_csv(inp / "plate_layout.tsv", sep="\t")
layout["row"] = [ord(w[0]) - ord("A") for w in layout["well"]]
layout["col"] = [int(w[1:]) - 1 for w in layout["well"]]
nrow, ncol = layout["row"].max() + 1, layout["col"].max() + 1
layout["edge"] = ((layout["row"] == 0) | (layout["row"] == nrow - 1) |
                  (layout["col"] == 0) | (layout["col"] == ncol - 1)).astype(float)

# ---- 1. illumination correction ---------------------------------------------------------------
dark = tifffile.imread(inp / "flatfield" / "dark.tif").astype(float)
flat = {}
for ch in (1, 2):
    f = ndi.gaussian_filter(tifffile.imread(inp / "flatfield" / f"flat_ch{ch}.tif").astype(float) - dark, 4)
    flat[ch] = f / f.mean()


def load(well, ch):
    return (tifffile.imread(inp / "images" / f"{well}_ch{ch}.tif").astype(float) - dark) / flat[ch]


nuc = {w: load(w, 1) for w in layout["well"]}

# ---- 2. nuclei ---------------------------------------------------------------------------------
def seg_classical(img, split):
    sm = ndi.gaussian_filter(img, 1.0)
    base = np.median(sm)
    noise = 1.4826 * np.median(np.abs(sm - base))
    thr = max(filters.threshold_otsu(sm), base + 8 * noise)
    bw = morphology.remove_small_objects(sm > thr, 20)
    if not split:
        return measure.label(bw)
    dist = ndi.gaussian_filter(ndi.distance_transform_edt(bw), 0.8)
    comp = measure.label(bw)
    pk = feature.peak_local_max(dist, min_distance=3, labels=comp, exclude_border=False)
    markers = np.zeros(bw.shape, dtype=int)
    markers[tuple(pk.T)] = np.arange(1, len(pk) + 1)
    return segmentation.watershed(-dist, markers, mask=bw)


if args.segmenter == "cellpose":
    import torch
    from cellpose import models
    torch.set_num_threads(args.threads)
    wells = list(layout["well"])
    model = models.CellposeModel(gpu=False, model_type="nuclei")
    if args.cellpose_norm == "plate":
        # one intensity scale for the whole plate: per-image percentile scaling (the default) blows up
        # fields that hold only a handful of nuclei, and the model then drops most of them
        lo = float(np.median([np.percentile(nuc[w], 1) for w in wells]))
        hi = float(np.median([np.percentile(nuc[w], 99) for w in layout.loc[layout["role"] == "vehicle", "well"]]))
        imgs, norm = [np.clip((nuc[w] - lo) / (hi - lo), 0, None) for w in wells], False
    else:
        imgs, norm = [nuc[w] for w in wells], True
    # nuclei are 10-15 px across, inside the working range of the model at native resolution;
    # diameter=17 (the training scale) switches the rescaling off, which is 5x faster on CPU
    res = model.eval(imgs, diameter=17.0, channels=[0, 0], batch_size=8, normalize=norm)[0]
    labels = {w: morphology.remove_small_objects(np.asarray(m).astype(int), 20) for w, m in zip(wells, res)}
else:
    labels = {w: seg_classical(nuc[w], args.segmenter == "watershed") for w in layout["well"]}

# ---- 3. per-well measurements ------------------------------------------------------------------
rows = []
for w in layout["well"]:
    lab = labels[w]
    ids = np.unique(lab[lab > 0])
    n = len(ids)
    m2 = load(w, 2)
    grown = segmentation.expand_labels(lab, 4)
    ring = np.where(lab > 0, 0, grown)
    far = ndi.distance_transform_edt(lab == 0) > 12
    bg_local = float(np.median(m2[far])) if far.sum() > 400 else float(np.percentile(m2, 5))
    per_cell = np.array([])
    if n:
        area = ndi.sum(np.ones_like(m2), ring, ids)
        per_cell = (ndi.sum(m2, ring, ids)[area > 0]) / area[area > 0]   # enclosed nuclei have no ring
    # focus: sharpest nuclear edges relative to the nuclear signal of the same field
    g = ndi.gaussian_gradient_magnitude(nuc[w], 1.0)
    amp = np.percentile(nuc[w], 99.9) - np.median(nuc[w])
    noise = 1.4826 * np.median(np.abs(nuc[w] - np.median(nuc[w])))
    focus = float(np.percentile(g, 99.9) / amp) if amp > 12 * noise else np.nan
    rows.append({"well": w, "n_cells": n, "raw_ring": float(np.mean(per_cell)) if len(per_cell) else np.nan,
                 "bg_local": bg_local, "field_mean": float(m2.mean()), "focus": focus})
meas = layout.merge(pd.DataFrame(rows), on="well")

# ---- 4. focus QC -------------------------------------------------------------------------------
if args.qc == "focus":
    ref = np.nanmedian(meas["focus"])
    meas["qc_pass"] = (~(meas["focus"] < 0.6 * ref)).astype(int)
else:
    meas["qc_pass"] = 1
ok = meas["qc_pass"] == 1
vehicle = (meas["role"] == "vehicle") & ok

# ---- 5. background and readouts ----------------------------------------------------------------
if args.background == "local":
    bg = meas["bg_local"]
elif args.background == "global":
    bg = float(np.median(meas.loc[vehicle, "bg_local"]))
else:
    bg = 0.0
meas["marker_per_cell"] = meas["raw_ring"] - bg
meas["marker_readout"] = meas["marker_per_cell"] if args.readout == "percell" else meas["field_mean"] - bg

# ---- 6. plate-position model from the vehicle wells --------------------------------------------
def expected_vehicle(col):
    v = meas.loc[vehicle]
    if args.plate == "none":
        return np.full(len(meas), float(v[col].mean()))
    X = np.column_stack([np.ones(len(v)), v["col"], v["row"], v["edge"]])
    beta, *_ = np.linalg.lstsq(X, np.log(v[col].astype(float)), rcond=None)
    Xa = np.column_stack([np.ones(len(meas)), meas["col"], meas["row"], meas["edge"]])
    return np.exp(Xa @ beta)


meas["resp_marker"] = meas["marker_readout"] / expected_vehicle("marker_readout")
meas["resp_count"] = meas["n_cells"] / expected_vehicle("n_cells")

# ---- 7. dose-response --------------------------------------------------------------------------
def potency(d, y):
    d, y = np.asarray(d, dtype=float), np.asarray(y, dtype=float)
    keep = np.isfinite(y)
    d, y = d[keep], y[keep]
    if len(d) < 4:
        return "inactive"

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

wm = meas[["well", "n_cells", "marker_per_cell", "qc_pass"]].copy()
wm["marker_per_cell"] = [f"{v:.2f}" if np.isfinite(v) else "NA" for v in wm["marker_per_cell"]]
wm.to_csv(out / "well_measurements.tsv", sep="\t", index=False)
meas.to_csv(out / "qc_table.tsv", sep="\t", index=False, float_format="%.5g")
