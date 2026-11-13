#!/usr/bin/env python3
"""Reference analysis of the 8-colour PBMC panel: spillover from the single-stain controls
(median of the positive minus median of the negative population), compensation, then a
hierarchical gate per sample with thresholds placed at density valleys of that sample.

The non-default flags replay plausible-but-wrong analyses for the validator's attack suite; the
oracle (solve.sh) uses the defaults.
"""
import argparse
from pathlib import Path

import flowio
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks

ap = argparse.ArgumentParser()
ap.add_argument("--input", default="/task/input")
ap.add_argument("--output", default="/task/output")
ap.add_argument("--comp", choices=["controls", "none", "embedded"], default="controls")
ap.add_argument("--no-doublet-gate", action="store_true")
ap.add_argument("--no-viability-gate", action="store_true")
ap.add_argument("--global-thresholds", action="store_true")
a = ap.parse_args()
inp, out = Path(a.input), Path(a.output)
out.mkdir(parents=True, exist_ok=True)

panel = pd.read_csv(inp / "panel.tsv", sep="\t")
samples = pd.read_csv(inp / "samples.tsv", sep="\t")
controls = pd.read_csv(inp / "controls.tsv", sep="\t", keep_default_na=False)
dets = list(panel["detector"])
fluors = list(panel["fluorochrome"])
det_of = dict(zip(panel["marker"], panel["detector"]))
COF = 150.0

POPS = [("cells", "all_events"), ("singlets", "cells"), ("live", "singlets"), ("lymphocytes", "live"),
        ("monocytes", "live"), ("T_cells", "lymphocytes"), ("CD4_T", "T_cells"), ("CD8_T", "T_cells"),
        ("DP_T", "T_cells"), ("DN_T", "T_cells"), ("B_cells", "lymphocytes"), ("NK_cells", "lymphocytes")]


def read_fcs(path):
    fd = flowio.FlowData(str(path))
    arr = np.asarray(fd.as_array(preprocess=False), dtype=float)
    return pd.DataFrame(arr, columns=list(fd.pnn_labels)), fd.text


def valley(y, bw, min_sep, min_prom=0.001):
    """Threshold at the density minimum between the dominant mode and the most prominent other
    mode at least `min_sep` away. None if the distribution has a single mode."""
    y = np.asarray(y, dtype=float)
    y = y[np.isfinite(y)]
    if len(y) < 100:
        return None
    lo, hi = np.percentile(y, [0.05, 99.95])
    pad = 0.05 * (hi - lo) + 1e-9
    lo, hi = lo - pad, hi + pad
    bins = 400
    h, edges = np.histogram(y, bins=bins, range=(lo, hi), density=True)
    w = (hi - lo) / bins
    d = gaussian_filter1d(h, max(bw / w, 1.0), mode="constant")
    peaks, _ = find_peaks(d, prominence=min_prom * d.max())
    if len(peaks) < 2:
        return None
    main = peaks[np.argmax(d[peaks])]
    best = None
    for p in peaks:
        if p == main or abs(p - main) * w < min_sep:
            continue
        i, j = sorted((int(p), int(main)))
        v = i + int(np.argmin(d[i:j + 1]))
        depth = d[p] - d[v]
        if best is None or depth > best[0]:
            best = (depth, v)
    if best is None:
        return None
    return float(edges[best[1]] + w / 2)


def cofactor(x):
    """arcsinh cofactor = 3 x the spread of the below-zero part of the data (compensation spreads
    unstained events symmetrically around zero), never below 150."""
    neg = x[x < 0]
    if len(neg) < max(20, 0.005 * len(x)):
        return COF
    return max(COF, 3.0 * float(np.sqrt(np.mean(neg ** 2))))


def fl_threshold(x):
    """Linear-scale threshold between the negative and the positive population of x (or inf)."""
    x = np.asarray(x, dtype=float)
    c = cofactor(x)
    t = valley(np.arcsinh(x / c), 0.12, 1.0)
    return np.inf if t is None else float(c * np.sinh(t))


# ----------------------------------------------------------------------------- spillover matrix
def spillover_from_controls():
    S = np.eye(len(dets))
    for r in controls.itertuples():
        if r.control_type != "single_stain":
            continue
        df, _ = read_fcs(inp / "fcs" / r.fcs_file)
        t = valley(np.log10(df["FSC-A"].clip(lower=1)), 0.04, 0.3)     # drop debris of cell-based controls
        if t is not None:
            df = df[np.log10(df["FSC-A"].clip(lower=1)) > t]
        i = dets.index(r.detector)
        y = df[r.detector].to_numpy()
        t = fl_threshold(y)
        if not np.isfinite(t):
            raise SystemExit(f"control {r.fcs_file}: no positive population found")
        D = df[dets].to_numpy()
        mp, mn = np.median(D[y > t], axis=0), np.median(D[y <= t], axis=0)
        S[i] = (mp - mn) / (mp[i] - mn[i])
        S[i, i] = 1.0
    return S


def embedded_matrix(text):
    key = next((k for k in text if k.lower().lstrip("$") in ("spillover", "spill")), None)
    if key is None:
        return np.eye(len(dets))
    parts = text[key].split(",")
    n = int(parts[0])
    names = parts[1:1 + n]
    M = np.array(parts[1 + n:1 + n + n * n], dtype=float).reshape(n, n)
    idx = [names.index(d) for d in dets]
    return M[np.ix_(idx, idx)]


S_ctrl = spillover_from_controls()


# ----------------------------------------------------------------------------- gating
def gate(df, S, thr):
    """Hierarchical gate. `thr` is a dict of thresholds; missing entries are derived from this
    sample and stored (so that a filled dict can be re-used as a fixed template)."""
    C = pd.DataFrame(df[dets].to_numpy() @ np.linalg.inv(S), columns=dets)

    def marker(name, mask):
        """Events positive for the marker; the threshold is derived on the events in `mask`."""
        x = C[det_of[name]].to_numpy()
        if name not in thr:
            thr[name] = fl_threshold(x[mask])
        return x > thr[name]

    def cut(name, values, mask, bw, sep):
        if name not in thr:
            t = valley(values[mask], bw, sep)
            thr[name] = np.inf if t is None else t
        return values > thr[name]

    n = len(df)
    lf = np.log10(df["FSC-A"].clip(lower=1).to_numpy())
    cells = cut("fsc", lf, np.ones(n, bool), 0.04, 0.3)
    if thr["fsc"] == np.inf:
        cells = np.ones(n, bool)
    if a.no_doublet_gate:
        singlets = cells.copy()
    else:
        r = df["FSC-A"].to_numpy() / df["FSC-H"].clip(lower=1).to_numpy()
        if "ratio" not in thr:
            med = np.median(r[cells])
            thr["ratio"] = med + 6 * 1.4826 * np.median(np.abs(r[cells] - med))
        singlets = cells & (r < thr["ratio"])
    if a.no_viability_gate:
        live = singlets.copy()
    else:
        live = singlets & ~marker("Viability", singlets)
    cd14 = marker("CD14", live)
    ls = np.log10(df["SSC-A"].clip(lower=1).to_numpy())
    high_ssc = cut("ssc", ls, live, 0.04, 0.25)
    mono = live & cd14
    lym = live & ~cd14 & ~high_ssc
    cd3 = marker("CD3", lym)
    t = lym & cd3
    cd4 = marker("CD4", t)
    cd8 = marker("CD8", t)
    non_t = lym & ~cd3
    cd19 = marker("CD19", non_t)
    b = non_t & cd19
    cd56 = marker("CD56", non_t & ~cd19)
    nk = non_t & ~cd19 & cd56
    return {"all_events": n, "cells": cells.sum(), "singlets": singlets.sum(), "live": live.sum(),
            "lymphocytes": lym.sum(), "monocytes": mono.sum(), "T_cells": t.sum(),
            "CD4_T": (t & cd4 & ~cd8).sum(), "CD8_T": (t & ~cd4 & cd8).sum(), "DP_T": (t & cd4 & cd8).sum(),
            "DN_T": (t & ~cd4 & ~cd8).sum(), "B_cells": b.sum(), "NK_cells": nk.sum()}


rows = []
template = {}
S_report = S_ctrl
for k, r in enumerate(samples.itertuples()):
    df, text = read_fcs(inp / "fcs" / r.fcs_file)
    if a.comp == "controls":
        S = S_ctrl
    elif a.comp == "embedded":
        S = embedded_matrix(text)
        S_report = S
    else:
        S = np.eye(len(dets))
    thr = template if a.global_thresholds else {}
    c = gate(df, S, thr)
    for pop, parent in POPS:
        rows.append((r.sample_id, pop, parent, int(c[pop]), 100.0 * c[pop] / max(int(c[parent]), 1)))

pd.DataFrame(S_report, index=pd.Index(fluors, name="fluorochrome"), columns=dets).to_csv(
    out / "spillover.csv", float_format="%.5f")
pd.DataFrame(rows, columns=["sample", "population", "parent", "count", "percent_of_parent"]).to_csv(
    out / "frequencies.tsv", sep="\t", index=False, float_format="%.3f")
print(f"wrote {len(rows)} rows for {len(samples)} samples; max off-diagonal spillover {np.max(S_report - np.eye(len(dets))):.3f}")
