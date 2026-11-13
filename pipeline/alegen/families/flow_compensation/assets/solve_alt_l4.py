#!/usr/bin/env python3
"""Independent second solution of the level-4 task. Differs from the reference pipeline in every step:
  reader            fcsparser (not flowio); detector voltages from its metadata dictionary
  control choice    by reagent lot of panel.tsv, then by detector settings
  spillover rows    least-squares slope of every detector against the primary detector over all
                    events of the single-stain tube (not medians of gated populations)
  gain factor       ratio of the positive-bead peaks of the two recordings, peaks located by a
                    two-component Gaussian mixture on the log signal (not background-subtracted medians)
  later matrix      column of the changed detector x gain factor; own row from the re-recorded tube
  gates             two-component Gaussian mixtures (not density valleys); lymphocytes/monocytes by a
                    2-D mixture on CD14 and SSC-A
  acceptance check  fine-tuning by matching population medians on the FMO tube, one coefficient per
                    population (not a regression over all lymphocytes)
  FMO percentile    nearest-rank definition (not numpy's interpolated percentile)
Used to show that the verifier accepts a different but correct method.
"""
import argparse
import math
import warnings
from pathlib import Path

import fcsparser
import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture

warnings.filterwarnings("ignore")
ap = argparse.ArgumentParser()
ap.add_argument("--input", default="/task/input")
ap.add_argument("--output", default="/task/output")
a = ap.parse_args()
inp, out = Path(a.input), Path(a.output)
out.mkdir(parents=True, exist_ok=True)

panel = pd.read_csv(inp / "panel.tsv", sep="\t", dtype=str)
samples = pd.read_csv(inp / "samples.tsv", sep="\t", dtype=str)
controls = pd.read_csv(inp / "controls.tsv", sep="\t", dtype=str, keep_default_na=False)
dets = list(panel["detector"])
fluors = list(panel["fluorochrome"])
det_of = dict(zip(panel["marker"], panel["detector"]))
lot_of = dict(zip(panel["fluorochrome"], panel["reagent_lot"]))
nd = len(dets)

_files = {}


def load(name):
    if name not in _files:
        meta, df = fcsparser.parse(str(inp / "fcs" / name), channel_naming="$PnN", reformat_meta=False)
        volt = {}
        for j in range(1, int(meta["$PAR"]) + 1):
            volt[str(meta[f"$P{j}N"])] = str(meta.get(f"$P{j}V", "")).strip()
        _files[name] = (df.astype(float).reset_index(drop=True), tuple(volt[d] for d in dets))
    return _files[name]


def regression_row(name, detector):
    X = load(name)[0][dets].to_numpy()
    X = X[(X < 250000).all(axis=1)]                    # no saturated events
    i = dets.index(detector)
    x = X[:, i] - X[:, i].mean()
    row = (x[:, None] * (X - X.mean(axis=0))).sum(axis=0) / (x * x).sum()
    row[i] = 1.0
    return row


def bead_peak(name, detector):
    """log10 position of the positive-bead peak in the primary detector (2-component mixture)."""
    y = np.log10(np.clip(load(name)[0][detector].to_numpy(), 10.0, None)).reshape(-1, 1)
    init = np.percentile(y, [20, 80]).reshape(2, 1)
    gm = GaussianMixture(2, means_init=init, random_state=0, max_iter=500).fit(y)
    return float(gm.means_.max())


# ---- controls: lot of the cocktail, grouped by detector settings ------------------------------
ss = controls[controls["control_type"] == "single_stain"]
by_fluor = {}
for r in ss.itertuples():
    if r.reagent_lot != lot_of[r.fluorochrome]:
        print(f"skipping {r.fcs_file}: lot {r.reagent_lot}, cocktail lot {lot_of[r.fluorochrome]}")
        continue
    by_fluor.setdefault(r.fluorochrome, {}).setdefault(load(r.fcs_file)[1], []).append(r)
sample_set = {r.sample_id: load(r.fcs_file)[1] for r in samples.itertuples()}
# session-start settings = the settings under which every fluorochrome has a control
all_sets = [set(v) for v in by_fluor.values()]
common = set.intersection(*all_sets)
if len(common) != 1 or len(by_fluor) != nd:
    raise SystemExit("cannot identify one complete control set")
base = next(iter(common))

S1 = np.eye(nd)
for i, f in enumerate(fluors):
    S1[i] = regression_row(by_fluor[f][base][0].fcs_file, dets[i])
mats = {base: S1}
for sett in set(sample_set.values()) - {base}:
    M = S1.copy()
    for i in range(nd):
        if sett[i] == base[i]:
            continue
        f = fluors[i]
        if sett not in by_fluor[f]:
            raise SystemExit(f"no {f} control for the changed setting of {dets[i]}")
        new, old = by_fluor[f][sett][0].fcs_file, by_fluor[f][base][0].fcs_file
        g = 10 ** (bead_peak(new, dets[i]) - bead_peak(old, dets[i]))
        print(f"{dets[i]}: {base[i]} V -> {sett[i]} V, gain factor {g:.4f}")
        M[:, i] *= g
        M[i] = regression_row(new, dets[i])
    mats[sett] = M


# ---- mixture-model gates ---------------------------------------------------------------------
def boundary(y, lo_q=5.0, hi_q=99.9):
    y = np.asarray(y, dtype=float).reshape(-1, 1)
    if len(y) < 100:
        return None
    m0 = np.percentile(y, [lo_q, hi_q])
    if m0[1] - m0[0] < 1e-9:
        return None
    sd0 = (m0[1] - m0[0]) / 8.0
    gm = GaussianMixture(2, means_init=m0.reshape(2, 1), weights_init=[0.8, 0.2],
                         precisions_init=np.full((2, 1, 1), 1.0 / sd0 ** 2),
                         random_state=0, max_iter=1000, tol=1e-6).fit(y)
    mu = gm.means_.ravel()
    sd = np.sqrt(gm.covariances_.ravel())
    lo, hi = int(np.argmin(mu)), int(np.argmax(mu))
    if mu[hi] - mu[lo] < 2.0 * max(sd[lo], sd[hi]) or gm.weights_.min() < 2e-4:
        return None
    grid = np.linspace(mu[lo], mu[hi], 2001).reshape(-1, 1)
    post = gm.predict_proba(grid)[:, hi]
    return float(grid[int(np.argmax(post >= 0.5)), 0])


def scaled(x, mask):
    cof = max(150.0, 2.0 * abs(np.percentile(x[mask], 1.0)))
    return np.arcsinh(x / cof)


def positive(x, mask, transform=True, **kw):
    y = scaled(x, mask) if transform else x
    t = boundary(y[mask], **kw)
    return np.zeros(len(y), bool) if t is None else y > t


def gates(name):
    df, sett = load(name)
    comp = np.linalg.solve(mats[sett].T, df[dets].to_numpy().T).T
    M = {m: comp[:, dets.index(d)] for m, d in det_of.items()}
    n = len(df)
    fsc_a, fsc_h, ssc = (df[c].to_numpy() for c in ("FSC-A", "FSC-H", "SSC-A"))
    lf = np.log10(np.clip(fsc_a, 1, None))
    t = boundary(lf, lo_q=2.0, hi_q=60.0)
    cells = lf > t if t is not None else np.ones(n, bool)
    lr = np.log(np.clip(fsc_a, 1, None) / np.clip(fsc_h, 1, None))
    singlets = cells & ~positive(lr, cells, transform=False, lo_q=30.0, hi_q=99.5)
    live = singlets & ~positive(M["Viability"], singlets)
    Z = np.column_stack([scaled(M["CD14"], live), 4.0 * np.log10(np.clip(ssc, 1, None))])[live]
    init = np.vstack([np.percentile(Z, 20, axis=0), np.percentile(Z, 99, axis=0)])
    gm = GaussianMixture(2, means_init=init, random_state=0, max_iter=1000).fit(Z)
    mono_comp = int(np.argmax(gm.means_[:, 0]))
    is_mono = np.zeros(n, bool)
    if abs(gm.means_[0, 0] - gm.means_[1, 0]) > 2.0:
        is_mono[np.flatnonzero(live)] = gm.predict(Z) == mono_comp
    mono = live & is_mono
    lym = live & ~is_mono
    cd3 = positive(M["CD3"], lym)
    tc = lym & cd3
    cd4 = positive(M["CD4"], tc)
    cd8 = positive(M["CD8"], tc)
    non_t = lym & ~cd3
    cd19 = positive(M["CD19"], non_t)
    cd56 = positive(M["CD56"], non_t & ~cd19)
    masks = {"cells": cells, "singlets": singlets, "live": live, "lymphocytes": lym, "monocytes": mono,
             "T_cells": tc, "CD4_T": tc & cd4 & ~cd8, "CD8_T": tc & ~cd4 & cd8, "DP_T": tc & cd4 & cd8,
             "DN_T": tc & ~cd4 & ~cd8, "B_cells": non_t & cd19, "NK_cells": non_t & ~cd19 & cd56}
    return masks, comp, sett, n


def nearest_rank(v, pct):
    v = np.sort(np.asarray(v, dtype=float))
    return float(v[max(0, math.ceil(pct / 100.0 * len(v)) - 1)])


# ---- FMO tubes: acceptance check, fine-tuning, thresholds per detector setting ------------------
ACT = {"activated_CD4_T": "CD4_T", "activated_CD8_T": "CD8_T"}
QC = ("CD4_T", "CD8_T", "DN_T", "B_cells", "NK_cells")
# population -> marker whose fluorochrome it carries and DN_T does not
CARRIES = {"CD4_T": "CD4", "CD8_T": "CD8", "B_cells": "CD19", "NK_cells": "CD56"}
thr = {}
for r in controls[controls["control_type"] == "fmo"].itertuples():
    j = dets.index(r.detector)
    for it in range(4):
        masks, comp, sett, _ = gates(r.fcs_file)
        x = comp[:, j]
        med = {p: float(np.median(x[masks[p]])) for p in QC}
        spread = max(med.values()) - min(med.values())
        print(f"{r.fcs_file} pass {it}: FMO medians " + "/".join(f"{med[p]:.1f}" for p in QC) + f", spread {spread:.1f}")
        if it == 3:
            break
        # a population that sits above DN_T is under-compensated for the fluorochrome it carries
        for pop, marker in CARRIES.items():
            i = dets.index(det_of[marker])
            mats[sett][i, j] += (med[pop] - med["DN_T"]) / float(np.median(comp[masks[pop], i]))
    thr[sett] = (j, {pop: nearest_rank(x[masks[par]], 99.0) for pop, par in ACT.items()})

parents = {"cells": "all_events", "singlets": "cells", "live": "singlets", "lymphocytes": "live",
           "monocytes": "live", "T_cells": "lymphocytes", "CD4_T": "T_cells", "CD8_T": "T_cells",
           "DP_T": "T_cells", "DN_T": "T_cells", "B_cells": "lymphocytes", "NK_cells": "lymphocytes",
           "activated_CD4_T": "CD4_T", "activated_CD8_T": "CD8_T"}
rows, spill = [], []
for r in samples.itertuples():
    masks, comp, sett, n = gates(r.fcs_file)
    if sett not in thr:
        raise SystemExit(f"no FMO tube for the detector settings of {r.sample_id}")
    j, t = thr[sett]
    c = {k: int(v.sum()) for k, v in masks.items()}
    c["all_events"] = n
    for pop, par in ACT.items():
        c[pop] = int((masks[par] & (comp[:, j] > t[pop])).sum())
    for pop, parent in parents.items():
        rows.append((r.sample_id, pop, parent, c[pop], round(100.0 * c[pop] / max(c[parent], 1), 3)))
    for i, f in enumerate(fluors):
        spill.append([r.sample_id, f] + [round(float(v), 6) for v in mats[sett][i]])

pd.DataFrame(spill, columns=["sample", "fluorochrome"] + dets).to_csv(out / "spillover.csv", index=False)
pd.DataFrame(rows, columns=["sample", "population", "parent", "count", "percent_of_parent"]).to_csv(
    out / "frequencies.tsv", sep="\t", index=False)
print(f"alt solution: {len(rows)} rows, {len(mats)} detector setting(s)")
