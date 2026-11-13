#!/usr/bin/env python3
"""Independent second solution. Differs from the reference pipeline in every step:
  reader        fcsparser (not flowio)
  spillover     least-squares slope of every detector against the primary detector over all events
                of the single-stain tube (not medians of gated positive/negative populations)
  compensation  linear solve per sample
  gates         two-component Gaussian mixtures; an event is assigned by the fitted decision boundary
                (not by density valleys); lymphocytes/monocytes by a 2-D mixture on CD14 and SSC-A
Used to show that the verifier accepts a different but correct method.
"""
import argparse
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

panel = pd.read_csv(inp / "panel.tsv", sep="\t")
samples = pd.read_csv(inp / "samples.tsv", sep="\t")
controls = pd.read_csv(inp / "controls.tsv", sep="\t", keep_default_na=False)
dets = list(panel["detector"])
det_of = dict(zip(panel["marker"], panel["detector"]))


def load(name):
    _, df = fcsparser.parse(str(inp / "fcs" / name), channel_naming="$PnN", reformat_meta=False)
    return df.astype(float).reset_index(drop=True)


# ---- spillover by regression ---------------------------------------------------------------
S = np.eye(len(dets))
for r in controls[controls["control_type"] == "single_stain"].itertuples():
    X = load(r.fcs_file)[dets].to_numpy()
    X = X[(X < 250000).all(axis=1)]                    # no saturated events
    i = dets.index(r.detector)
    x = X[:, i] - X[:, i].mean()
    S[i] = (x[:, None] * (X - X.mean(axis=0))).sum(axis=0) / (x * x).sum()
    S[i, i] = 1.0
pd.DataFrame(S, index=pd.Index(list(panel["fluorochrome"]), name="fluorochrome"), columns=dets).to_csv(
    out / "spillover.csv", float_format="%.5f")


# ---- mixture-model gates --------------------------------------------------------------------
def boundary(y, lo_q=5.0, hi_q=99.9):
    """Fit a two-component 1-D Gaussian mixture (deterministic start: one component at a low and
    one at a high quantile) and return the point between the two means where the posterior
    switches. None if the two components are not separated (single population)."""
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
    """arcsinh with a data-driven cofactor: the linear region is as wide as the negative spread of
    the events under study (compensation spreads the unstained population around zero)."""
    cof = max(150.0, 2.0 * abs(np.percentile(x[mask], 1.0)))
    return np.arcsinh(x / cof)


def positive(x, mask, transform=True, **kw):
    y = scaled(x, mask) if transform else x
    t = boundary(y[mask], **kw)
    return np.zeros(len(y), bool) if t is None else y > t


rows = []
for r in samples.itertuples():
    df = load(r.fcs_file)
    raw = df[dets].to_numpy()
    comp = np.linalg.solve(S.T, raw.T).T
    M = {m: comp[:, dets.index(d)] for m, d in det_of.items()}
    n = len(df)
    fsc_a, fsc_h, ssc = (df[c].to_numpy() for c in ("FSC-A", "FSC-H", "SSC-A"))

    lf = np.log10(np.clip(fsc_a, 1, None))
    t = boundary(lf, lo_q=2.0, hi_q=60.0)
    cells = lf > t if t is not None else np.ones(n, bool)
    lr = np.log(np.clip(fsc_a, 1, None) / np.clip(fsc_h, 1, None))
    singlets = cells & ~positive(lr, cells, transform=False, lo_q=30.0, hi_q=99.5)
    live = singlets & ~positive(M["Viability"], singlets)

    # lymphocytes vs monocytes: 2-D mixture on (CD14, log SSC-A) of the live single cells
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
    c = {"all_events": n, "cells": cells.sum(), "singlets": singlets.sum(), "live": live.sum(),
         "lymphocytes": lym.sum(), "monocytes": mono.sum(), "T_cells": tc.sum(),
         "CD4_T": (tc & cd4 & ~cd8).sum(), "CD8_T": (tc & ~cd4 & cd8).sum(), "DP_T": (tc & cd4 & cd8).sum(),
         "DN_T": (tc & ~cd4 & ~cd8).sum(), "B_cells": (non_t & cd19).sum(), "NK_cells": (non_t & ~cd19 & cd56).sum()}
    parents = {"cells": "all_events", "singlets": "cells", "live": "singlets", "lymphocytes": "live",
               "monocytes": "live", "T_cells": "lymphocytes", "CD4_T": "T_cells", "CD8_T": "T_cells",
               "DP_T": "T_cells", "DN_T": "T_cells", "B_cells": "lymphocytes", "NK_cells": "lymphocytes"}
    for pop, parent in parents.items():
        rows.append((r.sample_id, pop, parent, int(c[pop]), round(100.0 * c[pop] / max(int(c[parent]), 1), 3)))

pd.DataFrame(rows, columns=["sample", "population", "parent", "count", "percent_of_parent"]).to_csv(
    out / "frequencies.tsv", sep="\t", index=False)
print(f"alt solution: {len(rows)} rows")
