#!/usr/bin/env python3
"""Reference analysis of the level-7 session: the level-4 analysis (8-colour PBMC panel, detector
setting changed during the run, one bead control from another reagent lot, FMO-gated HLA-DR
read-out) preceded by a signal-versus-time check of every tube.

Steps
  1  read the keywords of every file; group the files by the voltages of the fluorescence detectors
  2  signal-versus-time check of every tube with enough events: median and 90th percentile of every
     fluorescence detector and of FSC-A / SSC-A in consecutive bins of 500 events ordered by Time; a
     run of bins that leaves the level of the rest of the file marks a disturbed stretch
  3  per detector, the factor between the signal inside and outside the stretch (median over the
     50th-95th percentiles); the detectors whose factor differs from 1 are the affected ones
  4  first and last event of the stretch (change points of an event-wise score: log ratio of the
     histograms of the affected detectors outside and inside the stretch); the factor of the
     affected detectors is refined (shift between the peaks of the log signal inside and outside);
     the events of the stretch are set aside (default), with a guard of half a bin on both sides
  5  choose, per fluorochrome, the single-stain control that was stained with the reagent lot of the
     sample cocktail and recorded with the session-start settings
  6  spillover rows from those controls (median of the positive minus median of the negative
     population) -> matrix of the session-start settings
  7  for a later detector setting: gain factor of the changed detector from the two recordings of
     the same bead tube (ratio of the background-subtracted bead signal), cross-checked with the
     ratio of the two spillover rows
  8  matrix of the later setting: the column of the changed detector multiplied by the gain factor
     (every fluorochrome's signal in that detector scales with it), the row of its own
     fluorochrome taken from the re-recorded control
  9-15 per tube: compensation, debris, singlets, viability, lymphocytes/monocytes, T/B/NK,
     CD4/CD8 gates (thresholds at density valleys of the undisturbed events of that tube)
  16 FMO tube of each setting gated like a sample; acceptance check of the SOP (medians of the
     compensated HLA-DR detector signal of CD4_T, CD8_T, DN_T, B_cells, NK_cells)
  17 fine-tuning on the FMO tube: the tube holds no reagent for the HLA-DR detector, so the median
     compensated signal of that detector must not depend on which other fluorochromes a lymphocyte
     population carries; a weighted least-squares fit across seven populations (median of the
     HLA-DR detector against the medians of the other fluorochromes) gives the errors of the spill
     coefficients into that detector, which are removed (3 passes); acceptance check repeated
  18 99th percentile of the compensated HLA-DR detector signal of CD4_T and CD8_T of the FMO tube
  19 HLA-DR + events of CD4_T / CD8_T of every sample with the FMO threshold of its setting
  20 population sizes. A population whose gates use neither an affected detector nor a detector
     that is compensated for one by 2% or more is counted on all events of the file. Every other
     population is estimated from the undisturbed events: its share of the parent population
     there, times the size of the parent in the file
  21 output tables

The non-default flags replay careless analyses and alternative treatments of the disturbed stretch
for the validator's battery; the oracle (solve.sh) uses the defaults.
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
ap.add_argument("--blocks", choices=["settings", "single", "row_only"], default="settings")
ap.add_argument("--tandem", choices=["lot", "first", "last"], default="lot")
ap.add_argument("--dr-gate", choices=["fmo", "valley", "unstained", "fmo_raw"], default="fmo")
ap.add_argument("--no-doublet-gate", action="store_true")
ap.add_argument("--no-viability-gate", action="store_true")
ap.add_argument("--global-thresholds", action="store_true")
ap.add_argument("--refine", choices=["auto", "on", "off"], default="auto",
                help="FMO-based acceptance check with fine-tuning; auto = on unless a careless flag is set")
ap.add_argument("--l3", action="store_true", help="level-3 reference method: one matrix from the first-listed "
                "controls, valley rule for every marker")
ap.add_argument("--time-qc", choices=["exclude", "after", "repair", "drop", "wrong", "off"], default="exclude",
                help="what is done with a disturbed stretch found by the signal-versus-time check. exclude: set "
                     "aside, populations that need the affected detectors are estimated from the other events; "
                     "after: as exclude, but everything from the start of the stretch to the end of the file is "
                     "set aside; repair: the affected detectors are divided by their factor inside the stretch, "
                     "all events are used; drop: the events are removed and the sizes of what is left are "
                     "reported; wrong: a stretch of the same length is set aside in another part of the file; "
                     "off: no check")
a = ap.parse_args()
if a.l3:
    a.blocks, a.dr_gate = "single", "valley"
    a.tandem = "first" if a.tandem == "lot" else a.tandem
careless = a.comp != "controls" or a.blocks != "settings" or a.tandem != "lot"
do_refine = a.comp == "controls" and (a.refine == "on" or (a.refine == "auto" and not careless))
inp, out = Path(a.input), Path(a.output)
out.mkdir(parents=True, exist_ok=True)

panel = pd.read_csv(inp / "panel.tsv", sep="\t", dtype=str)
samples = pd.read_csv(inp / "samples.tsv", sep="\t", dtype=str)
controls = pd.read_csv(inp / "controls.tsv", sep="\t", dtype=str, keep_default_na=False)
dets = list(panel["detector"])
fluors = list(panel["fluorochrome"])
det_of = dict(zip(panel["marker"], panel["detector"]))
lot_of = dict(zip(panel["fluorochrome"], panel["reagent_lot"])) if "reagent_lot" in panel else {}
COF = 150.0
FMO_PCT = 99.0
BIN, MIN_RUN = 500, 5            # signal-versus-time check: events per bin, shortest run of deviating bins
GUARD = BIN // 2                 # events set aside on both sides of a disturbed stretch

POPS = [("cells", "all_events"), ("singlets", "cells"), ("live", "singlets"), ("lymphocytes", "live"),
        ("monocytes", "live"), ("T_cells", "lymphocytes"), ("CD4_T", "T_cells"), ("CD8_T", "T_cells"),
        ("DP_T", "T_cells"), ("DN_T", "T_cells"), ("B_cells", "lymphocytes"), ("NK_cells", "lymphocytes"),
        ("activated_CD4_T", "CD4_T"), ("activated_CD8_T", "CD8_T")]
ACT = {"activated_CD4_T": "CD4_T", "activated_CD8_T": "CD8_T"}
QC_POPS = ["CD4_T", "CD8_T", "DN_T", "B_cells", "NK_cells"]     # populations of the SOP's acceptance check
# markers whose gates decide the membership in a population (beyond those of its parent)
GATE_MARKERS = {"cells": [], "singlets": [], "live": ["Viability"], "lymphocytes": ["CD14"], "monocytes": ["CD14"],
                "T_cells": ["CD3"], "CD4_T": ["CD4", "CD8"], "CD8_T": ["CD4", "CD8"], "DP_T": ["CD4", "CD8"],
                "DN_T": ["CD4", "CD8"], "B_cells": ["CD3", "CD19"], "NK_cells": ["CD3", "CD19", "CD56"],
                "activated_CD4_T": ["HLA-DR"], "activated_CD8_T": ["HLA-DR"]}

_cache = {}


def read_fcs(name):
    """-> (events as DataFrame, keywords with lower-case names without '$', voltage per detector)"""
    if name not in _cache:
        fd = flowio.FlowData(str(inp / "fcs" / name))
        arr = np.asarray(fd.as_array(preprocess=False), dtype=float)
        names = list(fd.pnn_labels)
        text = {str(k).lower().lstrip("$"): v for k, v in fd.text.items()}
        volts = {names[j]: str(text.get(f"p{j + 1}v", "")).strip() for j in range(len(names))}
        _cache[name] = (pd.DataFrame(arr, columns=names), text, volts)
    return _cache[name]


def settings(name):
    """Voltages of the fluorescence detectors: files with the same tuple share one spillover matrix."""
    v = read_fcs(name)[2]
    return tuple(v[d] for d in dets)


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


# ----------------------------------------------------------------------------- signal versus time
def longest_run(flag):
    """(first, last) index of the longest run of True values, or None."""
    best, start = None, None
    for i, f in enumerate(list(flag) + [False]):
        if f and start is None:
            start = i
        elif not f and start is not None:
            if best is None or i - start > best[1] - best[0] + 1:
                best = (start, i - 1)
            start = None
    return best


def scale_factor(x_in, x_out, r0):
    """Factor between two recordings of the same signal, one of them scaled by about r0: the shift
    that brings the densities of the log signal (above the noise floor) into register. A ratio of
    percentiles depends on how many cells of each kind the two sets happen to hold; the position of
    the peaks does not."""
    step, floor = 0.01, 300.0
    grid = np.arange(np.log(floor), np.log(262144.0) + step, step)

    def density(x):
        h = np.histogram(np.log(x[x > floor]), bins=grid)[0] / len(x)
        return gaussian_filter1d(h.astype(float), 5.0, mode="constant")

    h_in, h_out = density(x_in / r0), density(x_out)
    lags = np.arange(-30, 31)
    k = len(h_out)
    cc = np.array([np.dot(h_in[max(l, 0):k + min(l, 0)], h_out[max(-l, 0):k - max(l, 0)]) for l in lags])
    i = int(np.argmax(cc))
    if cc[i] <= 0 or i in (0, len(lags) - 1):
        return r0
    y0, y1, y2 = cc[i - 1], cc[i], cc[i + 1]
    peak = lags[i] + 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2)            # vertex of the parabola through three lags
    return float(r0 * np.exp(step * peak))


def find_stretch(df):
    """Signal-versus-time check of one tube. None if the signals are stable over the run; otherwise
    the disturbed stretch: `order` (rows in time order), `first` / `end` (positions in that order of
    its first event and of the first event after it), `factor` (signal inside / outside per
    parameter) and `affected` (parameters whose factor differs from 1 by more than 15%)."""
    params = dets + ["FSC-A", "SSC-A"]
    n = len(df)
    nb = n // BIN
    if nb < 4 * MIN_RUN:
        return None                                  # too short for a binned check (bead tubes)
    order = np.argsort(df["Time"].to_numpy(), kind="stable")
    X = df[params].to_numpy()[order]
    edges = np.linspace(0, n, nb + 1).astype(int)
    best = None
    for q in (50, 90):
        # level of every bin on a scale that is linear near zero and logarithmic above the cofactor
        v = np.arcsinh(np.array([np.percentile(X[edges[b]:edges[b + 1]], q, axis=0) for b in range(nb)]) / COF)
        base = np.median(v, axis=0)
        sd = 1.4826 * np.median(np.abs(v - base), axis=0)
        for j in range(len(params)):
            for sign in (-1.0, 1.0):
                flag = sign * (v[:, j] - base[j]) > max(4.0 * sd[j], 0.1)
                run = longest_run(flag)
                if run is None or run[1] - run[0] + 1 < MIN_RUN:
                    continue
                elsewhere = flag.copy()
                elsewhere[max(run[0] - 1, 0):run[1] + 2] = False
                if elsewhere.sum() > 0.1 * (nb - (run[1] - run[0] + 1)):
                    continue                         # a level that jumps all over the run is not a stretch
                size = (run[1] - run[0] + 1, float(np.mean(np.abs(v[run[0]:run[1] + 1, j] - base[j]))))
                if best is None or size > best[0]:
                    best = (size, run)
    if best is None:
        return None
    b0, b1 = best[1]
    inside = np.arange(edges[b0 + 1], edges[b1])                      # bins that lie wholly in the stretch
    outside = np.r_[0:edges[max(b0 - 1, 0)], edges[min(b1 + 2, nb)]:n]
    levels = np.arange(50, 96, 5)
    factor = {}
    for j, prm in enumerate(params):
        qi, qo = np.percentile(X[inside, j], levels), np.percentile(X[outside, j], levels)
        ok = (qo > 100.0) & (qi > 0.0)                                # percentiles above the noise floor
        factor[prm] = float(np.exp(np.median(np.log(qi[ok] / qo[ok])))) if ok.sum() >= 3 else 1.0
    affected = [prm for prm in params if abs(np.log(factor[prm])) > np.log(1.15)]
    if not affected:
        return None
    # change points: every event is scored by how much better its affected signals fit the rest of
    # the tube than the inside of the stretch (log ratio of the two histograms, summed over the
    # affected parameters); the running sum of the scores rises up to the first disturbed event and
    # falls up to the last one
    c = np.zeros(n)
    for prm in affected:
        u = np.arcsinh(X[:, params.index(prm)] / COF)
        k = np.clip(((u - u.min()) / (u.max() - u.min() + 1e-9) * 60).astype(int), 0, 59)
        p_in = (np.bincount(k[inside], minlength=60) + 1.0) / (len(inside) + 60.0)
        p_out = (np.bincount(k[outside], minlength=60) + 1.0) / (len(outside) + 60.0)
        c += np.log(p_out[k] / p_in[k])
    half = edges[(b0 + b1 + 1) // 2]                                  # an event in the middle of the run
    lo = edges[max(b0 - 3, 0)]
    first = int(lo + np.argmax(np.concatenate([[0.0], np.cumsum(c[lo:half])])))
    hi = edges[min(b1 + 4, nb)]
    end = int(half + np.argmin(np.concatenate([[0.0], np.cumsum(c[half:hi])])))
    for prm in affected:                                              # factor of the affected parameters, refined
        j = params.index(prm)
        factor[prm] = scale_factor(X[first:end, j], np.concatenate([X[:max(first - GUARD, 0), j], X[end + GUARD:, j]]),
                                   factor[prm])
    return {"order": order, "first": first, "end": end, "factor": factor, "affected": affected,
            "bins": (int(b0), int(b1), int(nb))}


_tubes, qc_time = {}, []


def tube(name):
    """Events of a tube after the signal-versus-time check: (events, mask of the events that are
    used for gates and estimates, the disturbed stretch or None)."""
    if name in _tubes:
        return _tubes[name]
    df = read_fcs(name)[0]
    keep = np.ones(len(df), bool)
    st = None if a.time_qc == "off" else find_stretch(df)
    if st is not None:
        n, order, first, end = len(df), st["order"], st["first"], st["end"]
        tm = df["Time"].to_numpy()[order]
        qc_time.append([name, n, first, end - 1, tm[first], tm[end - 1], round((end - first) / n, 4),
                        st["bins"][1] - st["bins"][0] + 1, st["bins"][2],
                        ";".join(f"{prm}:{st['factor'][prm]:.4f}" for prm in st["affected"]),
                        ";".join(f"{prm}:{f:.4f}" for prm, f in st["factor"].items() if prm not in st["affected"])])
        print(f"[time] {name}: events {first}-{end - 1} of {n} (Time {tm[first]:.0f}-{tm[end - 1]:.0f}) "
              + ", ".join(f"{prm} x{st['factor'][prm]:.3f}" for prm in st["affected"]) + f" -> {a.time_qc}")
        if a.time_qc == "repair":
            df = df.copy()
            rows = order[first:end]
            for prm in st["affected"]:
                df.loc[rows, prm] = df.loc[rows, prm] / st["factor"][prm]
        elif a.time_qc == "after":
            keep[order[max(first - GUARD, 0):]] = False
        elif a.time_qc == "wrong":
            length = end - first + 2 * GUARD
            if first - GUARD >= n - (end + GUARD):
                keep[order[:length]] = False          # more room before the stretch than after it
            else:
                keep[order[n - length:]] = False
        else:
            keep[order[max(first - GUARD, 0):end + GUARD]] = False
        if a.time_qc == "drop":
            df, keep, st = df[keep].reset_index(drop=True), np.ones(int(keep.sum()), bool), None
    _tubes[name] = (df, keep, st)
    return _tubes[name]


def kept(name):
    """Events of a control tube that pass the signal-versus-time check."""
    df, keep, _ = tube(name)
    return df if keep.all() else df[keep].reset_index(drop=True)


# ----------------------------------------------------------------------------- spillover matrices
def control_row(name, detector):
    """Spillover row of a single-stain tube and the background-subtracted signal of its positive
    population in the primary detector."""
    df = kept(name)
    lf = np.log10(df["FSC-A"].clip(lower=1))
    t = valley(lf, 0.04, 0.3)                      # drop debris of cell-based controls
    if t is not None:
        df = df[lf > t]
    i = dets.index(detector)
    y = df[detector].to_numpy()
    t = fl_threshold(y)
    if not np.isfinite(t):
        raise SystemExit(f"control {name}: no positive population found")
    D = df[dets].to_numpy()
    mp, mn = np.median(D[y > t], axis=0), np.median(D[y <= t], axis=0)
    row = (mp - mn) / (mp[i] - mn[i])
    row[i] = 1.0
    return row, float(mp[i] - mn[i])


ss = controls[controls["control_type"] == "single_stain"].reset_index(drop=True)
base = settings(ss.iloc[0]["fcs_file"])            # settings of the session-start control set
ss_set = [settings(f) for f in ss["fcs_file"]]


def pick(fluor, sett, like=None):
    """Single-stain control of `fluor` recorded with `sett`. Default: the one stained with the lot of
    the sample cocktail; --tandem first: the first one listed. `like`: same carrier and lot as that row."""
    if a.tandem == "last" and like is None:         # last-listed control, whatever its settings and lot
        cand = ss[ss["fluorochrome"] == fluor]
        return cand.iloc[-1] if len(cand) else None
    cand = ss[(ss["fluorochrome"] == fluor).to_numpy() & np.array([s == sett for s in ss_set])]
    has_lot = "reagent_lot" in ss.columns
    if like is not None:
        same = cand[cand["carrier"] == like["carrier"]]
        if has_lot:
            same = same[same["reagent_lot"] == like["reagent_lot"]]
        cand = same if len(same) else cand
    elif a.tandem == "lot" and has_lot and fluor in lot_of:
        match = cand[cand["reagent_lot"] == lot_of[fluor]]
        cand = match if len(match) else cand
    return cand.iloc[0] if len(cand) else None


S_base = np.eye(len(dets))
chosen = {}
for i, f in enumerate(fluors):
    r = pick(f, base)
    if r is None:
        raise SystemExit(f"no session-start single-stain control for {f}")
    chosen[f] = r
    S_base[i] = control_row(r["fcs_file"], r["detector"])[0]
    print(f"[matrix] {f:<12} from {r['fcs_file']} ({r['carrier']}, lot {r.get('reagent_lot', 'NA')})")

sample_set = {r.sample_id: settings(r.fcs_file) for r in samples.itertuples()}
mats = {base: S_base}
for sett in dict.fromkeys(sample_set.values()):
    if sett in mats:
        continue
    if a.blocks == "single":
        mats[sett] = S_base
        continue
    M = S_base.copy()
    for i, d in enumerate(dets):
        if sett[i] == base[i]:
            continue
        f = fluors[i]
        rer = pick(f, sett, like=chosen[f])
        if rer is None:
            raise SystemExit(f"voltage of {d} changed ({base[i]} -> {sett[i]}) and no control of {f} was "
                             "recorded with the new setting")
        row_new, amp_new = control_row(rer["fcs_file"], d)
        row_old, amp_old = control_row(chosen[f]["fcs_file"], d)
        g = amp_new / amp_old                       # same bead tube before / after: detector gain factor
        o = np.arange(len(dets)) != i
        g_rows = float((row_old[o] @ row_new[o]) / (row_new[o] @ row_new[o]))   # row of f scales with 1/g
        print(f"[gain] {d}: {base[i]} V -> {sett[i]} V, gain x{g:.4f} from the bead signal "
              f"(x{g_rows:.4f} from the ratio of the two {f} rows)")
        if a.blocks == "settings":
            M[:, i] *= g                            # spill of every fluorochrome into d scales with its gain
        M[i] = row_new
        M[i, i] = 1.0
    mats[sett] = M


def embedded_matrix(text):
    key = next((k for k in text if k in ("spillover", "spill")), None)
    if key is None:
        return np.eye(len(dets))
    parts = text[key].split(",")
    n = int(parts[0])
    names = parts[1:1 + n]
    M = np.array(parts[1 + n:1 + n + n * n], dtype=float).reshape(n, n)
    idx = [names.index(d) for d in dets]
    return M[np.ix_(idx, idx)]


def matrix_for(name):
    """Matrix used to compensate the tube `name` (and reported for it)."""
    if a.comp == "embedded":
        return embedded_matrix(read_fcs(name)[1])
    sett = settings(name)
    return mats.get(sett, S_base)


# ----------------------------------------------------------------------------- gating
def gate(df, S, thr, keep=None):
    """Hierarchical gate. `thr` is a dict of thresholds; missing entries are derived from the events
    `keep` of this tube and stored (so that a filled dict can be re-used as a fixed template).
    Returns masks over all events and the compensated signals."""
    C = pd.DataFrame(df[dets].to_numpy() @ np.linalg.inv(S), columns=dets)
    n = len(df)
    keep = np.ones(n, bool) if keep is None else keep

    def marker(name, mask):
        x = C[det_of[name]].to_numpy()
        if name not in thr:
            thr[name] = fl_threshold(x[mask & keep])
        return x > thr[name]

    def cut(name, values, mask, bw, sep):
        if name not in thr:
            t = valley(values[mask & keep], bw, sep)
            thr[name] = np.inf if t is None else t
        return values > thr[name]

    lf = np.log10(df["FSC-A"].clip(lower=1).to_numpy())
    cells = cut("fsc", lf, np.ones(n, bool), 0.04, 0.3)
    if thr["fsc"] == np.inf:
        cells = np.ones(n, bool)
    if a.no_doublet_gate:
        singlets = cells.copy()
    else:
        r = df["FSC-A"].to_numpy() / df["FSC-H"].clip(lower=1).to_numpy()
        if "ratio" not in thr:
            med = np.median(r[cells & keep])
            thr["ratio"] = med + 6 * 1.4826 * np.median(np.abs(r[cells & keep] - med))
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
    masks = {"cells": cells, "singlets": singlets, "live": live, "lymphocytes": lym, "monocytes": mono,
             "T_cells": t, "CD4_T": t & cd4 & ~cd8, "CD8_T": t & ~cd4 & cd8, "DP_T": t & cd4 & cd8,
             "DN_T": t & ~cd4 & ~cd8, "B_cells": b, "NK_cells": nk}
    return masks, C


S_gate = (lambda name: np.eye(len(dets))) if a.comp == "none" else matrix_for


def acceptance(C, masks, detector):
    """Medians of the SOP's acceptance check and their spread (nan = not evaluable because a
    population is missing, which is itself a failed check)."""
    cd = C[detector].to_numpy()
    med = [float(np.median(cd[masks[p]])) if masks[p].sum() >= 20 else float("nan") for p in QC_POPS]
    return med, (float("nan") if any(np.isnan(med)) else max(med) - min(med))


def fine_tune(df, S, masks, detector, passes=3):
    """Spill of the other fluorochromes into `detector`, corrected on a tube that holds no reagent
    for it. Fit across lymphocyte populations, not across single events: population medians are
    free of the measurement noise of the single event, which biases an event-level regression of
    one compensated signal on another (the noise of the regressor is part of both)."""
    groups = {k: masks[k] for k in ("CD4_T", "CD8_T", "DP_T", "DN_T", "B_cells", "NK_cells")}
    groups["lineage_negative"] = masks["lymphocytes"] & ~masks["T_cells"] & ~masks["B_cells"] & ~masks["NK_cells"]
    groups = [m for m in groups.values() if m.sum() >= 100]
    X = df[dets].to_numpy()
    j = dets.index(detector)
    w = np.sqrt([m.sum() for m in groups])
    S = S.copy()
    for _ in range(passes):
        C = X @ np.linalg.inv(S)
        med = np.array([np.median(C[m], axis=0) for m in groups])           # populations x detectors
        use = [i for i in range(len(dets)) if i != j and np.ptp(med[:, i]) > 1000]
        if len(use) + 1 > len(groups):
            break
        A = np.column_stack([med[:, use], np.ones(len(groups))]) * w[:, None]
        S[use, j] += np.linalg.lstsq(A, med[:, j] * w, rcond=None)[0][:-1]
    return S


# ---- FMO tubes: acceptance check, fine-tuning, thresholds of the HLA-DR gate --------------------
fmo = controls[controls["control_type"] == "fmo"].reset_index(drop=True)
fmo_thr, qc_rows = {}, []
for r in fmo.itertuples():
    df = kept(r.fcs_file)
    sett = settings(r.fcs_file)
    masks, C = gate(df, S_gate(r.fcs_file), {})
    if a.blocks == "single" and fmo_thr:
        continue                                    # one FMO gate for the whole session
    med0, spread0 = acceptance(C, masks, r.detector)
    if do_refine and masks["lymphocytes"].sum() >= 1000:
        mats[sett] = fine_tune(df, mats.get(sett, S_base), masks, r.detector)
        masks, C = gate(df, S_gate(r.fcs_file), {})
    med, spread = acceptance(C, masks, r.detector)
    x = (df if a.dr_gate == "fmo_raw" else C)[r.detector].to_numpy()
    fmo_thr[sett] = {pop: float(np.percentile(x[masks[par]], FMO_PCT)) if masks[par].sum() >= 50 else np.inf
                     for pop, par in ACT.items()}
    qc_rows.append([r.fcs_file, r.detector, sett[dets.index(r.detector)], spread0] + med +
                   [spread, fmo_thr[sett]["activated_CD4_T"], fmo_thr[sett]["activated_CD8_T"]])
    print(f"[fmo] {r.fcs_file}: acceptance check before fine-tuning: spread {spread0:.1f}; final medians of "
          f"compensated {r.detector} " + ", ".join(f"{p} {m:.1f}" for p, m in zip(QC_POPS, med))
          + f" -> spread {spread:.1f}; thresholds "
          f"{fmo_thr[sett]['activated_CD4_T']:.1f} / {fmo_thr[sett]['activated_CD8_T']:.1f}")
dr_det = fmo.iloc[0]["detector"] if len(fmo) else det_of.get("HLA-DR")
first_fmo = next(iter(fmo_thr.values())) if fmo_thr else None

unstained_thr = None
if a.dr_gate == "unstained":
    u = controls[controls["control_type"] == "unstained"].iloc[0]
    df = kept(u["fcs_file"])
    lf = np.log10(df["FSC-A"].clip(lower=1).to_numpy())
    t = valley(lf, 0.04, 0.3)
    keep = lf > t if t is not None else np.ones(len(df), bool)
    x = (df[dets].to_numpy() @ np.linalg.inv(S_gate(u["fcs_file"])))[:, dets.index(dr_det)]
    unstained_thr = float(np.percentile(x[keep], FMO_PCT))
    print(f"[unstained] HLA-DR threshold {unstained_thr:.1f}")

def needs_affected(st, M):
    """Populations whose own gates use a parameter that the disturbed stretch touches: an affected
    detector, or a detector whose compensated signal takes 2% or more of an affected one."""
    if any(prm not in dets for prm in st["affected"]):
        return {pop for pop, _ in POPS}              # scatter itself is off: nothing can be counted there
    aff = [dets.index(prm) for prm in st["affected"]]
    Minv = np.abs(np.linalg.inv(M))
    hit = {d for j, d in enumerate(dets) if j in aff or Minv[aff, j].max() >= 0.02}
    return {pop for pop, _ in POPS if any(det_of[m] in hit for m in GATE_MARKERS[pop])}


# ---- samples ------------------------------------------------------------------------------------
rows, spill_rows = [], []
template = {}
for r in samples.itertuples():
    df, keep, st = tube(r.fcs_file)
    thr = template if a.global_thresholds else {}
    masks, C = gate(df, S_gate(r.fcs_file), thr, keep)
    sett = sample_set[r.sample_id]
    x = (df if a.dr_gate == "fmo_raw" else C)[dr_det].to_numpy()
    for pop, par in ACT.items():
        if a.dr_gate in ("fmo", "fmo_raw"):
            t = (fmo_thr.get(sett) or first_fmo)[pop] if a.blocks != "single" else first_fmo[pop]
        elif a.dr_gate == "unstained":
            t = unstained_thr
        else:
            t = fl_threshold(x[masks[par] & keep])
        masks[pop] = masks[par] & (x > t)
    masks["all_events"] = np.ones(len(df), bool)
    c = {"all_events": len(df)}
    indirect = set() if keep.all() else needs_affected(st, S_gate(r.fcs_file))
    estimated = set()
    for pop, parent in POPS:
        if pop in indirect or parent in estimated:
            # share of the parent among the undisturbed events, times the size of the parent in the file
            share = (masks[pop] & keep).sum() / max((masks[parent] & keep).sum(), 1)
            c[pop] = int(round(share * c[parent]))
            estimated.add(pop)
        else:
            c[pop] = int(masks[pop].sum())
    if estimated:
        print(f"[time] {r.sample_id}: estimated from {int(keep.sum())} undisturbed events: "
              + ", ".join(pop for pop, _ in POPS if pop in estimated))
    for pop, parent in POPS:
        rows.append((r.sample_id, pop, parent, c[pop], 100.0 * c[pop] / max(c[parent], 1)))
    M = matrix_for(r.fcs_file)
    for i, f in enumerate(fluors):
        spill_rows.append([r.sample_id, f] + [f"{v:.5f}" for v in M[i]])

pd.DataFrame(spill_rows, columns=["sample", "fluorochrome"] + dets).to_csv(out / "spillover.csv", index=False)
pd.DataFrame(rows, columns=["sample", "population", "parent", "count", "percent_of_parent"]).to_csv(
    out / "frequencies.tsv", sep="\t", index=False, float_format="%.3f")
pd.DataFrame(qc_rows, columns=["fmo_file", "detector", "voltage", "spread_before_fine_tuning"] +
             ["median_" + p for p in QC_POPS] + ["spread", "threshold_CD4_T", "threshold_CD8_T"]).to_csv(
    out / "qc_fmo.tsv", sep="\t", index=False, float_format="%.2f")
pd.DataFrame(qc_time, columns=["fcs_file", "events", "first_disturbed_event", "last_disturbed_event", "time_first",
                               "time_last", "fraction_of_events", "bins_in_run", "bins", "factor_affected",
                               "factor_other"]).to_csv(out / "qc_time.tsv", sep="\t", index=False)
print(f"wrote {len(rows)} rows for {len(samples)} samples, {len(mats)} detector setting(s), "
      f"{len(qc_time)} tube(s) with a disturbed stretch")
