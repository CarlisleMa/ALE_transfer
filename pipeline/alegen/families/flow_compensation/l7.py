"""Level 7 of `flow_compensation`: the level-4 session plus one silent decision point.

  T_acquisition_disturbance  In two of the six samples (one recorded before and one after the voltage change
                             of the HLA-DR detector) a contiguous part of the run, 12-25% of the events in the
                             middle of the acquisition as given by `Time`, was recorded during a fluidics
                             disturbance. In that window the two detectors of the 561-nm laser (YG586-A: CD4
                             PE, YG780-A: HLA-DR PE-Cy7) record a fraction of the signal (0.30-0.45, with an
                             event-to-event jitter of 1-2%); scatter and the detectors of the other three
                             lasers are as they should be, and the event rate is 3-10% lower.

What an analysis that gates all events gets wrong: the compensation subtracts too little PE signal from the CD8
detector (B710-A), so the CD4 T cells of the window are counted as CD4+CD8+; and it subtracts the full
PerCP-Cy5.5 spill from a YG780-A signal that was recorded too low, so no CD8 T cell of the window passes the
HLA-DR gate. Nothing fails. No text mentions it: the instruction, the SOP and the notes are the level-4 texts
(`python -m alegen make --level 7` then removes the scaffolding and the incident reports as for level 6).

Truth by construction: the reference is that of the same session without the disturbance. `generate_l4`
simulates the recording, computes the reference from it (hidden labels; SOP rule on planted quantities) and
writes to the file what the function `disturb` below returns. Membership of an event does not depend on
whether it was recorded inside the window. Every other input file is the level-4 file of the same seed.

The window is drawn from the seed, among positions for which the reference percentages of the events that
remain after (a) removing the window and (b) removing everything from its start to the end of the file are
all within `KEEP_BANDS[0]` of the tolerance band (relaxed stepwise if no sample of a block has such a
position; the generator raises beyond the last step). That makes "an unbiased estimate from the undisturbed
part is within tolerance" a property of the instance instead of a matter of sampling luck. The generator also
checks that the QC named in the task card (median per time bin) separates the window from the rest of the file.

Levels 1-6 do not import this module.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from ...spec import Attack, Probe, Trap
from . import ASSETS, DETS, NF, POPS, TOL_ABS, TOL_REL, iCD4, iDR
from .l4 import ACT, ACT_ABS, ACT_POPS, ACT_REL, FMO_PCT, PARENT4, POP_NAMES4, generate_l4

AFFECTED = [iCD4, iDR]                # detectors of the yellow-green laser
WINDOW_FRAC = (0.12, 0.25)            # share of the events of the file
WINDOW_START = (0.25, 0.60)           # start of the window as a fraction of the run; it ends before 0.85
ATTENUATION = (0.30, 0.45)            # fraction of the signal that is recorded inside the window
JITTER = (0.01, 0.02)                 # log-sd of the event-to-event variation of that fraction
RATE = (0.90, 0.97)                   # event rate inside the window relative to the rest of the run
KEEP_BANDS = (0.5, 0.6, 0.7)          # see module docstring; the first level that a sample of the block meets
CANDIDATES = 600
QC_BIN = 500                          # events per time bin of the QC that reveals the window

L7_JUDGE = """

**Level-7 decision point (acquisition disturbance).** In two samples a contiguous time window of the file
(12-25% of its events) was recorded with the two detectors of the 561-nm laser reading 0.30-0.45 of the
signal. The reference is that of the same simulated session without the disturbance: lineage membership is
the hidden label of every event, disturbed or not, and the FMO-gated populations are the SOP rule applied to
the signals as they would have been recorded. The scorer is the level-4 scorer, unchanged; the reference of
a level-7 task is identical to that of the level-4 task of the same seed.

Three treatments reach the reference and are probes of the battery: the window is set aside and the
populations that need the affected detectors are estimated from the other events; everything from the start
of the window to the end of the file is set aside; the affected detectors are rescaled inside the window.
The generator only plants a window if the reference percentages of the events that remain, in both of the
first two cases, are within half the tolerance band for all 14 populations, so the estimate from the
undisturbed part is within tolerance by construction and not by sampling luck. Measured on twelve seeds
(reference pipeline, worst item as a fraction of its band): window set aside 0.48 lineage / 0.41 FMO-gated,
everything after the start set aside 0.64 / 0.57, rescaled 0.02 / 0.37.

Without the step (the level-4 reference pipeline on these inputs) six of the 84 scored percentages leave the
band, the same six on twelve seeds. In each of the two samples CD4_T is 8-14 points low (4.4-8.0 bands) and
DP_T 8-14 points high (15.8-28.8 bands; reported 9-17% instead of 1-3%), because the CD4 T cells of the
window are counted as CD4+CD8+, and activated_CD8_T is 1.0-5.8 points low (1.3-4.1 bands). The lymphocyte
counts are unaffected. Score 0.472 on every seed; the lineage effect alone would give 0.739, the effect on
activated_CD8_T alone 0.733.

**Scored count.** The only scored count is that of `lymphocytes`. The population is defined by scatter, the
viability dye and CD14, none of which is recorded by the disturbed detectors, so it can be counted on all
events of the file (both solutions report it within 0.3% of the planted count). A submission that removes
the window and reports the count of the events it kept is 19-26% low in two samples, loses two thirds of
the `lymphocyte_counts` component and passes with 0.933 if everything else is right (probe
`time_gate_counts_of_kept_events`).
"""


def _membership(data, labels, S, fmo):
    """Reference membership of every event in the 14 SOP populations, as in generate_l4."""
    comp = lambda d: (d[:, 3:3 + NF].astype(float) @ np.linalg.inv(S))[:, iDR]          # noqa: E731
    fd, fl = fmo
    xf, x = comp(fd), comp(data)
    m = {name: np.isin(labels, classes) for name, _, classes in POPS}
    for name, _, c in ACT_POPS:
        m[name] = (labels == c) & (x > float(np.percentile(xf[fl == c], FMO_PCT)))
    return m


def binned_qc(data, i0, i1):
    """The QC named in the task card, on the events as written: median (and 90th percentile) of every
    fluorescence detector in consecutive bins of QC_BIN events ordered by Time. Per detector and statistic:
    level outside and inside the planted window (median over the bins that lie wholly outside / inside), their
    ratio, the step in units of the bin-to-bin scatter outside the window (1.4826 x MAD of the outside bins),
    and how many inside / outside bins lie beyond the midpoint of the two levels."""
    n = len(data)
    nb = n // QC_BIN
    edges = np.linspace(0, n, nb + 1).astype(int)
    inside = (edges[:-1] >= i0) & (edges[1:] <= i1)
    outside = (edges[1:] <= i0) | (edges[:-1] >= i1)
    out = {}
    for name, q in (("median", 50), ("p90", 90)):
        v = np.array([np.percentile(data[edges[b]:edges[b + 1], 3:3 + NF].astype(float), q, axis=0)
                      for b in range(nb)])
        for j in range(NF):
            lo, hi = float(np.median(v[inside, j])), float(np.median(v[outside, j]))
            mid = float(np.sqrt(lo * hi)) if lo > 0 and hi > 0 else 0.5 * (lo + hi)
            beyond = (v[:, j] < mid) if lo < hi else (v[:, j] > mid)
            scatter = 1.4826 * float(np.median(np.abs(v[outside, j] - hi)))
            out.setdefault(DETS[j], {})[name] = {
                "outside": hi, "inside": lo, "ratio": lo / hi if hi else float("nan"),
                "step_in_robust_sd": (lo - hi) / scatter if scatter > 0 else float("nan"),
                "inside_bins_beyond_midpoint": int((beyond & inside).sum()),
                "outside_bins_beyond_midpoint": int((beyond & outside).sum())}
    return out, int(nb), int(inside.sum()), int(outside.sum())


def _find_window(drng, m, n):
    """First of CANDIDATES seeded window positions for which the reference percentages of the events that
    remain (window removed; everything from the window start removed) are within KEEP_BANDS[0] of the
    tolerance band for all 14 populations; otherwise the best one.
    -> (worst deviation in bands, first event, end event, candidates tried)"""
    cum = {p: np.concatenate([[0], np.cumsum(m[p])]) for p in POP_NAMES4}
    cum["all_events"] = np.arange(n + 1)

    def pct(i0, i1, p):                              # percent of parent among the events [0, i0) and [i1, n)
        c = cum[p][i0] + cum[p][n] - cum[p][i1]
        d = cum[PARENT4[p]][i0] + cum[PARENT4[p]][n] - cum[PARENT4[p]][i1]
        return 100.0 * c / max(d, 1)

    ref = {p: pct(0, 0, p) for p in POP_NAMES4}
    band = {p: max(ACT_ABS, ACT_REL * ref[p]) if p in ACT else max(TOL_ABS, TOL_REL * ref[p]) for p in POP_NAMES4}
    best = None
    for tried in range(1, CANDIDATES + 1):
        frac = float(drng.uniform(*WINDOW_FRAC))
        start = float(drng.uniform(WINDOW_START[0], min(WINDOW_START[1], 0.85 - frac)))
        i0 = int(round(start * n))
        i1 = i0 + int(round(frac * n))
        w = max(max(abs(pct(i0, i1, p) - ref[p]), abs(pct(i0, n, p) - ref[p])) / band[p] for p in POP_NAMES4)
        if best is None or w < best[0]:
            best = (float(w), i0, i1, tried)
        if w <= KEEP_BANDS[0]:
            break
    return best


def generate_l7(seed: int, level: int, workdir: Path):
    if level != 7:
        raise ValueError("l7.generate_l7 builds level 7 only")
    drng = np.random.default_rng(71000 + seed)       # own stream: the level-4 session of the seed is unchanged

    # ---- pass 1: simulate the session as at level 4 and look for a usable window in every sample ---------
    found: dict[str, tuple] = {}
    block_of: dict[str, int] = {}

    def survey(sid, block, blocks, data, labels, S, fmo):
        block_of.update(blocks)
        found[sid] = _find_window(drng, _membership(data, labels, S, fmo), len(data))
        return data

    generate_l4(seed, level, workdir, disturb=survey)
    chosen: dict[str, dict] = {}
    for bl in (1, 2):
        members = [s for s in found if block_of[s] == bl]
        usable, keep_band = [], None
        for keep_band in KEEP_BANDS:
            usable = [s for s in members if found[s][0] <= keep_band]
            if usable:
                break
        if not usable:
            raise RuntimeError(f"level 7, seed {seed}: no sample of block {bl} has a window that keeps the remaining "
                               f"events within {KEEP_BANDS[-1]} band (best {min(found[s][0] for s in members):.2f})")
        sid = usable[int(drng.integers(0, len(usable)))]
        chosen[sid] = {"attenuation": float(drng.uniform(*ATTENUATION)), "jitter": float(drng.uniform(*JITTER)),
                       "rate": float(drng.uniform(*RATE)), "z": drng.standard_normal(found[sid][2] - found[sid][1]),
                       "keep_band": keep_band, "usable": len(usable), "members": len(members)}

    # ---- pass 2: the same session; the two chosen files are written with the disturbance -----------------
    plants: dict[str, dict] = {}

    def disturb(sid, block, blocks, data, labels, S, fmo):
        if sid not in chosen:
            return data
        c = chosen[sid]
        w, i0, i1, tried = found[sid]
        n = len(data)
        out = data.astype(np.float64)
        f = c["attenuation"] * np.exp(c["jitter"] * c["z"])          # one factor per event, both detectors
        for j in AFFECTED:
            out[i0:i1, 3 + j] *= f
        gaps = np.diff(out[:, -1], prepend=0.0)
        gaps[i0:i1] /= c["rate"]                                     # events arrive a little more slowly
        out[:, -1] = np.cumsum(gaps)
        out = out.astype(np.float32)
        if np.any(np.diff(out[:, -1]) < 0):
            raise RuntimeError("level 7: Time is not monotonic after the disturbance")
        qc, nb, nb_in, nb_out = binned_qc(out, i0, i1)
        clean_step = [DETS[j] for j in AFFECTED
                      if qc[DETS[j]]["median"]["inside_bins_beyond_midpoint"] == nb_in
                      and qc[DETS[j]]["median"]["outside_bins_beyond_midpoint"] == 0]
        if not clean_step:
            raise RuntimeError(f"level 7, seed {seed}: the median per time bin does not separate the window of {sid} "
                               "from the rest of the file in either affected detector")
        plants[sid] = {"block": int(block), "n_events": int(n), "first_event": int(i0), "end_event": int(i1),
                       "fraction_of_events": round((i1 - i0) / n, 4), "start_fraction": round(i0 / n, 4),
                       "time_first": float(out[i0, -1]), "time_last": float(out[i1 - 1, -1]),
                       "detectors": [DETS[j] for j in AFFECTED], "attenuation": c["attenuation"],
                       "jitter_log_sd": c["jitter"], "event_rate_factor": c["rate"],
                       "candidates_tried": int(tried), "kept_events_worst_band": round(w, 3),
                       "kept_events_band_required": c["keep_band"], "usable_samples_in_block": c["usable"],
                       "samples_in_block": c["members"], "qc_bin_events": QC_BIN, "qc_bins": nb,
                       "qc_bins_inside": nb_in, "qc_bins_outside": nb_out,
                       "qc_detectors_with_clean_median_step": clean_step, "qc": qc}
        return out

    b = generate_l4(seed, level, workdir, disturb=disturb)
    if len(plants) != 2 or sorted(p["block"] for p in plants.values()) != [1, 2]:
        raise RuntimeError("level 7: expected one disturbed sample per voltage block")
    truth = json.loads(b.tests_files["reference/truth.json"])
    truth["plants"]["acquisition_disturbance"] = plants
    b.tests_files["reference/truth.json"] = json.dumps(truth)

    # ---- solutions ---------------------------------------------------------------------------------------
    b.solution_files = {"solve.sh": ASSETS / "solve_l7.sh", "flow_pipeline_l7.py": ASSETS / "flow_pipeline_l7.py",
                        "solve_alt_l7.py": ASSETS / "solve_alt_l7.py"}

    # ---- battery: the level-4 battery on the level-7 inputs, plus the level-7 probes and attacks ---------
    run4, run7 = "python3 /alegen_script/flow_pipeline_l4.py ", "python3 /alegen_script/flow_pipeline_l7.py "
    l4_py = {"flow_pipeline_l4.py": ASSETS / "flow_pipeline_l4.py"}
    l7_py = {"flow_pipeline_l7.py": ASSETS / "flow_pipeline_l7.py"}
    ceil = b.pass_threshold - 0.15
    attacks = [
        Attack("l6_reference", "the reference pipeline of levels 4-6, unchanged, on the level-7 inputs: every "
               "decision point of the lower levels is handled, the files are not checked along Time and all events "
               "are gated", ceil, script=run4.strip(), script_files=l4_py),
        Attack("l6_second_solution", "the independent second solution of levels 4-6, unchanged, on the level-7 "
               "inputs (mixture-model gates instead of density valleys, no check along Time)", ceil,
               script="python3 /alegen_script/solve_alt_l4.py",
               script_files={"solve_alt_l4.py": ASSETS / "solve_alt_l4.py"}),
        Attack("exclude_wrong_window", "T_acquisition_disturbance mis-handled: the level-7 reference pipeline, but "
               "the time gate (right length) is set on another part of the file: the start of the run if there is "
               "more room before the disturbed stretch than after it, otherwise the end", ceil,
               script=run7 + "--time-qc wrong", script_files=l7_py),
    ]
    for a in b.attacks:
        if a.script is None:
            attacks.append(a)                        # truth-derived files: the reference is that of level 4
        elif a.name.startswith("l3_reference"):      # lower-level reference methods run as they were written
            attacks.append(Attack(a.name, a.description.replace("level-4 inputs", "level-7 inputs"), a.max_score,
                                  script=a.script, script_files=a.script_files))
        else:                                        # "ignores exactly one element": the level-7 step is done
            if not a.script.startswith(run4):
                raise RuntimeError(f"level 7: unexpected script of the level-4 attack {a.name}")
            attacks.append(Attack(a.name, a.description + " (signal-versus-time check done, disturbed stretch set "
                                  "aside)", a.max_score, script=run7 + a.script[len(run4):], script_files=l7_py))
    b.attacks = attacks
    probes = []
    for p in b.probes:
        if p.name == "independent_solution":
            p = Probe(p.name, "second implementation: disturbed stretch located by an event-level rank scan (no bins, "
                      "no medians), cut out of the tube, counts scaled back to the file; fcsparser reader, regression "
                      "spillover, gain factor from mixture-model bead peaks, mixture-model gates", p.min_score,
                      script=p.script)
        probes.append(p)
    probes += [
        Probe("time_gate_everything_after", "reference pipeline, but everything from the start of the disturbed "
              "stretch to the end of the file is set aside (35-60% of the events remain for the populations that "
              "need the affected detectors)", 0.99, script=run7 + "--time-qc after", script_files=l7_py),
        Probe("repair_instead_of_exclude", "reference pipeline, but the stretch is not set aside: inside it the "
              "affected detectors are divided by their factor (position of the peaks of the log signal inside "
              "against outside), then all events are gated", 0.99, script=run7 + "--time-qc repair",
              script_files=l7_py),
        Probe("time_gate_counts_of_kept_events", "reference pipeline, but the disturbed stretch is removed from the "
              "file and the sizes of what is left are reported. Floor = pass threshold instead of 0.99: the "
              "percentages are right, the lymphocyte count of the two samples is low by the share of the window, "
              "which costs two thirds of the lymphocyte_counts component (0.933)", b.pass_threshold,
              script=run7 + "--time-qc drop", script_files=l7_py),
    ]
    b.probes = probes

    # ---- task card ---------------------------------------------------------------------------------------
    def where(sid):
        p = plants[sid]
        return (f"{sid} (events {p['first_event']}-{p['end_event'] - 1} of {p['n_events']}, "
                f"{100 * p['fraction_of_events']:.0f}% of the file, Time {p['time_first']:.0f}-{p['time_last']:.0f}, "
                f"signal x{p['attenuation']:.2f}, event rate x{p['event_rate_factor']:.2f})")

    def evidence(sid):
        p = plants[sid]
        q = p["qc"]
        lev = lambda d, s: (f"{d} {q[d][s]['outside']:.0f} -> {q[d][s]['inside']:.0f} (x{q[d][s]['ratio']:.2f}, "  # noqa: E731
                            f"{abs(q[d][s]['step_in_robust_sd']):.0f} SD)")
        stray = lambda d: (f"{q[d]['median']['outside_bins_beyond_midpoint']} of {p['qc_bins_outside']}")  # noqa: E731
        others = max(abs(q[d]["median"]["step_in_robust_sd"]) for j, d in enumerate(DETS) if j not in AFFECTED)
        ya, yb = DETS[iDR], DETS[iCD4]
        return (f"{sid}: bin medians {lev(ya, 'median')} and {lev(yb, 'median')} in {p['qc_bins_inside']} "
                f"consecutive bins of {p['qc_bins']} (SD = robust standard deviation of the bin medians outside the "
                f"window); bins outside the window whose median is nearer to the level inside: {stray(ya)} ({ya}), "
                f"{stray(yb)} ({yb}); 90th percentile per bin {lev(yb, 'p90')}; the bin medians of the other six "
                f"detectors differ between inside and outside by at most {others:.1f} SD")

    sids = sorted(plants, key=lambda s: plants[s]["block"])
    b.traps = list(b.traps) + [Trap(
        "T_acquisition_disturbance",
        "A contiguous part of two runs was recorded during a fluidics disturbance: " + " and ".join(
            where(s) for s in sids) + ". Inside the window the two detectors of the 561-nm laser, YG586-A (CD4 PE) "
        "and YG780-A (HLA-DR PE-Cy7), read a fraction of the signal (event-to-event jitter 1-2%); scatter and the "
        "other six detectors are unaffected. Gated with all events, the PE spill into B710-A is under-compensated "
        "inside the window, so its CD4 T cells are counted as CD4+CD8+ (CD4_T low and DP_T high by about the share "
        "of the window), and the full PerCP-Cy5.5 spill is subtracted from a YG780-A signal that was recorded too "
        "low, so no CD8 T cell of the window passes the HLA-DR gate (activated_CD8_T low by the share of the "
        "window). Nothing fails. No text mentions time, stability or the disturbance; `Time` is one of the "
        "recorded parameters.",
        "signal-versus-Time QC as done by flowAI, PeacoQC or a manual Time gate: median of every fluorescence "
        f"detector in consecutive bins of {QC_BIN} events ordered by Time. " + ". ".join(evidence(s) for s in sids)
        + ". The same window shows in any CD4-versus-CD8 plot of the T cells of these samples as a third population "
        "(CD4 lower by the factor of the window, CD8 intermediate) that exists only in that Time range.",
        "Approach / skipped QC step")]
    b.judge_explanation_md = b.judge_explanation_md.rstrip() + L7_JUDGE
    b.capabilities = list(b.capabilities) + ["acquisition QC along the Time parameter"]
    b.oracle_steps, b.oracle_runtime_min = 21, 0.1
    return b
