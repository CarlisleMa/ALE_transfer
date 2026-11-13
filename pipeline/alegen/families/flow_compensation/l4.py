"""Level 4 of `flow_compensation`: compensation that has to be right, not just done (P + Q).

The level-3 task (identity/odd keyword -> matrix from the controls, aggregates, dead cells, one
re-acquired sample with a lowered FSC voltage) plus four elements taken from the small-model pilot
(README section 6: every QC-discoverable decision point was handled in 2 min 40 s):

  P_embedded_matrix  every file carries a well-formed, non-identity $SPILLOVER from the experiment
                     template of an older panel set-up; it is wrong for today's detectors.
  P_voltage_change   the PMT voltage of the HLA-DR detector (YG780-A) was raised part-way through the
                     run. Only the PE-Cy7 bead control and the HLA-DR FMO tube were recorded again.
                     Samples before/after need different matrices; the new control gives the PE-Cy7
                     row, but every other fluorochrome's spill INTO that detector scaled as well.
  P_tandem_lot       the PerCP-Cy5.5 bead control was stained from another lot than the sample
                     cocktail and spills less into YG780-A; a PBMC single stain from the cocktail lot
                     is provided.
  Q_fmo_gate         two populations without a valley (HLA-DR+ CD4 / CD8 T cells, dim continuous
                     signal) whose gate the SOP defines on the FMO control (99th percentile, same
                     parent, compensated data). The FMO tubes are also the in-data control of the
                     compensation: the SOP states an acceptance criterion on them (medians of five
                     lymphocyte populations within 50 units). A compensation that meets it lands in
                     the scored band; every careless matrix misses it by a factor 16 or more, and
                     can be repaired on the FMO tube itself.

Truth by construction: the per-sample spillover matrix is the planted one (sample lot, detector
settings of the sample's block). Lineage populations are counted from the hidden class labels. The
two FMO-gated populations are operational by definition; their reference is the SOP rule applied
to the planted quantities (signals compensated with the planted matrix, parent = hidden label,
threshold = 99th percentile of the label-defined parent in the FMO tube of the same block). No
gating or matrix estimate of the reference pipeline enters the reference.

Kept in a separate module so that levels 1-3 stay byte-identical.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from ... import dockerlib
from ...spec import Attack, Probe, TaskBuild, Trap
from . import (AF_BEAD, ASSETS, BC, BEAD, CD4T, CD8T, DEAD, DEBRIS, DETS, DNT, DOUBLET, DPT, FLUORS, FSC_V,
               MARKERS, MONO, NF, NK, OTHER, PANEL, POPS, SSC_V, TOL_ABS, TOL_COUNT, TOL_REL, _assemble,
               _cell_control, _planted_spillover, _singles, iCD8, iCD19, iCD56, iDR, iVIA, write_fcs)

DR_DIM, DR_SIG = 700.0, 0.75          # HLA-DR PE-Cy7 on activated T cells: median amount, log-sd
FMO_PCT = 99.0                        # percentile of the FMO control that defines HLA-DR +
QC_UNITS = 50                         # SOP acceptance criterion on the FMO control (detector units)
ACT_POPS = [("activated_CD4_T", "CD4_T", CD4T), ("activated_CD8_T", "CD8_T", CD8T)]
POPS4 = [(n, p) for n, p, _ in POPS] + [(n, p) for n, p, _ in ACT_POPS]
POP_NAMES4 = [p[0] for p in POPS4]
PARENT4 = dict(POPS4)
LINEAGE = [p[0] for p in POPS]
ACT = [p[0] for p in ACT_POPS]
# tolerances of the two FMO-gated populations (percentage points / fraction of the reference value)
ACT_ABS, ACT_REL = 0.75, 0.06


# --------------------------------------------------------------------------------------- simulator
def _detect_gain(rng, T, af, S, q, en, g):
    """Detector model of levels 1-3 with the gain of the HLA-DR detector multiplied by g: optical
    signal and photon noise scale with the gain, electronic noise does not."""
    mean = T @ S + af
    x = mean + np.sqrt(np.maximum(mean, 0.0) * q) * rng.standard_normal(mean.shape)
    x[:, iDR] *= g
    return x + en * rng.standard_normal(mean.shape)


def _singles4(rng, cls, stain, a4, a8, fmo):
    T, fsc, ssc, af = _singles(rng, cls, stain)
    is_t = (cls >= CD4T) & (cls <= DNT)
    T[is_t, iDR] = 0.0                # level 4 replaces the bright HLA-DR+ T subset of levels 1-3
    p = np.where((cls == CD8T) | (cls == DPT), a8, a4)
    act = is_t & (rng.random(len(cls)) < p)
    T[act, iDR] = DR_DIM * stain[iDR] * np.exp(DR_SIG * rng.standard_normal(int(act.sum())))
    if fmo:
        T[:, iDR] = 0.0               # no HLA-DR PE-Cy7 in the tube: no specific and no sticky signal
    return T, fsc, ssc, af, act


def _sample4(rng, S, q, en, g, fmo=False):
    """A level-3 style PBMC tube (debris, aggregates, dead cells) with dim HLA-DR on a seed-dependent
    fraction of the T cells. Returns scatter, detector signals, class labels, latent activation."""
    n = int(rng.integers(36000, 44001))
    n_deb = int(round(n * rng.uniform(0.06, 0.15)))
    n_cells = n - n_deb
    n_dbl = int(round(n_cells * rng.uniform(0.06, 0.12)))
    n_single = n_cells - n_dbl
    n_dead = int(round(n_single * rng.uniform(0.08, 0.18)))
    n_live = n_single - n_dead
    n_mono = int(round(n_live * rng.uniform(0.10, 0.25)))
    n_lym = n_live - n_mono
    p_l = np.array([rng.uniform(0.55, 0.75), rng.uniform(0.06, 0.15), rng.uniform(0.06, 0.18),
                    rng.uniform(0.015, 0.04)])
    n_t, n_b, n_nk, n_other = rng.multinomial(n_lym, p_l / p_l.sum())
    p_t = np.array([rng.uniform(0.50, 0.68), rng.uniform(0.22, 0.38), rng.uniform(0.008, 0.025),
                    rng.uniform(0.03, 0.07)])
    n4, n8, ndp, ndn = rng.multinomial(n_t, p_t / p_t.sum())
    a4, a8 = float(rng.uniform(0.03, 0.12)), float(rng.uniform(0.10, 0.35))
    counts = {DEBRIS: n_deb, DEAD: n_dead, CD4T: n4, CD8T: n8, DPT: ndp, DNT: ndn, BC: n_b, NK: n_nk,
              OTHER: n_other, MONO: n_mono}
    stain = np.exp(0.08 * rng.standard_normal(NF))
    cls = np.concatenate([np.full(c, k) for k, c in counts.items()])
    T, fsc, ssc, af, act = _singles4(rng, cls, stain, a4, a8, fmo)
    fsch = fsc * 0.92 * (1 + 0.02 * rng.standard_normal(len(cls)))
    cell_classes = [k for k in counts if k != DEBRIS]
    p = np.array([counts[k] for k in cell_classes], dtype=float)
    a = rng.choice(cell_classes, size=n_dbl, p=p / p.sum())
    b = rng.choice(cell_classes, size=n_dbl, p=p / p.sum())
    Ta, fa, sa, afa, _ = _singles4(rng, a, stain, a4, a8, fmo)
    Tb, fb, sb, afb, _ = _singles4(rng, b, stain, a4, a8, fmo)
    big, small = np.maximum(fa, fb), np.minimum(fa, fb)
    fd = big + small * rng.uniform(0.8, 1.0, n_dbl)
    fhd = big * 0.92 * rng.uniform(1.0, 1.05, n_dbl) * (1 + 0.02 * rng.standard_normal(n_dbl))
    T = np.vstack([T, Ta + Tb])
    af = np.vstack([af, afa + afb])
    fsc = np.concatenate([fsc, fd])
    fsch = np.concatenate([fsch, fhd])
    ssc = np.concatenate([ssc, sa + sb])
    cls = np.concatenate([cls, np.full(n_dbl, DOUBLET)])
    act = np.concatenate([act, np.zeros(n_dbl, bool)])
    fl = _detect_gain(rng, T, af, S, q, en, g)
    return fsc, fsch, ssc, fl, cls, act


def _bead_control4(rng, f, S, q, en, g, n=6000):
    T = np.zeros((n, NF))
    posb = rng.random(n) < 0.5
    T[posb, f] = BEAD[f] * np.exp(0.12 * rng.standard_normal(int(posb.sum())))
    af = np.exp(0.15 * rng.standard_normal(n))[:, None] * AF_BEAD[None, :] * np.exp(0.1 * rng.standard_normal((n, NF)))
    fsc = rng.normal(70000, 2500, n)
    ssc = 60000 * np.exp(0.08 * rng.standard_normal(n))
    fsch = fsc * 0.92 * (1 + 0.02 * rng.standard_normal(n))
    dbl = rng.random(n) < 0.03
    T[dbl] *= 2.0
    af[dbl] *= 2.0
    fsc[dbl] *= 1.9
    ssc[dbl] *= 2.0
    return _assemble(rng, fsc, fsch, ssc, _detect_gain(rng, T, af, S, q, en, g))[0]


def _cell_single_stain(rng, f, S, q, en, n=30000):
    """PBMC aliquot stained with one antibody conjugate only (from the vial used for the cocktail)."""
    cls = rng.choice([DEBRIS, DEAD, CD4T, CD8T, DPT, DNT, BC, NK, OTHER, MONO], size=n,
                     p=[0.08, 0.04, 0.36, 0.17, 0.01, 0.03, 0.07, 0.08, 0.02, 0.14])
    T, fsc, ssc, af = _singles(rng, cls, np.ones(NF))
    keep = T[:, f].copy()
    T[:] = 0.0
    T[:, f] = keep
    fsch = fsc * 0.92 * (1 + 0.02 * rng.standard_normal(n))
    return _assemble(rng, fsc, fsch, ssc, _detect_gain(rng, T, af, S, q, en, 1.0))[0]


# --------------------------------------------------------------------------------------- texts
INSTRUCTION = """\
You are the cytometry data analyst of an immunology core facility. A study team acquired PBMC
samples with an 8-colour immunophenotyping panel and needs, for every sample, the spillover matrix
that applies to it and the size of every population in the lab's gating SOP.

Task directory:
- `/task`

Inputs (read-only) under `/task/input`:
- `fcs/*.fcs` - FCS 3.1 list-mode files as exported by the acquisition software: {n_samples} PBMC
  samples and the {n_ctrl} control tubes of the session
- `samples.tsv` - sample sheet (`sample_id`, `fcs_file`, `donor`, `acquisition_order`)
- `controls.tsv` - control tubes (`fcs_file`, `control_type`, `fluorochrome`, `detector`, `carrier`,
  `reagent_lot`). For a single-stain control the columns name the one reagent in the tube; for an
  FMO control they name the one reagent of the panel that was left out.
- `panel.tsv` - staining panel of the samples (`detector`, `fluorochrome`, `marker`, `reagent_lot`)
- `gating_sop.md` - the lab's gating SOP: population names, parent of each population, definitions,
  the control-based HLA-DR gate and the compensation acceptance check
- `acquisition_notes.md` - the operator's notes from the acquisition session

Available software: Python 3.12 (on `PATH`) with numpy, scipy, pandas, scikit-learn, FlowKit,
FlowIO and fcsparser. There is no internet access.

Your task:
1. For every sample in `samples.tsv`, determine the spillover matrix that describes the signals in
   its file.
2. For every sample in `samples.tsv`, identify every population of `gating_sop.md` as defined there
   and report its size.

Required outputs under `/task/output`:
- `spillover.csv` - comma-separated. Header: `sample`, `fluorochrome`, followed by the eight
  fluorescence detectors named as in `panel.tsv`. Exactly one row for every `sample_id` of
  `samples.tsv` and every fluorochrome of `panel.tsv` (names as written there). The eight rows of a
  sample are the spillover matrix of that sample: the entry in row *f*, column *d* is the signal
  that fluorochrome *f*, as it is present on the cells of that sample, produces in detector *d* at
  the detector settings with which that sample was recorded, as a fraction of the signal it
  produces in its own primary detector at the same settings (background-free; the entry of each
  fluorochrome in its own detector is 1).
- `frequencies.tsv` - tab-separated with the columns `sample`, `population`, `parent`, `count`,
  `percent_of_parent`; exactly one row for every `sample_id` of `samples.tsv` and every population
  of `gating_sop.md` (names exactly as written there). Extra columns are allowed.
  - `parent`: the parent population named in the SOP.
  - `count`: number of events of the sample's FCS file that belong to the population (integer).
  - `percent_of_parent`: 100 x `count` / (count of the parent population). For the parent
    `all_events` the denominator is the total number of events in the FCS file.

Constraints:
- Do not modify anything under `/task/input`; write deliverables only under `/task/output`
  (additional files there are ignored).
- Every sample in the sample sheet belongs to the study and must be reported.
- Do not ask for confirmation. Execute directly.
"""

SOP = """\
# Gating SOP IMM-FC-012 (PBMC lineage panel with T-cell activation read-out, 8 colours)

Marker names are those of `panel.tsv`. "+" means that the cell carries the marker, i.e. its signal
from that reagent is clearly above the signal of cells that do not carry it; "-" means that it does not.
Populations are nested: an event can only belong to a population if it belongs to the parent.

| population | parent | definition |
|---|---|---|
| cells | all_events | intact cells of any kind; sub-cellular debris is excluded |
| singlets | cells | events that consist of exactly one cell |
| live | singlets | Viability dye - (membrane-intact cells) |
| lymphocytes | live | mononuclear cells with low side scatter, CD14 - |
| monocytes | live | CD14 + |
| T_cells | lymphocytes | CD3 + |
| CD4_T | T_cells | CD4 + CD8 - |
| CD8_T | T_cells | CD4 - CD8 + |
| DP_T | T_cells | CD4 + CD8 + |
| DN_T | T_cells | CD4 - CD8 - |
| B_cells | lymphocytes | CD3 - CD19 + |
| NK_cells | lymphocytes | CD3 - CD19 - CD56 + |
| activated_CD4_T | CD4_T | HLA-DR + by the FMO gate of section 2 |
| activated_CD8_T | CD8_T | HLA-DR + by the FMO gate of section 2 |

`all_events` denotes every event recorded in the file.

## 1. Controls of the HLA-DR read-out

HLA-DR on T cells is a dim, continuous signal without a separate positive peak. Its gate is
therefore not set on the stained sample but on a fluorescence-minus-one control ("FMO HLA-DR"):
PBMC of the laboratory's reference pool, stained with the complete panel except the HLA-DR reagent
and processed like a sample.

## 2. FMO gate for HLA-DR + (populations activated_CD4_T and activated_CD8_T)

1. The FMO tube is analysed exactly like a sample: same compensation as the samples that were
   recorded with the same fluorescence-detector voltages, same gates down to CD4_T and CD8_T.
2. For each of the two parent populations P (CD4_T, CD8_T), the threshold is the {pct:g}th percentile of
   the compensated signal of the HLA-DR detector (linear scale, detector units) over the events of
   P in the FMO tube.
3. An event of P in a sample is HLA-DR + if its compensated HLA-DR-detector signal is greater than
   the threshold of P from the FMO tube that was recorded with the same fluorescence-detector
   voltages as the sample.

## 3. Compensation acceptance for the HLA-DR read-out

The FMO tube contains no HLA-DR reagent, so whatever remains in the HLA-DR detector after
compensation is background, and lymphocyte populations that differ only in the other reagents they
carry must not differ in it. The compensation used for a group of samples is accepted for
reporting if, in the FMO tube recorded with the same fluorescence-detector voltages and compensated
with the same matrix, the median compensated HLA-DR-detector signals of CD4_T, CD8_T, DN_T, B_cells
and NK_cells lie within {qc} detector units of each other (largest minus smallest of the five
medians). A compensation that does not meet this criterion has to be corrected before the HLA-DR
populations are reported.
"""

JUDGE = """\
**Where the truth comes from.** Every event is simulated from a hidden class label (debris, aggregate
of two cells, dead cell, CD4/CD8/DP/DN T cell, B cell, NK cell, other lymphocyte, monocyte). The
reference `truth.json` holds, per sample, the planted spillover matrix of that sample (sample reagent
lots, detector settings of the sample's acquisition block; the block-2 matrix is the block-1 matrix
with the HLA-DR detector column multiplied and the PE-Cy7 row divided by the planted gain factor) and
the number of events of each lineage population counted from the labels. The two FMO-gated
populations are defined operationally by the SOP; their reference is that rule applied to planted
quantities only: signals compensated with the planted matrix, parent population = hidden label,
threshold = 99th percentile of the label-defined parent in the FMO tube of the same block. No gate or
matrix estimate of the reference pipeline is used as reference.

**Hard gates (score 0).** Missing/unparseable file; `spillover.csv` without the columns `sample`,
`fluorochrome` and the eight detectors, not exactly one row per sample x fluorochrome, non-numeric
entries; `frequencies.tsv` not holding exactly one row per sample x SOP population, wrong `parent`,
`count` not a non-negative integer, `percent_of_parent` outside [0, 100].

**Components.**

| weight | component | metric | full credit | zero |
|---|---|---|---|---|
| 0.20 | spillover_matrices | maximum absolute error over the 8 x 8 entries of every sample's matrix | <= 0.02 | >= 0.05 |
| 0.25 | lineage_every_sample | per sample, fraction of the 12 lineage populations whose `percent_of_parent` is within max(0.5 points, 3% of the reference); **minimum** over samples | 1.0 | <= 10/12 |
| 0.10 | lineage_overall | fraction of all sample x lineage-population items within the same tolerance | 1.0 | <= 0.5 |
| 0.10 | lymphocyte_counts | fraction of samples whose `lymphocytes` count is within 3% of the planted count | 1.0 | <= 0.5 |
| 0.25 | activation_every_sample | per sample, fraction of the two FMO-gated populations within max(0.75 points, 6% of the reference); **minimum** over samples | 1.0 | <= 0.5 |
| 0.10 | activation_overall | fraction of all sample x FMO-gated-population items within the same tolerance | 1.0 | 0 |

**Why these tolerances.** Lineage populations and matrix: as at level 3 (two unrelated correct
strategies agree with the planted values to a small fraction of the band). FMO-gated populations: the
band is tied to the SOP's own acceptance criterion. In a sweep of deliberately mis-set coefficients on
eight seeds, a compensation whose FMO medians agreed within the stated 50 units stayed within 0.58 of
the band (25 block results), between 50 and 100 units one of 47 block results left the band, and each
careless matrix misses the criterion by a factor 16 or more (the residual signal widens the FMO
background of CD8 T cells and loses a large part of the dim positives). Measured values are in
docs/families/flow_compensation.md ("Level 4").

**What is deliberately not scored.** Counts of populations other than `lymphocytes`, the SOP's
acceptance check itself (it is the analyst's own control; only its consequences are scored), any
plot, any additional file.
"""


def _spill_long(sids, mats, fmt="{:.5f}") -> str:
    out = ["sample,fluorochrome," + ",".join(DETS)]
    for s in sids:
        M = mats[s]
        for i in range(NF):
            out.append(s + "," + FLUORS[i] + "," + ",".join(fmt.format(M[i, j]) for j in range(NF)))
    return "\n".join(out) + "\n"


def _freq_tsv(sids, table) -> str:
    out = ["sample\tpopulation\tparent\tcount\tpercent_of_parent"]
    for s in sids:
        for pop in POP_NAMES4:
            c, pct = table[s][pop]
            out.append(f"{s}\t{pop}\t{PARENT4[pop]}\t{int(c)}\t{pct:.3f}")
    return "\n".join(out) + "\n"


def _spill_keyword(M: np.ndarray) -> str:
    """$SPILLOVER value: n, detector names, n x n coefficients row by row (4 decimals, trailing zeros dropped)."""
    vals = [f"{M[i, j]:.4f}".rstrip("0").rstrip(".") or "0" for i in range(NF) for j in range(NF)]
    return f"{NF}," + ",".join(DETS) + "," + ",".join(vals)


def generate_l4(seed: int, level: int, workdir: Path, disturb=None) -> TaskBuild:
    """`disturb` is None at levels 4-6. Level 7 (l7.py) passes a function that receives the recording of
    every sample (events in time order, hidden labels, planted matrix, FMO tube of its block) and returns
    the events that are written to the file. It must not draw from `rng`. The reference is computed from
    the recording as simulated here, never from what the function returns."""
    rng = np.random.default_rng(64000 + seed)
    S = _planted_spillover(rng)                      # sample lots, detector settings of block 1
    q = rng.uniform(6.0, 12.0, size=NF)
    en = rng.uniform(20.0, 35.0, size=NF)
    n_samples = 6
    donors = sorted(int(x) for x in rng.permutation(np.arange(101, 990))[:n_samples])
    sids = [f"PBMC_D{d}" for d in donors]
    shift_idx = int(rng.integers(1, n_samples))
    fsc_gain = float(rng.uniform(0.20, 0.27))
    shift_v = int(round(FSC_V * fsc_gain ** (1 / 7.0)))
    day = int(rng.integers(3, 27))
    month = ['FEB', 'MAR', 'APR', 'MAY', 'JUN', 'SEP', 'OCT'][int(rng.integers(0, 7))]
    date = f"{day:02d}-{month}-2025"

    # ---- level-4 plants ------------------------------------------------------------------------
    n_block1 = int(rng.integers(2, 5))               # samples recorded before the voltage change
    g = float(rng.uniform(1.8, 2.6))                 # gain factor of the HLA-DR detector in block 2
    dr_v1 = PANEL[iDR][3]
    dr_v2 = int(round(dr_v1 * g ** (1 / float(rng.uniform(6.8, 7.9)))))
    S2 = S.copy()
    S2[:, iDR] *= g
    S2[iDR, :] /= g
    S2[iDR, iDR] = 1.0
    # PerCP-Cy5.5 lot of the bead control: less spill into the HLA-DR detector than the cocktail lot
    S_beadlot = S.copy()
    tandem_delta = float(rng.uniform(0.10, 0.14))
    S_beadlot[iCD8, iDR] = S[iCD8, iDR] - tandem_delta
    S_beadlot[iCD8, iCD56] = S[iCD8, iCD56] * rng.uniform(0.6, 0.9)
    S_beadlot[iCD8, iCD19] = S[iCD8, iCD19] * rng.uniform(0.6, 0.9)
    S_beadlot = np.round(S_beadlot, 4)
    # template matrix of the older panel set-up: other detector gains, hence other coefficients
    h = np.exp(rng.uniform(-0.35, 0.35, size=NF))
    h[iCD8] = rng.uniform(1.6, 2.4)
    h[iDR] = rng.uniform(0.7, 1.0)
    S_old = S * h[None, :] / h[:, None] * np.exp(0.12 * rng.standard_normal((NF, NF)))
    np.fill_diagonal(S_old, 1.0)
    S_old = np.round(S_old, 4)
    lots = {f: str(int(x)) for f, x in zip(FLUORS, rng.permutation(np.arange(2100000, 9900000))[:NF])}
    bead_lot_percp = str(int(lots[FLUORS[iCD8]]) - int(rng.integers(40000, 900000)))
    tmpl_year, svc_month = 2024, ['NOV', 'DEC', 'JAN'][int(rng.integers(0, 3))]
    tmpl_date = f"{int(rng.integers(3, 27)):02d}-{['MAR', 'APR', 'MAY'][int(rng.integers(0, 3))]}-{tmpl_year}"
    svc_date = f"{int(rng.integers(3, 27)):02d}-{svc_month}-{2025 if svc_month == 'JAN' else 2024}"

    fcs_dir = workdir / "fcs"
    if fcs_dir.exists():
        for f in fcs_dir.glob("*.fcs"):
            f.unlink()
    fcs_dir.mkdir(parents=True, exist_ok=True)
    spill_kw = _spill_keyword(S_old)
    clock = [9 * 60 + int(rng.integers(5, 40))]

    def emit(fname, data, tube, src, stained: list[int], fscv=FSC_V, drv=dr_v1):
        pns = ["", "", ""] + [MARKERS[i] if i in stained else "" for i in range(NF)] + [""]
        volts = [str(p[3]) for p in PANEL]
        volts[iDR] = str(drv)
        pnv = [str(fscv), str(fscv), str(SSC_V)] + volts + [""]
        start = clock[0]
        dur = max(1, int(round(float(data[-1, -1]) / 100 / 60)))
        clock[0] += dur + int(rng.integers(1, 4))
        kw = {"$CYT": "LSRFortessa X-20", "$CYTSN": "X65811-R0428", "$DATE": date,
              "$BTIM": f"{start // 60:02d}:{start % 60:02d}:{int(rng.integers(0, 60)):02d}",
              "$ETIM": f"{(start + dur) // 60:02d}:{(start + dur) % 60:02d}:{int(rng.integers(0, 60)):02d}",
              "$FIL": fname, "$INST": "Immunology Core", "$OP": "core-operator-2", "$SRC": src,
              "$SYS": "Windows 10", "$TIMESTEP": "0.01", "TUBE NAME": tube,
              "EXPERIMENT NAME": "IMM-FC-012_PBMC_lineage", "CREATOR": "BD FACSDiva Software Version 9.0",
              "FSC ASF": "0.78", "APPLY COMPENSATION": "TRUE", "$SPILLOVER": spill_kw}
        write_fcs(fcs_dir / fname, data, pns, pnv, kw)

    def tag(f):
        return FLUORS[f].replace(" ", "").replace(".", "").replace("-", "")

    # ---- tubes recorded before the samples -----------------------------------------------------
    controls = []   # fcs_file, control_type, fluorochrome, detector, carrier, reagent_lot
    data = _cell_control(rng, S, q, en, dye=False)
    emit("ctrl_unstained.fcs", data, "Unstained Control", "PBMC pool", [])
    controls.append(("ctrl_unstained.fcs", "unstained", "NA", "NA", "cells", "NA"))
    for f in range(NF):
        fname = f"ctrl_{tag(f)}.fcs"
        if f == iVIA:
            data = _cell_control(rng, S, q, en, dye=True)
            carrier, src, lot = "cells", "PBMC pool, 1:1 live / heat-killed", lots[FLUORS[f]]
        else:
            data = _bead_control4(rng, f, S_beadlot if f == iCD8 else S, q, en, 1.0)
            carrier, src = "beads", "antibody-capture beads"
            lot = bead_lot_percp if f == iCD8 else lots[FLUORS[f]]
        emit(fname, data, f"{FLUORS[f]} Stained Control", src, [f])
        controls.append((fname, "single_stain", FLUORS[f], DETS[f], carrier, lot))
    percp_cells = f"ctrl_{tag(iCD8)}_cells.fcs"
    emit(percp_cells, _cell_single_stain(rng, iCD8, S, q, en), f"{FLUORS[iCD8]} Stained Control (cells)",
         "PBMC pool", [iCD8])
    all_but_dr = [i for i in range(NF) if i != iDR]
    fmo = {}
    fsc, fsch, ssc, fl, cls, _ = _sample4(rng, S, q, en, 1.0, fmo=True)
    fmo[1] = _assemble(rng, fsc, fsch, ssc, fl, cls)
    emit("fmo_HLADR.fcs", fmo[1][0], "FMO HLA-DR", "PBMC pool", all_but_dr)

    # ---- samples; the HLA-DR detector voltage is raised after n_block1 samples -----------------
    block = {s: (1 if k < n_block1 else 2) for k, s in enumerate(sids)}
    S_of = {1: S, 2: S2}
    rerun_pecy7, rerun_fmo = f"ctrl_{tag(iDR)}_rerun.fcs", "fmo_HLADR_rerun.fcs"
    sample_data, latent = {}, {}
    n_events: dict[str, int] = {}
    truth_counts: dict[str, dict[str, int]] = {}
    for k, sid in enumerate(sids):
        if k == n_block1:
            emit(rerun_pecy7, _bead_control4(rng, iDR, S, q, en, g), f"{FLUORS[iDR]} Stained Control",
                 "antibody-capture beads", [iDR], drv=dr_v2)
            fsc, fsch, ssc, fl, cls, _ = _sample4(rng, S, q, en, g, fmo=True)
            fmo[2] = _assemble(rng, fsc, fsch, ssc, fl, cls)
            emit(rerun_fmo, fmo[2][0], "FMO HLA-DR", "PBMC pool", all_but_dr, drv=dr_v2)
        b = block[sid]
        fsc, fsch, ssc, fl, cls, act = _sample4(rng, S, q, en, 1.0 if b == 1 else g)
        gs = fsc_gain if k == shift_idx else 1.0
        n = len(cls)
        order_data, lab = _assemble(rng, fsc * gs, fsch * gs, ssc, np.column_stack([fl, act.astype(float)]), cls)
        act_lab = order_data[:, 3 + NF] > 0.5
        data = np.delete(order_data, 3 + NF, axis=1)
        written = data if disturb is None else disturb(sid=sid, block=b, blocks=block, data=data, labels=lab,
                                                       S=S_of[b], fmo=fmo[b])
        emit(f"{sid}.fcs", written, sid, f"donor {donors[k]}", list(range(NF)),
             fscv=shift_v if k == shift_idx else FSC_V, drv=dr_v1 if b == 1 else dr_v2)
        assert data.shape == (n, len(PANEL) + 4)
        sample_data[sid] = (data, lab)
        n_events[sid] = int(n)
        truth_counts[sid] = {name: int(np.isin(lab, classes).sum()) for name, _, classes in POPS}
        latent[sid] = {name: int((act_lab & (lab == c)).sum()) for name, _, c in ACT_POPS}
    # single-stain rows are listed in panel order (the session-start bead set first), then the extra tubes
    controls.append((percp_cells, "single_stain", FLUORS[iCD8], DETS[iCD8], "cells", lots[FLUORS[iCD8]]))
    controls.append((rerun_pecy7, "single_stain", FLUORS[iDR], DETS[iDR], "beads", lots[FLUORS[iDR]]))
    controls.append(("fmo_HLADR.fcs", "fmo", FLUORS[iDR], DETS[iDR], "cells", "NA"))
    controls.append((rerun_fmo, "fmo", FLUORS[iDR], DETS[iDR], "cells", "NA"))

    # ---- reference of the FMO-gated populations: the SOP rule on planted quantities ------------
    def comp_dr(data, M):
        return (data[:, 3:3 + NF].astype(float) @ np.linalg.inv(M))[:, iDR]

    fmo_thr, fmo_medians = {}, {}
    for b in (1, 2):
        d, lab = fmo[b]
        x = comp_dr(d, S_of[b])
        fmo_thr[b] = {name: float(np.percentile(x[lab == c], FMO_PCT)) for name, _, c in ACT_POPS}
        fmo_medians[b] = {nm: float(np.median(x[lab == c])) for nm, c in (
            ("CD4_T", CD4T), ("CD8_T", CD8T), ("DN_T", DNT), ("B_cells", BC), ("NK_cells", NK))}
    for sid in sids:
        d, lab = sample_data[sid]
        x = comp_dr(d, S_of[block[sid]])
        for name, _, c in ACT_POPS:
            truth_counts[sid][name] = int((x[lab == c] > fmo_thr[block[sid]][name]).sum())

    def pct_table(counts_of):
        tab = {}
        for s in sids:
            c = dict(counts_of[s])
            c["all_events"] = n_events[s]
            tab[s] = {p: (c[p], 100.0 * c[p] / max(c[PARENT4[p]], 1)) for p in POP_NAMES4}
        return tab

    truth_tab = pct_table(truth_counts)
    mats = {s: S_of[block[s]] for s in sids}
    truth = {"fluorochromes": FLUORS, "detectors": DETS, "samples": sids,
             "spillover": {s: mats[s].tolist() for s in sids},
             "populations": [{"name": n, "parent": p} for n, p in POPS4],
             "lineage": LINEAGE, "fmo_gated": ACT,
             "n_events": n_events, "counts": truth_counts,
             "percent_of_parent": {s: {p: truth_tab[s][p][1] for p in POP_NAMES4} for s in sids},
             "tolerance": {"percent_abs": TOL_ABS, "percent_rel": TOL_REL, "count_rel": TOL_COUNT,
                           "fmo_abs": ACT_ABS, "fmo_rel": ACT_REL},
             # not used for scoring; kept for the design notes
             "plants": {"block": block, "gain": g, "tandem_delta": tandem_delta, "fmo_threshold": fmo_thr,
                        "fmo_medians_planted_matrix": fmo_medians, "latent_activated": latent,
                        "embedded": S_old.tolist(), "bead_lot": S_beadlot.tolist(), "fsc_shift_sample": sids[shift_idx]}}

    # ---- text inputs ---------------------------------------------------------------------------
    panel_tsv = "detector\tfluorochrome\tmarker\treagent_lot\n" + "".join(
        f"{d}\t{f}\t{m}\t{lots[f]}\n" for d, f, m, _ in PANEL)
    controls_tsv = "fcs_file\tcontrol_type\tfluorochrome\tdetector\tcarrier\treagent_lot\n" + "".join(
        "\t".join(c) + "\n" for c in controls)
    samples_tsv = "sample_id\tfcs_file\tdonor\tacquisition_order\n" + "".join(
        f"{s}\t{s}.fcs\tD{donors[k]}\t{k + 1}\n" for k, s in enumerate(sids))
    last1, first2 = sids[n_block1 - 1], sids[n_block1]
    notes_md = "\n".join([
        f"# Acquisition notes - experiment IMM-FC-012_PBMC_lineage ({date})", "",
        "- Instrument: LSRFortessa X-20 (405/488/561/640 nm), FACSDiva 9.0, daily CS&T passed.",
        "- Parameters recorded: FSC-A, FSC-H, SSC-A, eight fluorescence areas (see panel.tsv), Time. "
        "Detector voltages are stored in the files.",
        "- Staining: surface antibody cocktail plus Zombie Aqua fixable viability dye, 1% PFA fixation. "
        "The reagent lots of the cocktail are listed in panel.tsv, the lot used for each control tube in "
        "controls.tsv.",
        "- Material: cryopreserved PBMC, thawed on the day of staining; trypan blue viability after thaw "
        "between 80% and 92% depending on the donor.",
        "- Tubes were vortexed briefly and run at medium flow rate (300-350 events/s) without straining "
        "through a filter cap.",
        f"- Experiment template: created on {tmpl_date}, when the panel was first set up on this instrument; "
        f"its compensation values date from that set-up. The detector voltages of the panel were re-optimised "
        f"after the laser service of {svc_date}. The compensation setup of the acquisition software was not run "
        "for this session; the experiment was recorded with the template compensation and exported as FCS 3.1 "
        "(the export stores the detector signals as recorded and the template values in the file header).",
        "- Controls recorded at the start of the session with the panel voltages: one unstained PBMC tube; "
        "antibody-capture beads (positive and negative beads in the same tube) for the seven antibody "
        "conjugates, stained the day before from the vials that were open at the time; Zombie Aqua on a 1:1 mix "
        "of live and heat-killed PBMC; an aliquot of the reference pool stained with CD8 PerCP-Cy5.5 alone, taken "
        f"from the vial that was used for the cocktail ({percp_cells}); the HLA-DR FMO tube (reference pool, "
        "complete cocktail without HLA-DR PE-Cy7).",
        f"- After {last1} the HLA-DR signal of the T cells was judged too close to the axis, and the voltage of "
        f"YG780-A was raised from {dr_v1} V to {dr_v2} V. {first2} and all later tubes were recorded with the new "
        f"voltage. Recorded again with the new voltage before continuing: the PE-Cy7 capture-bead tube (same tube, "
        f"{rerun_pecy7}) and the HLA-DR FMO tube (same tube, {rerun_fmo}). No other control was recorded again.",
        f"- {sids[shift_idx]}: large events were piling up in the top FSC channel during the first attempt, "
        f"so the tube was re-acquired after the FSC voltage had been lowered from {FSC_V} V to {shift_v} V. "
        "Only the re-acquisition was exported. The fluorescence detectors were not touched for this.",
    ]) + "\n"

    # ---- battery: truth-derived probes / attacks -----------------------------------------------
    prng = np.random.default_rng(seed + 4711)

    def tol(p, pct):
        return max(ACT_ABS, ACT_REL * pct) if p in ACT else max(TOL_ABS, TOL_REL * pct)

    def perturbed(scale_fn, count_scale):
        tab = {}
        for s in sids:
            tab[s] = {}
            for p in POP_NAMES4:
                c, pct = truth_tab[s][p]
                v = float(np.clip(pct + scale_fn(s, p) * tol(p, pct), 0.0, 100.0))
                tab[s][p] = (int(round(c * (1 + count_scale * prng.choice([-1, 1])))), v)
        return tab

    def inside(s, p):
        pct = truth_tab[s][p][1]
        sign = -1.0 if pct > 99.0 else float(prng.choice([-1, 1]))
        return 0.6 * sign

    def outside(s, p):
        return -1.4 if truth_tab[s][p][1] > 50.0 else 1.4

    shifted_pop = {s: LINEAGE[int(prng.integers(3, len(LINEAGE)))] for s in sids}

    def one_outside(s, p):
        if p != shifted_pop[s]:
            return 0.0
        return -1.5 if truth_tab[s][p][1] > 50.0 else 1.5

    def act_outside(s, p):
        return 1.5 if p in ACT else 0.0

    def one_act_outside(s, p):
        return -1.5 if (p == "activated_CD8_T" and s == sids[-1]) else 0.0

    def jitter(M, amp):
        J = M + prng.uniform(-amp, amp, size=M.shape)
        np.fill_diagonal(J, 1.0)
        return J

    mats_in = {s: jitter(mats[s], 0.008) for s in sids}
    mats_out = {s: mats[s].copy() for s in sids}
    mats_out[sids[-1]][iCD56, iCD19] += 0.06
    mats_block1 = {s: S for s in sids}
    textbook = {"cells": 90.0, "singlets": 95.0, "live": 90.0, "lymphocytes": 80.0, "monocytes": 20.0,
                "T_cells": 70.0, "CD4_T": 60.0, "CD8_T": 30.0, "DP_T": 1.0, "DN_T": 5.0, "B_cells": 10.0,
                "NK_cells": 12.0, "activated_CD4_T": 5.0, "activated_CD8_T": 15.0}
    prior_tab = {}
    for s in sids:
        c = {"all_events": 40000.0}
        prior_tab[s] = {}
        for p in POP_NAMES4:
            c[p] = c[PARENT4[p]] * textbook[p] / 100.0
            prior_tab[s][p] = (int(round(c[p])), textbook[p])
    eye = {s: np.eye(NF) for s in sids}
    truth_freq = _freq_tsv(sids, truth_tab)

    pipeline_py = {"flow_pipeline_l4.py": ASSETS / "flow_pipeline_l4.py"}
    run = "python3 /alegen_script/flow_pipeline_l4.py "
    ceil = 0.75      # pass_threshold - 0.15
    probes = [
        Probe("planted_truth", "submission written directly from the planted matrices, the event labels and "
              "the SOP rule applied to them", 0.999,
              files={"spillover.csv": _spill_long(sids, mats), "frequencies.tsv": truth_freq}),
        Probe("within_tolerance", "matrix entries +-0.008, every percentage moved by 0.6 x its tolerance, "
              "counts +-1.5%", 0.99,
              files={"spillover.csv": _spill_long(sids, mats_in),
                     "frequencies.tsv": _freq_tsv(sids, perturbed(inside, 0.015))}),
        Probe("independent_solution", "second implementation: fcsparser reader, regression spillover, gain "
              "factor from mixture-model bead peaks, mixture-model gates", 0.99,
              script="SOLUTION=2 bash /solution/solve.sh"),
    ]
    attacks = [
        Attack("l3_reference", "the level-3 reference method on the level-4 inputs: one matrix for all samples "
               "from the first-listed single-stain control of each fluorochrome, level-3 gates, HLA-DR + by the "
               "density-valley rule used for every other marker (no FMO gate, no acceptance check)", ceil,
               script=run + "--l3", script_files=pipeline_py),
        Attack("l3_reference_last_control", "as l3_reference, but with the last-listed single-stain control of "
               "each fluorochrome (what the unmodified level-3 script does when a fluorochrome is listed twice)",
               ceil, script=run + "--l3 --tandem last", script_files=pipeline_py),
        Attack("careless_embedded_matrix", "P_embedded_matrix ignored: samples and FMO tubes compensated with "
               "the $SPILLOVER keyword of the files, which is also reported; FMO gate applied, acceptance check "
               "not acted on", ceil,
               script=run + "--comp embedded", script_files=pipeline_py),
        Attack("careless_one_matrix", "P_voltage_change ignored: the session-start matrix (correct controls) and "
               "the session-start FMO tube are used for every sample; acceptance check not acted on", ceil,
               script=run + "--blocks single", script_files=pipeline_py),
        Attack("careless_row_only", "P_voltage_change half handled: for the later samples only the PE-Cy7 row "
               "is replaced from the re-recorded bead tube (and the re-recorded FMO is used); the spill of the "
               "other fluorochromes into the HLA-DR detector is left at its session-start value; acceptance check "
               "not acted on", ceil,
               script=run + "--blocks row_only", script_files=pipeline_py),
        Attack("careless_bead_tandem", "P_tandem_lot ignored: PerCP-Cy5.5 row from the bead control stained "
               "with the other lot (detector settings handled correctly); acceptance check not acted on", ceil,
               script=run + "--tandem first", script_files=pipeline_py),
        Attack("careless_no_fmo", "Q_fmo_gate ignored: correct matrices, HLA-DR threshold = 99th percentile of "
               "the unstained control instead of the FMO control", ceil,
               script=run + "--dr-gate unstained", script_files=pipeline_py),
        Attack("careless_fmo_uncompensated", "Q_fmo_gate mis-applied: correct matrices, FMO percentile and "
               "sample gate taken on the uncompensated HLA-DR-detector signal (lineage gates on compensated "
               "data)", ceil,
               script=run + "--dr-gate fmo_raw", script_files=pipeline_py),
        Attack("textbook_prior", "identity matrix and textbook PBMC percentages, no data analysed", 0.1,
               files={"spillover.csv": _spill_long(sids, eye), "frequencies.tsv": _freq_tsv(sids, prior_tab)}),
        Attack("block1_matrix_true_frequencies", "true frequencies, planted block-1 matrix reported for every "
               "sample", 0.85,
               files={"spillover.csv": _spill_long(sids, mats_block1), "frequencies.tsv": truth_freq}),
        Attack("all_percent_beyond_tolerance", "true matrices and counts, every percentage 1.4 x tolerance off", 0.45,
               files={"spillover.csv": _spill_long(sids, mats),
                      "frequencies.tsv": _freq_tsv(sids, perturbed(outside, 0.0))}),
        Attack("one_population_beyond_tolerance", "true matrices and counts, one lineage population per sample "
               "1.5 x tolerance off", 0.89,
               files={"spillover.csv": _spill_long(sids, mats),
                      "frequencies.tsv": _freq_tsv(sids, perturbed(one_outside, 0.0))}),
        Attack("activation_beyond_tolerance", "true matrices and lineage, both FMO-gated populations of every "
               "sample 1.5 x tolerance off", 0.70,
               files={"spillover.csv": _spill_long(sids, mats),
                      "frequencies.tsv": _freq_tsv(sids, perturbed(act_outside, 0.0))}),
        Attack("one_activation_item_beyond_tolerance", "true matrices and lineage, activated_CD8_T of one sample "
               "1.5 x tolerance off", 0.80,
               files={"spillover.csv": _spill_long(sids, mats),
                      "frequencies.tsv": _freq_tsv(sids, perturbed(one_act_outside, 0.0))}),
        Attack("matrix_entry_beyond_tolerance", "true frequencies, one spillover coefficient of one sample off "
               "by 0.06", 0.85,
               files={"spillover.csv": _spill_long(sids, mats_out), "frequencies.tsv": truth_freq}),
        # level-3 decision points, still present
        Attack("naive_matrix_not_applied", "matrices computed and reported, samples and FMO gated on the "
               "uncompensated signals", 0.45, script=run + "--comp none", script_files=pipeline_py),
        Attack("naive_no_doublet_gate", "correct compensation, no aggregate exclusion (singlets = cells)", 0.45,
               script=run + "--no-doublet-gate", script_files=pipeline_py),
        Attack("naive_no_viability_gate", "correct compensation and singlets, dead cells not excluded", 0.45,
               script=run + "--no-viability-gate", script_files=pipeline_py),
        Attack("naive_global_thresholds", "correct compensation and gates, but one gate template set on the "
               "first sample and applied unchanged to all samples", 0.60,
               script=run + "--global-thresholds", script_files=pipeline_py),
        Attack("naive_default_pipeline", "embedded matrix, no doublet gate, no viability gate, one gate template, "
               "no FMO", 0.20,
               script=run + "--comp embedded --no-doublet-gate --no-viability-gate --global-thresholds "
                            "--dr-gate valley", script_files=pipeline_py),
    ]
    b1 = [s for s in sids if block[s] == 1]
    b2 = [s for s in sids if block[s] == 2]
    traps = [
        Trap("compensation_required",
             "Nothing was compensated at acquisition; uncompensated, CD3-BV421 spills into the viability "
             "detector, CD4-PE into the CD8 detector and CD56-APC into the CD19 detector at the level of true "
             "positives.",
             "acquisition_notes.md: compensation setup not run; diagonal populations in any V450-A vs V525-A or "
             "YG586-A vs B710-A plot", "Understanding / domain knowledge gap"),
        Trap("doublets",
             "6-12% of cell events are aggregates of two cells; their summed fluorescence inflates DP_T and the "
             "lymphocyte count.",
             "gating_sop.md lists a singlets population; FSC-A vs FSC-H shows a second population off the "
             "diagonal; notes say tubes were not strained", "Approach / skipped QC step"),
        Trap("dead_cells",
             "8-18% of the single cells are dead and take up every antibody non-specifically.",
             "gating_sop.md defines live as viability dye negative; notes report 80-92% post-thaw viability",
             "Approach / skipped QC step"),
        Trap("per_sample_shift",
             f"{sids[shift_idx]} was recorded with the FSC voltage lowered (gain x{fsc_gain:.2f}).",
             "acquisition_notes.md states the re-acquisition; $P1V/$P2V differ in that file",
             "Verification / per-sample QC"),
        Trap("P_embedded_matrix",
             "Every file carries a well-formed, non-identity $SPILLOVER (APPLY COMPENSATION = TRUE) from the "
             f"experiment template of {tmpl_date}; the voltages were re-optimised later, so it is wrong for this "
             f"session (PerCP-Cy5.5 -> YG780-A {S_old[iCD8, iDR]:.3f} instead of {S[iCD8, iDR]:.3f}).",
             "acquisition_notes.md gives the dates of the template and of the voltage re-optimisation and says "
             "that the compensation setup was not run; the single-stain controls of the session give other "
             "coefficients; with the embedded matrix the lineage gates break down (a population of the SOP's "
             "acceptance check is missing in the FMO tube)",
             "Prior decoupled from truth / trusting embedded metadata"),
        Trap("P_voltage_change",
             f"The YG780-A voltage was raised from {dr_v1} V to {dr_v2} V (gain x{g:.2f}) after {last1}: "
             f"{', '.join(b1)} need the session-start matrix, {', '.join(b2)} a second one. Only the PE-Cy7 bead "
             "tube and the FMO tube were recorded again; the spill of every other fluorochrome into YG780-A "
             "scaled by the same gain factor and has to be derived (bead brightness ratio of the two PE-Cy7 "
             "tubes, ratio of the two PE-Cy7 rows, or the FMO tube itself).",
             "acquisition_notes.md states the voltage change and which tubes were recorded again; $P9V differs "
             "between the files; the SOP's acceptance check on the re-recorded FMO tube fails when only the "
             "PE-Cy7 row is replaced (CD8_T median far above DN_T)",
             "Prior decoupled from truth / one matrix per experiment"),
        Trap("P_tandem_lot",
             f"The PerCP-Cy5.5 bead control was stained from lot {bead_lot_percp}, the samples from lot "
             f"{lots[FLUORS[iCD8]]}; the bead lot spills {tandem_delta:.3f} less into YG780-A. The cell-based "
             "single stain from the cocktail vial gives the row that applies to the samples.",
             "reagent_lot columns of panel.tsv and controls.tsv differ for that bead control and agree for the "
             "cell-based one; acquisition_notes.md says when and from which vials the beads and the extra tube "
             "were stained; the two PerCP-Cy5.5 controls give different rows; the SOP's acceptance check on the "
             "FMO tube fails with the bead-derived row",
             "Prior decoupled from truth / beads as default compensation control"),
        Trap("Q_fmo_gate",
             "HLA-DR on T cells has no valley; the SOP defines the gate on the FMO control (99th percentile, "
             "same parent, compensated) and states an acceptance criterion for the compensation on the same "
             f"tube (medians of five lymphocyte populations within {QC_UNITS} units). A gate from the unstained "
             "control, from uncompensated data or from a valley search gives a plausible but out-of-band "
             "percentage; a residual spillover of several hundred units widens the FMO background and loses dim "
             "positives. The default matrices (embedded, bead lot, session-start column) miss the criterion by a "
             "factor 16 or more; the analyst has to find the cause or tune the matrix on the FMO tube.",
             "gating_sop.md sections 1-3; the FMO tubes are listed in controls.tsv; the acceptance criterion "
             "lets the analyst measure the residual spillover on cells and says that a failing compensation has "
             "to be corrected",
             "Quality bar with an in-data control"),
    ]

    files = {f"fcs/{p.name}": p for p in sorted(fcs_dir.glob("*.fcs"))}
    extra = "RUN pip install --no-deps fcsparser==0.2.8\n"
    return TaskBuild(
        task_id=f"flow-compensation-s{seed:04d}-l{level}", family="flow_compensation", seed=seed, level=level,
        title="8-colour PBMC flow cytometry: per-sample spillover matrices, lineage and FMO-gated activation "
              "frequencies",
        summary="Derive the spillover matrix that applies to every sample of an 8-colour PBMC session whose "
                "detector settings and control reagents changed, and report the lineage and FMO-gated HLA-DR "
                "population frequencies of every sample from raw FCS 3.1 files.",
        subdomain="3.3 Cell & Imaging Biology",
        software=["FlowKit", "FlowIO", "fcsparser", "scikit-learn", "Python", "numpy/scipy/pandas"],
        instruction_md=INSTRUCTION.format(n_samples=n_samples, n_ctrl=len(controls)),
        agent_dockerfile=dockerlib.agent_dockerfile(
            conda=["numpy=2.3.5", "scipy=1.16.3", "pandas=2.3.3", "scikit-learn=1.7.2", "flowio=1.4.0",
                   "six=1.17.0", "pip"],
            pip=["flowkit==1.3.2"], extra=extra),
        verifier_dockerfile=dockerlib.verifier_dockerfile(pip=[]),
        input_files={**files, "samples.tsv": samples_tsv, "controls.tsv": controls_tsv, "panel.tsv": panel_tsv,
                     "gating_sop.md": SOP.format(pct=FMO_PCT, qc=QC_UNITS), "acquisition_notes.md": notes_md},
        solution_files={"solve.sh": ASSETS / "solve_l4.sh", "flow_pipeline_l4.py": ASSETS / "flow_pipeline_l4.py",
                        "solve_alt_l4.py": ASSETS / "solve_alt_l4.py"},
        tests_files={"score.py": ASSETS / "score_l4.py", "reference/truth.json": json.dumps(truth)},
        traps=traps, attacks=attacks, probes=probes,
        capabilities=["instrument-file metadata inspection", "compensation from single-stain controls",
                      "reasoning about detector gain and spillover", "control selection by reagent lot",
                      "FMO-based gating", "self-check against an in-data control",
                      "hierarchical gating with per-sample QC", "exact output contracts"],
        judge_explanation_md=JUDGE, pass_threshold=0.9,
        oracle_steps=17, oracle_runtime_min=0.05, cpus=4, memory_mb=8000, storage_mb=10000,
    )
