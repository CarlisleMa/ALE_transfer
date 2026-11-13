"""Family: 8-colour PBMC flow cytometry, raw FCS -> spillover matrix -> gated population frequencies.

Truth-by-construction: a forward simulator draws, for every event, a hidden class (debris, cell
aggregate, dead cell, CD4/CD8/DP/DN T cell, B cell, NK cell, other lymphocyte, monocyte), its
fluorochrome amounts and its scatter, mixes the fluorochromes through a planted spillover matrix,
adds photon and electronic noise and writes FCS 3.1 list-mode files. The reference is the planted
matrix and the per-sample counts of the hidden classes; no analysis output is used as truth.

Decision points (level 1: none; level 2: the first two; level 3: all four):
  compensation_required  exported files carry the acquisition default (identity) $SPILLOVER although
                         nothing was compensated; the matrix has to come from the single-stain controls
  doublets               6-12% of the cell events are aggregates of two cells (off the FSC-A/FSC-H diagonal)
  dead_cells             8-18% of the single cells are dead and bind every antibody non-specifically
  per_sample_shift       one sample was recorded with a lowered FSC voltage (gain x0.20-0.27)
Level 4 (all four plus P_embedded_matrix, P_voltage_change, P_tandem_lot, Q_fmo_gate) is in l4.py.
Level 7 (level 4 plus T_acquisition_disturbance, a stretch of two runs that only the Time axis shows) is in l7.py.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from ... import dockerlib
from ...spec import Attack, Probe, TaskBuild, Trap

ASSETS = Path(__file__).parent / "assets"

# detector, fluorochrome, marker, PMT voltage
PANEL = [
    ("V450-A", "BV421", "CD3", 478),
    ("V525-A", "Zombie Aqua", "Viability", 522),
    ("B530-A", "FITC", "CD14", 506),
    ("B710-A", "PerCP-Cy5.5", "CD8", 612),
    ("YG586-A", "PE", "CD4", 498),
    ("YG780-A", "PE-Cy7", "HLA-DR", 641),
    ("R670-A", "APC", "CD56", 587),
    ("R780-A", "APC-Cy7", "CD19", 655),
]
DETS = [p[0] for p in PANEL]
FLUORS = [p[1] for p in PANEL]
MARKERS = [p[2] for p in PANEL]
iCD3, iVIA, iCD14, iCD8, iCD4, iDR, iCD56, iCD19 = range(8)
NF = 8
PNN = ["FSC-A", "FSC-H", "SSC-A"] + DETS + ["Time"]
FSC_V, SSC_V = 420, 265
RANGE = 262144

# planted spillover: (source fluorochrome, receiving detector) -> uniform range; everything else U(0, 0.003)
SPILL_RANGES = {
    (iCD3, iVIA): (0.18, 0.32), (iCD3, iCD14): (0.004, 0.015),
    (iVIA, iCD3): (0.03, 0.08), (iVIA, iCD14): (0.04, 0.10), (iVIA, iCD4): (0.004, 0.015),
    (iCD14, iVIA): (0.02, 0.06), (iCD14, iCD4): (0.06, 0.14), (iCD14, iCD8): (0.02, 0.05),
    (iCD8, iDR): (0.15, 0.30), (iCD8, iCD56): (0.03, 0.08), (iCD8, iCD19): (0.02, 0.06),
    (iCD4, iCD8): (0.18, 0.32), (iCD4, iCD14): (0.004, 0.015), (iCD4, iDR): (0.02, 0.06),
    (iDR, iCD19): (0.06, 0.15), (iDR, iCD4): (0.01, 0.04), (iDR, iCD8): (0.02, 0.06),
    (iCD56, iCD19): (0.20, 0.38), (iCD56, iCD8): (0.02, 0.06), (iCD56, iDR): (0.004, 0.015),
    (iCD19, iCD56): (0.04, 0.10), (iCD19, iDR): (0.04, 0.10),
}
# median fluorochrome amount of a positive cell, in units of the primary detector
POS = np.array([40000.0, 12000.0, 15000.0, 8000.0, 30000.0, 8000.0, 20000.0, 5000.0])
BEAD = np.array([60000.0, 0.0, 35000.0, 25000.0, 70000.0, 30000.0, 50000.0, 20000.0])
AF_CELL = np.array([70.0, 110.0, 120.0, 60.0, 80.0, 40.0, 50.0, 40.0])
AF_BEAD = np.array([45.0, 60.0, 70.0, 40.0, 50.0, 35.0, 40.0, 35.0])

DEBRIS, DOUBLET, DEAD, CD4T, CD8T, DPT, DNT, BC, NK, OTHER, MONO = range(11)
POPS = [  # name, parent, classes
    ("cells", "all_events", list(range(1, 11))),
    ("singlets", "cells", list(range(2, 11))),
    ("live", "singlets", list(range(3, 11))),
    ("lymphocytes", "live", list(range(3, 10))),
    ("monocytes", "live", [MONO]),
    ("T_cells", "lymphocytes", [CD4T, CD8T, DPT, DNT]),
    ("CD4_T", "T_cells", [CD4T]),
    ("CD8_T", "T_cells", [CD8T]),
    ("DP_T", "T_cells", [DPT]),
    ("DN_T", "T_cells", [DNT]),
    ("B_cells", "lymphocytes", [BC]),
    ("NK_cells", "lymphocytes", [NK]),
]
POP_NAMES = [p[0] for p in POPS]
PARENT = {p[0]: p[1] for p in POPS}
# scoring tolerances (kept in the generator so probes/attacks and the verifier use the same numbers)
TOL_ABS, TOL_REL, TOL_COUNT = 0.5, 0.03, 0.03


# --------------------------------------------------------------------------------------- FCS 3.1
def write_fcs(path: Path, data: np.ndarray, pns: list[str], pnv: list[str], keywords: dict[str, str]) -> None:
    """Minimal FCS 3.1 list-mode writer (float32, little endian, form-feed delimiter as FACSDiva)."""
    n, p = data.shape
    assert p == len(PNN)
    d = "\x0c"
    kv: list[tuple[str, str]] = [
        ("$BEGINANALYSIS", "0"), ("$ENDANALYSIS", "0"), ("$BEGINSTEXT", "0"), ("$ENDSTEXT", "0"),
        ("$BYTEORD", "1,2,3,4"), ("$DATATYPE", "F"), ("$MODE", "L"), ("$NEXTDATA", "0"),
        ("$PAR", str(p)), ("$TOT", str(n)),
    ]
    for k, v in keywords.items():
        assert v != "" and d not in v
        kv.append((k, v))
    for i in range(p):
        j = i + 1
        kv += [(f"$P{j}N", PNN[i]), (f"$P{j}B", "32"), (f"$P{j}E", "0,0"), (f"$P{j}R", str(RANGE))]
        if pns[i]:
            kv.append((f"$P{j}S", pns[i]))
        if pnv[i]:
            kv.append((f"$P{j}V", pnv[i]))
    raw = np.ascontiguousarray(data, dtype="<f4").tobytes()
    text_start = 256

    def text(b: int, e: int) -> bytes:
        items = [("$BEGINDATA", str(b)), ("$ENDDATA", str(e))] + kv
        return (d + d.join(k + d + v for k, v in items) + d).encode("utf-8")

    b = text_start + len(text(0, 0))
    for _ in range(5):
        t = text(b, b + len(raw) - 1)
        nb = text_start + len(t)
        if nb == b:
            break
        b = nb
    t = text(b, b + len(raw) - 1)
    assert text_start + len(t) == b
    header = "FCS3.1    " + "".join(f"{v:>8d}" for v in (text_start, b - 1, b, b + len(raw) - 1, 0, 0))
    assert len(header) == 58
    with open(path, "wb") as fh:
        fh.write(header.encode("ascii") + b" " * (text_start - 58) + t + raw)


# --------------------------------------------------------------------------------------- simulator
def _planted_spillover(rng: np.random.Generator) -> np.ndarray:
    S = rng.uniform(0.0, 0.003, size=(NF, NF))
    for (i, j), (lo, hi) in SPILL_RANGES.items():
        S[i, j] = rng.uniform(lo, hi)
    np.fill_diagonal(S, 1.0)
    return np.round(S, 4)


def _singles(rng: np.random.Generator, cls: np.ndarray, stain: np.ndarray):
    """True fluorochrome amounts, scatter and autofluorescence of single events of the given classes."""
    n = len(cls)
    T = np.zeros((n, NF))

    def pos(mask, i, sig, scale=1.0):
        T[mask, i] = POS[i] * stain[i] * scale * np.exp(sig * rng.standard_normal(int(mask.sum())))

    T[:, iVIA] = 60.0 * np.exp(0.4 * rng.standard_normal(n))          # dim dye uptake by live cells
    is_t = (cls >= CD4T) & (cls <= DNT)
    pos(is_t, iCD3, 0.25)
    pos((cls == CD4T) | (cls == DPT), iCD4, 0.25)
    pos((cls == CD8T) | (cls == DPT), iCD8, 0.28)
    pos(cls == BC, iCD19, 0.28)
    pos(cls == BC, iDR, 0.35)
    pos(cls == NK, iCD56, 0.35)
    pos(cls == MONO, iCD14, 0.30)
    pos(cls == MONO, iDR, 0.35, 1.5)
    u = rng.random(n)
    pos((cls == OTHER) & (u < 0.5), iDR, 0.5, 0.75)
    pos(((cls == CD8T) & (u < 0.12)) | ((cls == DNT) & (u < 0.20)), iCD56, 0.30, 0.3)   # CD56+ T cells
    pos(is_t & (u > 0.92), iDR, 0.40, 0.375)                                           # activated T cells
    dead = cls == DEAD
    nd = int(dead.sum())
    if nd:
        pos(dead, iVIA, 0.45)
        sticky = np.exp(0.8 * rng.standard_normal(nd))                 # per-cell non-specific binding
        for i in range(NF):
            if i != iVIA:
                T[dead, i] = 0.06 * POS[i] * sticky * np.exp(0.5 * rng.standard_normal(nd))
    deb = cls == DEBRIS
    T[deb, iVIA] = 400.0 * np.exp(0.8 * rng.standard_normal(int(deb.sum())))

    fsc = np.empty(n)
    ssc = np.empty(n)
    lym = (cls >= CD4T) & (cls <= OTHER)
    mono = cls == MONO
    for mask, fm, fs, sm, ss in ((lym, 90000, 6000, 25000, 0.20), (mono, 128000, 9000, 90000, 0.18),
                                 (dead, 66000, 6000, 30000, 0.25)):
        k = int(mask.sum())
        fsc[mask] = rng.normal(fm, fs, k)
        ssc[mask] = sm * np.exp(ss * rng.standard_normal(k))
    k = int(deb.sum())
    fsc[deb] = 10000.0 * np.exp(0.4 * rng.standard_normal(k))
    ssc[deb] = 8000.0 * np.exp(0.6 * rng.standard_normal(k))

    cellf = np.exp(0.25 * rng.standard_normal(n))
    cellf[mono] *= 2.5
    cellf[deb] *= 0.3
    af = cellf[:, None] * AF_CELL[None, :] * np.exp(0.1 * rng.standard_normal((n, NF)))
    return T, fsc, ssc, af


def _detect(rng, T, af, S, q, en):
    mean = T @ S + af
    return mean + np.sqrt(np.maximum(mean, 0.0) * q) * rng.standard_normal(mean.shape) \
        + en * rng.standard_normal(mean.shape)


def _assemble(rng, fsc, fsch, ssc, fl, labels=None):
    n = len(fsc)
    order = rng.permutation(n)
    t = np.sort(rng.uniform(0, n / 3.2, size=n))                 # ~320 events/s, $TIMESTEP 0.01 s
    data = np.column_stack([fsc, fsch, ssc, fl])[order]
    data[:, :3] = np.clip(data[:, :3], 1.0, RANGE - 1)
    data[:, 3:] = np.clip(data[:, 3:], -RANGE / 100, RANGE - 1)
    data = np.column_stack([data, t]).astype(np.float32)
    return data, (labels[order] if labels is not None else None)


def _sample(rng, level, S, q, en):
    n = int(rng.integers(36000, 44001))
    n_deb = int(round(n * rng.uniform(0.06, 0.15)))
    n_cells = n - n_deb
    n_dbl = int(round(n_cells * rng.uniform(0.06, 0.12))) if level >= 2 else 0
    n_single = n_cells - n_dbl
    f_dead = rng.uniform(0.08, 0.18) if level >= 3 else rng.uniform(0.004, 0.012)
    n_dead = int(round(n_single * f_dead))
    n_live = n_single - n_dead
    n_mono = int(round(n_live * rng.uniform(0.10, 0.25)))
    n_lym = n_live - n_mono
    p_l = np.array([rng.uniform(0.55, 0.75), rng.uniform(0.06, 0.15), rng.uniform(0.06, 0.18),
                    rng.uniform(0.015, 0.04)])
    n_t, n_b, n_nk, n_other = rng.multinomial(n_lym, p_l / p_l.sum())
    p_t = np.array([rng.uniform(0.50, 0.68), rng.uniform(0.22, 0.38), rng.uniform(0.008, 0.025),
                    rng.uniform(0.03, 0.07)])
    n4, n8, ndp, ndn = rng.multinomial(n_t, p_t / p_t.sum())
    counts = {DEBRIS: n_deb, DEAD: n_dead, CD4T: n4, CD8T: n8, DPT: ndp, DNT: ndn, BC: n_b, NK: n_nk,
              OTHER: n_other, MONO: n_mono}
    stain = np.exp(0.08 * rng.standard_normal(NF))               # donor / staining variation
    cls = np.concatenate([np.full(c, k) for k, c in counts.items()])
    T, fsc, ssc, af = _singles(rng, cls, stain)
    fsch = fsc * 0.92 * (1 + 0.02 * rng.standard_normal(len(cls)))
    if n_dbl:
        cell_classes = [k for k in counts if k != DEBRIS]
        p = np.array([counts[k] for k in cell_classes], dtype=float)
        a = rng.choice(cell_classes, size=n_dbl, p=p / p.sum())
        b = rng.choice(cell_classes, size=n_dbl, p=p / p.sum())
        Ta, fa, sa, afa = _singles(rng, a, stain)
        Tb, fb, sb, afb = _singles(rng, b, stain)
        big, small = np.maximum(fa, fb), np.minimum(fa, fb)
        fd = big + small * rng.uniform(0.8, 1.0, n_dbl)
        fhd = big * 0.92 * rng.uniform(1.0, 1.05, n_dbl) * (1 + 0.02 * rng.standard_normal(n_dbl))
        T = np.vstack([T, Ta + Tb])
        af = np.vstack([af, afa + afb])
        fsc = np.concatenate([fsc, fd])
        fsch = np.concatenate([fsch, fhd])
        ssc = np.concatenate([ssc, sa + sb])
        cls = np.concatenate([cls, np.full(n_dbl, DOUBLET)])
    fl = _detect(rng, T, af, S, q, en)
    return fsc, fsch, ssc, fl, cls


def _bead_control(rng, f, S, q, en, n=6000):
    T = np.zeros((n, NF))
    posb = rng.random(n) < 0.5
    T[posb, f] = BEAD[f] * np.exp(0.12 * rng.standard_normal(int(posb.sum())))
    af = np.exp(0.15 * rng.standard_normal(n))[:, None] * AF_BEAD[None, :] * np.exp(0.1 * rng.standard_normal((n, NF)))
    fsc = rng.normal(70000, 2500, n)
    ssc = 60000 * np.exp(0.08 * rng.standard_normal(n))
    fsch = fsc * 0.92 * (1 + 0.02 * rng.standard_normal(n))
    dbl = rng.random(n) < 0.03                                    # bead aggregates
    T[dbl] *= 2.0
    af[dbl] *= 2.0
    fsc[dbl] *= 1.9
    ssc[dbl] *= 2.0
    return _assemble(rng, fsc, fsch, ssc, _detect(rng, T, af, S, q, en))[0]


def _cell_control(rng, S, q, en, dye: bool, n=8000):
    cls = rng.choice([DEBRIS, OTHER, MONO, DEAD], size=n, p=[0.06, 0.45, 0.07, 0.42] if dye else [0.06, 0.80, 0.13, 0.01])
    T, fsc, ssc, af = _singles(rng, cls, np.ones(NF))
    keep = T[:, iVIA].copy()
    T[:] = 0.0                                                    # no antibodies in this tube
    if dye:
        T[:, iVIA] = keep * np.where(cls == DEAD, 1.25, 1.0)
    fsch = fsc * 0.92 * (1 + 0.02 * rng.standard_normal(n))
    return _assemble(rng, fsc, fsch, ssc, _detect(rng, T, af, S, q, en))[0]


# --------------------------------------------------------------------------------------- texts
INSTRUCTION = """\
You are the cytometry data analyst of an immunology core facility. A study team acquired PBMC
samples with an 8-colour immunophenotyping panel and needs the spillover matrix of the panel and
the size of every population in the lab's gating SOP for every sample.

Task directory:
- `/task`

Inputs (read-only) under `/task/input`:
- `fcs/*.fcs` - FCS 3.1 list-mode files as exported by the acquisition software: {n_samples} PBMC
  samples, {n_ss} single-stain controls and one unstained control
- `samples.tsv` - sample sheet (`sample_id`, `fcs_file`, `donor`, `acquisition_order`)
- `controls.tsv` - control tubes (`fcs_file`, `control_type`, `fluorochrome`, `detector`, `carrier`)
- `panel.tsv` - staining panel (`detector`, `fluorochrome`, `marker`)
- `gating_sop.md` - the lab's gating SOP: population names, parent of each population, definition
- `acquisition_notes.md` - the operator's notes from the acquisition session

Available software: Python 3.12 (on `PATH`) with numpy, scipy, pandas, scikit-learn, FlowKit,
FlowIO and fcsparser. There is no internet access.

Your task:
1. Determine the spillover matrix of the panel on this instrument.
2. For every sample in `samples.tsv`, identify every population of `gating_sop.md` and report its size.

Required outputs under `/task/output`:
- `spillover.csv` - comma-separated. Header: `fluorochrome` followed by the eight fluorescence
  detectors named as in `panel.tsv`. One row per fluorochrome of `panel.tsv` (first field = the
  fluorochrome name as in `panel.tsv`). The entry in row *f*, column *d* is the signal that
  fluorochrome *f* produces in detector *d*, as a fraction of the signal it produces in its own
  primary detector (background-free; the entry of each fluorochrome in its own detector is 1).
- `frequencies.tsv` - tab-separated with the columns `sample`, `population`, `parent`, `count`,
  `percent_of_parent`; exactly one row for every `sample_id` of `samples.tsv` and every population
  of `gating_sop.md` (names exactly as written there). Extra columns are allowed.
  - `parent`: the parent population named in the SOP.
  - `count`: number of events of the sample's FCS file that belong to the population (integer).
  - `percent_of_parent`: 100 x `count` / (count of the parent population). For the parent
    `all_events` the denominator is the total number of events in the FCS file.

Constraints:
- Do not modify anything under `/task/input`; write deliverables only under `/task/output`.
- Every sample in the sample sheet belongs to the study and must be reported.
- Do not ask for confirmation. Execute directly.
"""

SOP = """\
# Gating SOP IMM-FC-012 (PBMC lineage panel, 8 colours)

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

HLA-DR is recorded for the activation read-out of a separate SOP and is not part of this list.
`all_events` denotes every event recorded in the file.
"""

JUDGE = """\
**Where the truth comes from.** Every event is simulated from a hidden class label (debris, aggregate
of two cells, dead cell, CD4/CD8/DP/DN T cell, B cell, NK cell, other lymphocyte, monocyte). The
reference `truth.json` holds the planted spillover matrix and, per sample, the number of events of
each SOP population counted from those labels. No gating or compensation result is used as reference.

**Hard gates (score 0).** Missing/unparseable file; `spillover.csv` not labelled with the eight
fluorochromes and the eight detectors of the panel (the transposed layout is accepted), non-numeric
entries; `frequencies.tsv` not holding exactly one row per sample x SOP population, wrong `parent`,
`count` not a non-negative integer, `percent_of_parent` outside [0, 100].

**Components.**

| weight | component | metric | full credit | zero |
|---|---|---|---|---|
| 0.25 | spillover_matrix | maximum absolute error over the 64 matrix entries | <= 0.02 | >= 0.05 |
| 0.40 | frequencies_every_sample | per sample, fraction of the 12 populations whose `percent_of_parent` is within max(0.5 percentage points, 3% of the true value) of the planted value; **minimum** over samples | 1.0 | <= 0.75 |
| 0.20 | frequencies_overall | fraction of all sample x population items within the same tolerance | 1.0 | <= 0.5 |
| 0.15 | lymphocyte_counts | fraction of samples whose `lymphocytes` count is within 3% of the planted count | 1.0 | <= 0.5 |

**Why these tolerances.** After correct compensation the positive and negative populations are
separated by more than 1.5 decades, so two different correct strategies (density-valley thresholds
on flowio data; mixture-model boundaries on fcsparser data) agree with the planted labels to within
a small fraction of the band (measured values in docs/families/flow_compensation.md), while
each ignored decision point moves several populations per sample by many times the band. The
matrix is recovered from the bead/cell controls to better than 0.01 by either the median-difference
or the regression estimator; the identity matrix is off by 0.2-0.4. The minimum over samples is
used for the main frequency component because a sample that is wholly wrong is a wrong deliverable
and must not be averaged away by the other samples.

**What is deliberately not scored.** Counts of populations other than `lymphocytes` (they are
implied by the percentages), any plot, any compensated FCS export.
"""


def _spill_csv(M: np.ndarray, fmt: str = "{:.5f}") -> str:
    return "fluorochrome," + ",".join(DETS) + "\n" + "".join(
        FLUORS[i] + "," + ",".join(fmt.format(M[i, j]) for j in range(NF)) + "\n" for i in range(NF))


def _freq_tsv(sids, table) -> str:
    """table[sample][population] = (count, percent)."""
    out = ["sample\tpopulation\tparent\tcount\tpercent_of_parent"]
    for s in sids:
        for pop in POP_NAMES:
            c, pct = table[s][pop]
            out.append(f"{s}\t{pop}\t{PARENT[pop]}\t{int(c)}\t{pct:.3f}")
    return "\n".join(out) + "\n"


def generate(seed: int, level: int, workdir: Path) -> TaskBuild:
    if level == 7:                    # level 7 lives in l7.py: the level-4 session plus an acquisition disturbance
        from .l7 import generate_l7
        return generate_l7(seed, level, workdir)
    if level >= 4:                    # level 4 lives in l4.py; levels 1-3 below are unchanged
        from .l4 import generate_l4
        return generate_l4(seed, level, workdir)
    rng = np.random.default_rng(52000 + seed)
    S = _planted_spillover(rng)
    q = rng.uniform(6.0, 12.0, size=NF)          # detector units per photoelectron (photon noise)
    en = rng.uniform(20.0, 35.0, size=NF)        # electronic noise sd
    n_samples = 6
    donors = sorted(int(x) for x in rng.permutation(np.arange(101, 990))[:n_samples])
    sids = [f"PBMC_D{d}" for d in donors]
    shift_idx = int(rng.integers(1, n_samples)) if level >= 3 else -1
    fsc_gain = float(rng.uniform(0.20, 0.27))
    shift_v = int(round(FSC_V * fsc_gain ** (1 / 7.0)))
    day = int(rng.integers(3, 27))
    date = f"{day:02d}-{['FEB', 'MAR', 'APR', 'MAY', 'JUN', 'SEP', 'OCT'][int(rng.integers(0, 7))]}-2025"

    fcs_dir = workdir / "fcs"
    if fcs_dir.exists():
        for f in fcs_dir.glob("*.fcs"):
            f.unlink()
    fcs_dir.mkdir(parents=True, exist_ok=True)
    ident = f"{NF}," + ",".join(DETS) + "," + ",".join("1" if i == j else "0" for i in range(NF) for j in range(NF))
    clock = [9 * 60 + int(rng.integers(5, 40))]

    def emit(fname, data, tube, src, stained: list[int], fscv=FSC_V):
        pns = ["", "", ""] + [MARKERS[i] if i in stained else "" for i in range(NF)] + [""]
        pnv = [str(fscv), str(fscv), str(SSC_V)] + [str(p[3]) for p in PANEL] + [""]
        start = clock[0]
        dur = max(1, int(round(float(data[-1, -1]) / 100 / 60)))
        clock[0] += dur + int(rng.integers(1, 4))
        kw = {"$CYT": "LSRFortessa X-20", "$CYTSN": "X65811-R0428", "$DATE": date,
              "$BTIM": f"{start // 60:02d}:{start % 60:02d}:{int(rng.integers(0, 60)):02d}",
              "$ETIM": f"{(start + dur) // 60:02d}:{(start + dur) % 60:02d}:{int(rng.integers(0, 60)):02d}",
              "$FIL": fname, "$INST": "Immunology Core", "$OP": "core-operator-2", "$SRC": src,
              "$SYS": "Windows 10", "$TIMESTEP": "0.01", "TUBE NAME": tube,
              "EXPERIMENT NAME": "IMM-FC-012_PBMC_lineage", "CREATOR": "BD FACSDiva Software Version 9.0",
              "FSC ASF": "0.78", "APPLY COMPENSATION": "FALSE"}
        if level >= 2:
            kw["$SPILLOVER"] = ident
        write_fcs(fcs_dir / fname, data, pns, pnv, kw)

    # ---- controls (recorded first) -------------------------------------------------------------
    controls = []   # fcs_file, control_type, fluorochrome, detector, carrier
    data = _cell_control(rng, S, q, en, dye=False)
    emit("ctrl_unstained.fcs", data, "Unstained Control", "PBMC pool", [])
    controls.append(("ctrl_unstained.fcs", "unstained", "NA", "NA", "cells"))
    for f in range(NF):
        tag = FLUORS[f].replace(" ", "").replace(".", "").replace("-", "")
        fname = f"ctrl_{tag}.fcs"
        if f == iVIA:
            data = _cell_control(rng, S, q, en, dye=True)
            carrier, src = "cells", "PBMC pool, 1:1 live / heat-killed"
        else:
            data = _bead_control(rng, f, S, q, en)
            carrier, src = "beads", "antibody-capture beads"
        emit(fname, data, f"{FLUORS[f]} Stained Control", src, [f])
        controls.append((fname, "single_stain", FLUORS[f], DETS[f], carrier))

    # ---- samples -------------------------------------------------------------------------------
    truth_counts: dict[str, dict[str, int]] = {}
    n_events: dict[str, int] = {}
    for k, sid in enumerate(sids):
        fsc, fsch, ssc, fl, cls = _sample(rng, level, S, q, en)
        g = fsc_gain if k == shift_idx else 1.0
        data, lab = _assemble(rng, fsc * g, fsch * g, ssc, fl, cls)
        emit(f"{sid}.fcs", data, sid, f"donor {donors[k]}", list(range(NF)), fscv=shift_v if k == shift_idx else FSC_V)
        n_events[sid] = int(len(lab))
        truth_counts[sid] = {name: int(np.isin(lab, classes).sum()) for name, _, classes in POPS}

    def pct_table(counts_of):
        tab = {}
        for s in sids:
            c = dict(counts_of[s])
            c["all_events"] = n_events[s]
            tab[s] = {p: (c[p], 100.0 * c[p] / max(c[PARENT[p]], 1)) for p in POP_NAMES}
        return tab

    truth_tab = pct_table(truth_counts)
    truth = {"fluorochromes": FLUORS, "detectors": DETS, "spillover": S.tolist(), "samples": sids,
             "populations": [{"name": n, "parent": p} for n, p, _ in POPS],
             "n_events": n_events,
             "counts": truth_counts,
             "percent_of_parent": {s: {p: truth_tab[s][p][1] for p in POP_NAMES} for s in sids},
             "tolerance": {"percent_abs": TOL_ABS, "percent_rel": TOL_REL, "count_rel": TOL_COUNT}}

    # ---- text inputs ---------------------------------------------------------------------------
    panel_tsv = "detector\tfluorochrome\tmarker\n" + "".join(f"{d}\t{f}\t{m}\n" for d, f, m, _ in PANEL)
    controls_tsv = "fcs_file\tcontrol_type\tfluorochrome\tdetector\tcarrier\n" + "".join(
        "\t".join(c) + "\n" for c in controls)
    samples_tsv = "sample_id\tfcs_file\tdonor\tacquisition_order\n" + "".join(
        f"{s}\t{s}.fcs\tD{donors[k]}\t{k + 1}\n" for k, s in enumerate(sids))
    notes = [
        f"# Acquisition notes - experiment IMM-FC-012_PBMC_lineage ({date})", "",
        "- Instrument: LSRFortessa X-20 (405/488/561/640 nm), FACSDiva 9.0, daily CS&T passed.",
        "- Parameters recorded: FSC-A, FSC-H, SSC-A, eight fluorescence areas (see panel.tsv), Time. "
        "Detector voltages are stored in the files.",
        "- Staining: surface antibody cocktail plus Zombie Aqua fixable viability dye, 1% PFA fixation.",
        ("- Material: cryopreserved PBMC, thawed on the day of staining; trypan blue viability after thaw "
         "between 80% and 92% depending on the donor." if level >= 3 else
         "- Material: freshly isolated PBMC (Ficoll), stained within 3 h; trypan blue viability > 98% for all donors."),
        ("- Tubes were vortexed briefly and run at medium flow rate (300-350 events/s) without straining "
         "through a filter cap." if level >= 2 else
         "- Tubes were passed through a 35 um strainer cap immediately before acquisition and run at low "
         "flow rate (300-350 events/s)."),
        "- Controls recorded at the start of the session with the panel voltages: one unstained PBMC tube; "
        "antibody-capture beads (positive and negative beads in the same tube) for the seven antibody "
        "conjugates; Zombie Aqua on a 1:1 mix of live and heat-killed PBMC.",
        ("- The compensation setup of the acquisition software was not run; the experiment was recorded "
         "with the template compensation and exported as FCS 3.1 without any correction of the data."
         if level >= 2 else
         "- No compensation was calculated or applied in the acquisition software; the files hold the "
         "signals as recorded (FCS 3.1 export)."),
    ]
    if level >= 3:
        notes.append(
            f"- {sids[shift_idx]}: large events were piling up in the top FSC channel during the first attempt, "
            f"so the tube was re-acquired after the FSC voltage had been lowered from {FSC_V} V to {shift_v} V. "
            "Only the re-acquisition was exported. All other settings were left unchanged and the controls "
            "were not recorded again.")
    notes_md = "\n".join(notes) + "\n"

    # ---- battery: truth-derived probes / attacks -----------------------------------------------
    prng = np.random.default_rng(seed + 4711)

    def tol(p):
        return max(TOL_ABS, TOL_REL * p)

    def perturbed(scale_fn, count_scale):
        tab = {}
        for s in sids:
            tab[s] = {}
            for p in POP_NAMES:
                c, pct = truth_tab[s][p]
                v = float(np.clip(pct + scale_fn(s, p) * tol(pct), 0.0, 100.0))
                tab[s][p] = (int(round(c * (1 + count_scale * prng.choice([-1, 1])))), v)
        return tab

    def inside(s, p):
        pct = truth_tab[s][p][1]
        sign = -1.0 if pct > 99.0 else float(prng.choice([-1, 1]))
        return 0.6 * sign

    def outside(s, p):
        pct = truth_tab[s][p][1]
        return -1.4 if pct > 50.0 else 1.4

    shifted_pop = {s: POP_NAMES[int(prng.integers(3, len(POP_NAMES)))] for s in sids}

    def one_outside(s, p):
        if p != shifted_pop[s]:
            return 0.0
        return -1.5 if truth_tab[s][p][1] > 50.0 else 1.5

    S_in = S + prng.uniform(-0.008, 0.008, size=S.shape)
    np.fill_diagonal(S_in, 1.0)
    S_out = S.copy()
    S_out[iCD56, iCD19] += 0.06
    textbook = {"cells": 90.0, "singlets": 95.0, "live": 90.0, "lymphocytes": 80.0, "monocytes": 20.0,
                "T_cells": 70.0, "CD4_T": 60.0, "CD8_T": 30.0, "DP_T": 1.0, "DN_T": 5.0, "B_cells": 10.0,
                "NK_cells": 12.0}
    prior_tab = {}
    for s in sids:
        c = {"all_events": 40000.0}
        prior_tab[s] = {}
        for p in POP_NAMES:
            c[p] = c[PARENT[p]] * textbook[p] / 100.0
            prior_tab[s][p] = (int(round(c[p])), textbook[p])

    pipeline_py = {"flow_pipeline.py": ASSETS / "flow_pipeline.py"}
    run = "python3 /alegen_script/flow_pipeline.py "
    probes = [
        Probe("planted_truth", "submission written directly from the planted matrix and event labels", 0.999,
              files={"spillover.csv": _spill_csv(S), "frequencies.tsv": _freq_tsv(sids, truth_tab)}),
        Probe("within_tolerance", "matrix entries +-0.008, every percentage moved by 0.6 x its tolerance, "
              "counts +-1.5%", 0.99,
              files={"spillover.csv": _spill_csv(S_in), "frequencies.tsv": _freq_tsv(sids, perturbed(inside, 0.015))}),
        Probe("independent_solution", "second implementation: fcsparser reader, regression spillover, "
              "mixture-model gates", 0.99, script="SOLUTION=2 bash /solution/solve.sh"),
    ]
    attacks = [
        Attack("textbook_prior", "identity matrix and textbook PBMC percentages, no data analysed", 0.1,
               files={"spillover.csv": _spill_csv(np.eye(NF)), "frequencies.tsv": _freq_tsv(sids, prior_tab)}),
        Attack("all_percent_beyond_tolerance", "true matrix and counts, every percentage 1.4 x tolerance off", 0.45,
               files={"spillover.csv": _spill_csv(S), "frequencies.tsv": _freq_tsv(sids, perturbed(outside, 0.0))}),
        Attack("one_population_beyond_tolerance", "true matrix and counts, one lineage population per sample "
               "1.5 x tolerance off", 0.85,
               files={"spillover.csv": _spill_csv(S), "frequencies.tsv": _freq_tsv(sids, perturbed(one_outside, 0.0))}),
        Attack("matrix_entry_beyond_tolerance", "true frequencies, one spillover coefficient off by 0.06", 0.80,
               files={"spillover.csv": _spill_csv(S_out), "frequencies.tsv": _freq_tsv(sids, truth_tab)}),
        Attack("naive_matrix_not_applied", "matrix computed from the controls and reported, but the samples "
               "are gated on the uncompensated signals", 0.40, script=run + "--comp none", script_files=pipeline_py),
    ]
    traps = []
    if level >= 2:
        attacks += [
            Attack("naive_embedded_matrix", "samples compensated with the $SPILLOVER keyword of the files "
                   "(acquisition default), which is also reported as the matrix", 0.20,
                   script=run + "--comp embedded", script_files=pipeline_py),
            Attack("naive_no_doublet_gate", "correct compensation, no aggregate exclusion (singlets = cells)", 0.50,
                   script=run + "--no-doublet-gate", script_files=pipeline_py),
        ]
        traps += [
            Trap("compensation_required",
                 "Files carry an identity $SPILLOVER (acquisition default) although nothing was compensated. "
                 "Uncompensated, CD3-BV421 spills into the viability detector, CD4-PE into the CD8 detector and "
                 "CD56-APC into the CD19 detector at the level of true positives.",
                 "acquisition_notes.md: compensation setup not run, template compensation; keywords $SPILLOVER = "
                 "identity and APPLY COMPENSATION = FALSE; diagonal single-stain controls / diagonal populations "
                 "in any V450-A vs V525-A or YG586-A vs B710-A plot",
                 "Understanding / domain knowledge gap"),
            Trap("doublets",
                 "6-12% of cell events are aggregates of two cells; their fluorescence is the sum of both cells, "
                 "which roughly doubles to quadruples the CD4+CD8+ T-cell percentage and inflates the lymphocyte count "
                 "by 7-11%.",
                 "gating_sop.md lists a singlets population; FSC-A vs FSC-H shows a second population off the "
                 "diagonal; notes say tubes were not strained",
                 "Approach / skipped QC step"),
        ]
    if level >= 3:
        attacks += [
            Attack("naive_no_viability_gate", "correct compensation and singlets, dead cells not excluded "
                   "(live = singlets)", 0.45, script=run + "--no-viability-gate", script_files=pipeline_py),
            Attack("naive_global_thresholds", "correct compensation and gates, but one gate template set on the "
                   "first sample and applied unchanged to all samples", 0.65,
                   script=run + "--global-thresholds", script_files=pipeline_py),
            Attack("naive_default_pipeline", "embedded matrix, no doublet gate, no viability gate, one gate "
                   "template", 0.20,
                   script=run + "--comp embedded --no-doublet-gate --no-viability-gate --global-thresholds",
                   script_files=pipeline_py),
        ]
        traps += [
            Trap("dead_cells",
                 "8-18% of the single cells are dead; they take up every antibody non-specifically (a smear from "
                 "negative to positive in every channel) and, if not removed, inflate the lymphocyte count by ~25% "
                 "and shift T/CD4/CD8/DN/B/NK percentages by several points.",
                 "gating_sop.md defines live as viability dye negative; notes report 80-92% post-thaw viability; "
                 "the V525-A histogram of compensated singlets is bimodal",
                 "Approach / skipped QC step"),
            Trap("per_sample_shift",
                 f"{sids[shift_idx]} was recorded with the FSC voltage lowered (gain x{fsc_gain:.2f}); a debris/"
                 "cell gate set on another sample discards its lymphocytes.",
                 "acquisition_notes.md states the re-acquisition; $P1V/$P2V differ in that file; its FSC-A "
                 "histogram is shifted against the other samples",
                 "Verification / per-sample QC"),
        ]
    else:
        if level >= 2:
            attacks.append(Attack("naive_default_pipeline", "embedded matrix, no doublet gate", 0.20,
                                  script=run + "--comp embedded --no-doublet-gate", script_files=pipeline_py))

    files = {f"fcs/{p.name}": p for p in sorted(fcs_dir.glob("*.fcs"))}
    extra = "RUN pip install --no-deps fcsparser==0.2.8\n"
    return TaskBuild(
        task_id=f"flow-compensation-s{seed:04d}-l{level}", family="flow_compensation", seed=seed, level=level,
        title="8-colour PBMC flow cytometry: spillover matrix and gated population frequencies",
        summary="Derive the spillover matrix of an 8-colour PBMC panel from single-stain controls and report "
                "the SOP population frequencies of every sample from raw FCS 3.1 files.",
        subdomain="3.3 Cell & Imaging Biology",
        software=["FlowKit", "FlowIO", "fcsparser", "scikit-learn", "Python", "numpy/scipy/pandas"],
        instruction_md=INSTRUCTION.format(n_samples=n_samples, n_ss=NF),
        agent_dockerfile=dockerlib.agent_dockerfile(
            conda=["numpy=2.3.5", "scipy=1.16.3", "pandas=2.3.3", "scikit-learn=1.7.2", "flowio=1.4.0",
                   "six=1.17.0", "pip"],
            pip=["flowkit==1.3.2"], extra=extra),
        verifier_dockerfile=dockerlib.verifier_dockerfile(pip=[]),
        input_files={**files, "samples.tsv": samples_tsv, "controls.tsv": controls_tsv, "panel.tsv": panel_tsv,
                     "gating_sop.md": SOP, "acquisition_notes.md": notes_md},
        solution_files={"solve.sh": ASSETS / "solve.sh", "flow_pipeline.py": ASSETS / "flow_pipeline.py",
                        "solve_alt.py": ASSETS / "solve_alt.py"},
        tests_files={"score.py": ASSETS / "score.py", "reference/truth.json": json.dumps(truth)},
        traps=traps, attacks=attacks, probes=probes,
        capabilities=["instrument-file metadata inspection", "compensation from single-stain controls",
                      "hierarchical gating with per-sample QC", "exact output contracts"],
        judge_explanation_md=JUDGE, pass_threshold=0.9,
        oracle_steps=8, oracle_runtime_min=0.5, cpus=4, memory_mb=8000, storage_mb=10000,
    )
