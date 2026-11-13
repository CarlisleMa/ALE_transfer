"""Level 7 of `amplicon_16s`: the level-3 batch plus one silent decision point.

  M_control_swap   The reads of one extraction blank and of one study specimen changed places before
                   sequencing (the index pairs of the two libraries were exchanged): the FASTQ pair filed
                   under the blank's `sample_id` holds the specimen's library and the pair filed under the
                   specimen's `sample_id` holds the blank's. The sample sheet (`sample_type`, `specimen`,
                   `qpcr_16s_copies_per_ul`) is what the wet lab wrote down and is the level-3 sheet.

What an analysis that takes the roles from the sheet gets wrong: contaminant identification treats a stool
community as a blank and a blank as a specimen. Nothing fails; the contaminant list is wrong (reagent ASVs
missed, ASVs of the specimen flagged) and so are the profiles the investigator derives from it. No text
mentions it: instruction, SOP and notes are the level-3 texts (`python -m alegen make --level 7` then removes
the pointer sentences as for level 6).

Truth by construction, unchanged in meaning. The level-3 generator is run unchanged for the same seed (same
strains, same read pairs, same sheet); this module then files the two libraries under each other's name and
permutes the planted truth accordingly:
  * sequences, genus, contaminant status, scored subset: unchanged (they do not depend on the file name);
  * `counts`: defined by the reads of each file, so the two columns change places;
  * `samples[].type`: the kind of template the file holds (the scorer takes the minimum over the libraries
    that hold a specimen; the columns of the two blanks are not scored, as at level 3). The sheet's entry
    is kept next to it as `sheet_type`.
The two FASTQ pairs are written again under the other name (same reads, same order), so that the name
stored in the gzip header is the name of the file. Read names carry no sample information.

Scorer: `assets/score_l7.py`, a thin wrapper around the level-3 verifier, which is shipped unchanged next to
it as `tests/score_l3.py`. The count table is accepted under either filing of the two libraries: as the files
are (the stored truth) or by the sample each library really is (the same two columns the other way round).
The level-3 logic is run for both readings and the better total counts; nothing else in the scoring differs
from level 3 (maintainer decision of 2026-11-12, see docs/families/amplicon_16s.md).

The pair is drawn from the seed (own random stream) among the 2 blanks x the specimens with a reagent share
below 5 %. A draw is used only if, on the planted counts, (1) the QC named in the task card separates the two
libraries from the rest of the batch, (2) the reference contaminant test with the sheet's roles gives a
submission that scores at most pass - 0.15, (3) the same test with the real roles returns exactly the planted
contaminants and (4) so does the test with both libraries left out. (1) makes the decision point discoverable,
(2) makes it matter, (3) and (4) make the two correct routes a property of the instance.

Levels 1-6 do not import this module.
"""
from __future__ import annotations

import gzip
import json
import math
from pathlib import Path

import numpy as np

from ...spec import Attack, Probe, Trap
from . import ASSETS, PASS, generate

RICH_SHARE = 0.005                # an ASV "is there" from 0.5 % of the read pairs of a library
MAX_REAGENT_SHARE = 0.05          # the specimen that changes places is a community, not reagent background
# requirements on the planted counts (see module docstring, point 1)
BC_SAME_MAX, BC_NEXT_MIN, BC_DIFF_MIN = 0.15, 0.33, 0.75
RICH_BLANK_MIN, RICH_SPECIMEN_MAX = 0.75, 0.5          # relative to the poorest of the other specimens
COPIES_FACTOR_MIN, COPIES_SCATTER_MAX = 10.0, 2.5

L7_JUDGE = """

**Level-7 decision point (a blank and a specimen that changed places).** The FASTQ pair filed under one
extraction blank holds the library of a study specimen (a stool community with less than 5 % reagent
background) and the pair filed under that specimen holds the library of the blank; the sample sheet is
unchanged. The planted truth is that of the level-3 batch of the same seed with nothing redefined:
sequences, genus and contaminant status do not depend on a file name; the count of an ASV in a sample is
the number of read pairs of the library filed under that `sample_id` (instruction; SOP item 2), so the two
count columns change places; a reagent contaminant is defined by where its template DNA comes from (SOP
item 4), so the contaminant set is the same. Gates, components, weights and tolerances are those of level 3:
`tests/score_l3.py` is the level-3 verifier, byte for byte, and `tests/score.py` only runs it. `sample_counts`
is the minimum over the libraries that hold a specimen, i.e. the seven unaffected specimens and the file
labelled as a blank; the columns of the two libraries that are blanks are not scored, as at level 3. No
deliverable asks for sample types.

**The count table is accepted under either filing of the two libraries.** An analyst who found the exchange
can report the two columns as the files are (each `sample_id` carries the read pairs of the FASTQ pair filed
under it; this is how the truth is stored) or by the sample each library really is (the same two columns the
other way round, which is what the investigator will want to read). Once the exchange is known, "the read
pairs of that sample" (instruction) and "of that library" (SOP item 2) can be read either way, no text says
which is meant, and none may be added at this level. The verifier therefore runs the level-3 logic twice,
once on the table as submitted and once with the values of the two columns exchanged (both together, nothing
else touched), and the better total is the score; `info.count_table_reading` in details.json records which
reading was used and both totals. A table in which only one of the two columns was moved is not one of the
two filings. If the blank's read pairs were moved and the specimen's were not, both columns hold the blank
library, the specimen's read pairs are nowhere, and the table loses `sample_counts` under both readings
(attack `one_column_moved_specimen_reads_lost`, 0.70). The opposite case cannot be told apart from a correct
table and is not rejected: if the specimen's read pairs were moved and the blank's were not, both columns
hold the specimen's library, which is right for the column that is scored under either reading, and the
column that would hold the blank library is, like every blank column since level 3, not scored (measured:
1.00 on seeds 1 and 2).

What the sheet's roles cost (`l6_reference`, the reference pipeline of levels 3-6 unchanged; seed 1 / seed
2): all five reagent ASVs are missed, because the two labelled blanks together are four fifths stool
community and the reagent ASVs are not enriched in them (planted counts: their share in the labelled blanks
is 0.69-1.03 / 0.97-1.26 times their share in the other libraries, the test asks for 2; with the real
roles 4.9-6.6 / 5.5-6.7), and 3 / 5 ASVs of the specimen are listed instead. `contaminants` drops to 0
(MCC -0.10 / -0.13) and `sample_counts` to 0 (lowest agreement 0.32 / 0.35: the three low-biomass specimens
keep their reagent reads, other specimens lose the listed ASVs; six / five of the eight specimen profiles
are below 0.80). Score 0.55 on both seeds; `asv_sequences`, `genus` and `mock_counts` are untouched. The
second solution of levels 3-6 on the same inputs lists nothing (with the qPCR value of a stool specimen in
the row of a blank library the absolute abundance of the reagent ASVs is not the same in blanks and
specimens) and scores 0.55 as well. On the planted counts of all 120 possible exchanges of twelve seeds the
reference test with the sheet's roles misses all five reagent ASVs (enrichment 0.58-1.52), lists one to
nine ASVs of the specimen, and the submission scores 0.55.

Routes that reach the truth: the two libraries are used in the roles their content supports (oracle: found
from the profiles alone, Bray-Curtis dissimilarity and richness; second solution: found from the mass
balance against the qPCR column, without sample types and without comparing profiles; both report the table
as the files are, probe `counts_filed_under_the_other_sample_id` is the same result filed by sample), or both
are left out
of the contaminant identification while their counts are still reported (probe `drop_the_odd_blank`: the
remaining blank is enough, the reagent ASVs are enriched 3.8-7.2 fold in it on the planted counts of the 120
exchanges, and the generator only plants an exchange for which this route returns exactly the planted
contaminants).

Not accepted (0.70, `sample_counts` lost, everything else right): removing the two libraries from the count
table as well (`drop_the_odd_blank_and_its_counts`; the profile of one specimen is gone under both readings).
The sheet-role submissions above gain nothing from the second reading: with the two columns exchanged the
scored column holds the blank library, so their better total is the one of the table as submitted (0.55).
"""


def _poisson_sf(k: int, mu: float) -> float:
    """P(X >= k) for X ~ Poisson(mu). (The generator host has numpy only; the solutions use scipy.)"""
    if k <= 0:
        return 1.0
    if k <= mu:                                   # the lower tail is the small one
        x = np.arange(0, k)
        logs = x * math.log(mu) - mu - np.array([math.lgamma(v + 1.0) for v in x])
        return max(0.0, 1.0 - float(np.exp(logs).sum()))
    x = np.arange(k, k + 400 + int(12 * math.sqrt(mu + k)))
    logs = x * math.log(mu) - mu - np.array([math.lgamma(v + 1.0) for v in x])
    top = float(logs.max())
    return float(math.exp(top) * np.exp(logs - top).sum()) if top > -700 else 0.0


def _enrichment(C: np.ndarray, neg: list[int], real: list[int]) -> np.ndarray:
    """The contaminant test of assets/amplicon_pipeline.py (levels 3-6) on a count table."""
    n_neg = C[:, neg].sum(axis=1)
    big_n = max(int(n_neg.sum()), 1)
    share_real = C[:, real].sum(axis=1) / max(int(C[:, real].sum()), 1)
    p = np.array([_poisson_sf(int(k), big_n * max(float(s), 1e-9)) for k, s in zip(n_neg, share_real)])
    in_all = (C[:, neg] > 0).all(axis=1)
    return in_all & (n_neg / big_n >= 2.0 * share_real) & (p < 1e-6)


def _constancy(C: np.ndarray, neg: list[int], real: list[int], copies: np.ndarray) -> np.ndarray:
    """The contaminant test of assets/solve_alt.py (levels 3-6) on a count table."""
    ab = C / np.maximum(C.sum(axis=0), 1) * copies
    return (C[:, neg] >= 2).all(axis=1) & (ab[:, neg].mean(axis=1) >= 0.3 * ab[:, real].mean(axis=1))


def _ramp(x: float, lo: float, hi: float) -> float:
    return max(0.0, min(1.0, (x - lo) / (hi - lo)))


def _score(C: np.ndarray, con: np.ndarray, scored: np.ndarray, types: list[str], flags: np.ndarray) -> dict:
    """What assets/score.py gives a submission that reports every planted sequence with its planted counts and
    genus and lists `flags` as contaminants, for the table as the files are (the second reading of the level-7
    verifier puts the blank library into the scored column and cannot be better). Used to select the instance
    only; the gate is the validator's run of the real scorer on the same submission (attack
    `sheet_roles_everything_else_right`)."""
    def agreement(j):
        obs_total = int(C[~flags, j].sum())
        true_total = int(C[~con, j].sum())
        return int(C[~con & ~flags, j].sum()) / max(obs_total, true_total, 1)

    study = min(agreement(j) for j, t in enumerate(types) if t == "study")
    mock = min(agreement(j) for j, t in enumerate(types) if t == "mock")
    tp, fn = int((flags & con & scored).sum()), int((~flags & con & scored).sum())
    fp, tn = int((flags & ~con & scored).sum()), int((~flags & ~con & scored).sum())
    den = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = (tp * tn - fp * fn) / den if den else 0.0
    score = (0.35 + 0.30 * _ramp(study, 0.60, 0.80) + 0.15 * _ramp(mcc, 0.50, 0.85) + 0.10
             + 0.10 * _ramp(mock, 0.60, 0.80))
    return {"score": round(score, 4), "min_study_agreement": round(study, 4), "mock_agreement": round(mock, 4),
            "mcc": round(mcc, 4), "reagent_missed": fn, "specimen_flagged": fp}


def _assess(C3: np.ndarray, total3: list[int], jb: int, jx: int, sheet_types: list[str], copies: np.ndarray,
            con: np.ndarray, scored: np.ndarray) -> dict:
    """Everything the generator needs to know about the exchange of blank jb and specimen jx, on the planted
    counts as filed (columns jb and jx of the level-3 counts exchanged)."""
    S = C3.shape[1]
    content = list(range(S))
    content[jb], content[jx] = jx, jb
    C = C3[:, content]
    real_types = [sheet_types[k] for k in content]
    pairs = [int(total3[k]) for k in content]
    rel = C / np.maximum(C.sum(axis=0), 1)
    bc = 1.0 - np.minimum(rel[:, :, None], rel[:, None, :]).sum(axis=0)
    rich = (rel >= RICH_SHARE).sum(axis=0)
    share = rel[con].sum(axis=0)
    amount = share * copies                                      # reagent share x qPCR copies of the row
    sheet_neg = [j for j in range(S) if sheet_types[j] == "negative_control"]
    fellows = [j for j in sheet_neg if j != jb]                  # the blank(s) that are blanks
    peers = [j for j in range(S) if sheet_types[j] == "study" and j != jx]   # the specimens that are specimens
    rest = [j for j in range(S) if j not in (jb, jx)]
    centre = float(np.median(amount[rest]))
    ev = {
        "read_pairs": {"labelled_blank": pairs[jb], "labelled_specimen": pairs[jx],
                       "other_blanks": [pairs[j] for j in fellows],
                       "other_specimens": [min(pairs[j] for j in peers), max(pairs[j] for j in peers)]},
        "asvs_above_0.5pct": {"labelled_blank": int(rich[jb]), "labelled_specimen": int(rich[jx]),
                              "other_blanks": [int(rich[j]) for j in fellows],
                              "other_specimens": [int(min(rich[j] for j in peers)), int(max(rich[j] for j in peers))]},
        "reagent_share": {"labelled_blank": round(float(share[jb]), 4), "labelled_specimen": round(float(share[jx]), 4),
                          "other_blanks": [round(float(share[j]), 4) for j in fellows],
                          "other_specimens": [round(float(min(share[j] for j in peers)), 4),
                                              round(float(max(share[j] for j in peers)), 4)]},
        "reagent_share_x_copies": {"labelled_blank": round(float(amount[jb]), 2),
                                   "labelled_specimen": round(float(amount[jx]), 1),
                                   "other_libraries": [round(float(amount[rest].min()), 1), round(centre, 1),
                                                       round(float(amount[rest].max()), 1)],
                                   "labelled_blank_factor_below": round(centre / max(float(amount[jb]), 1e-9), 1),
                                   "labelled_specimen_factor_above": round(float(amount[jx]) / centre, 1)},
        "bray_curtis": {"labelled_blank_to_other_blank": round(float(min(bc[jb, k] for k in fellows)), 3),
                        "labelled_blank_to_other_specimens": [round(float(min(bc[jb, j] for j in peers)), 3),
                                                              round(float(np.median([bc[jb, j] for j in peers])), 3)],
                        "labelled_specimen_to_other_blank": round(float(min(bc[jx, k] for k in fellows)), 3),
                        "other_specimens_to_other_blank_min": round(float(min(bc[j, k] for j in peers for k in fellows)), 3)},
    }
    real_neg = [j for j in range(S) if real_types[j] == "negative_control"]
    as_sheet = _enrichment(C, sheet_neg, [j for j in range(S) if j not in sheet_neg])
    as_real = _enrichment(C, real_neg, [j for j in range(S) if j not in real_neg])
    dropped = _enrichment(C, fellows, [j for j in range(S) if j not in sheet_neg and j != jx])
    alt_sheet = _constancy(C, sheet_neg, [j for j in range(S) if j not in sheet_neg], copies)

    def ratio(neg, real):                                        # enrichment of the reagent ASVs, range
        r = (C[con][:, neg].sum(axis=1) / max(int(C[:, neg].sum()), 1)) / np.maximum(
            C[con][:, real].sum(axis=1) / max(int(C[:, real].sum()), 1), 1e-12)
        return [round(float(r.min()), 2), round(float(r.max()), 2)]

    sheet_score = _score(C, con, scored, real_types, as_sheet)
    ev["sheet_roles"] = dict(sheet_score, reagent_enrichment=ratio(sheet_neg, [j for j in range(S) if j not in sheet_neg]),
                             second_method_flags=int(alt_sheet.sum()),
                             second_method_reagent_flagged=int((alt_sheet & con).sum()))
    ev["real_roles"] = {"exact": bool((as_real == con).all()),
                        "reagent_enrichment": ratio(real_neg, [j for j in range(S) if j not in real_neg])}
    ev["both_left_out"] = {"exact": bool((dropped == con).all()),
                           "reagent_enrichment": ratio(fellows, [j for j in range(S) if j not in sheet_neg and j != jx])}
    b_, a_ = ev["bray_curtis"], ev["reagent_share_x_copies"]
    checks = {
        "blank_is_as_rich_as_a_specimen": bool(rich[jb] >= RICH_BLANK_MIN * min(rich[j] for j in peers)
                                               and rich[jx] <= RICH_SPECIMEN_MAX * min(rich[j] for j in peers)),
        "profiles_separate": bool(b_["labelled_blank_to_other_blank"] >= BC_DIFF_MIN
                                  and b_["labelled_specimen_to_other_blank"] <= BC_SAME_MAX
                                  and b_["other_specimens_to_other_blank_min"] >= BC_NEXT_MIN
                                  and b_["labelled_blank_to_other_specimens"][0] < b_["labelled_blank_to_other_blank"]),
        "copies_relation_broken_for_the_two_only": bool(
            a_["labelled_blank_factor_below"] >= COPIES_FACTOR_MIN and a_["labelled_specimen_factor_above"] >= COPIES_FACTOR_MIN
            and amount[rest].max() <= COPIES_SCATTER_MAX * centre and amount[rest].min() >= centre / COPIES_SCATTER_MAX),
        "sheet_roles_fail": bool(sheet_score["score"] <= PASS - 0.15 + 1e-9),
        "real_roles_exact": ev["real_roles"]["exact"],
        "both_left_out_exact": ev["both_left_out"]["exact"],
    }
    return {"content": content, "counts": C, "real_types": real_types, "pairs": pairs, "evidence": ev,
            "checks": checks, "sheet_flags": as_sheet, "usable": all(checks.values())}


def _refile(src: Path, dst: Path) -> None:
    """The reads of `src` under the name `dst` (the gzip header stores the file name, mtime 0 as at level 3)."""
    with gzip.open(src, "rb") as fh:
        data = fh.read()
    dst.parent.mkdir(parents=True, exist_ok=True)
    with gzip.GzipFile(filename=str(dst), mode="wb", compresslevel=6, mtime=0) as fh:
        fh.write(data)


def _swap_columns(tsv: str, a: str, b: str) -> str:
    """A count table with the values of the columns a and b exchanged (header unchanged)."""
    lines = tsv.splitlines()
    head = lines[0].split("\t")
    ia, ib = head.index(a), head.index(b)
    out = [lines[0]]
    for line in lines[1:]:
        f = line.split("\t")
        f[ia], f[ib] = f[ib], f[ia]
        out.append("\t".join(f))
    return "\n".join(out) + "\n"


def _copy_column(tsv: str, src: str, dst: str) -> str:
    """A count table in which column dst carries the values of column src (src itself unchanged)."""
    lines = tsv.splitlines()
    head = lines[0].split("\t")
    i_src, i_dst = head.index(src), head.index(dst)
    out = [lines[0]]
    for line in lines[1:]:
        f = line.split("\t")
        f[i_dst] = f[i_src]
        out.append("\t".join(f))
    return "\n".join(out) + "\n"


def generate_l7(seed: int, level: int, workdir: Path):
    if level != 7:
        raise ValueError("l7.generate_l7 builds level 7 only")
    b = generate(seed, 3, workdir)                   # the level-3 batch of the seed, generated as at level 3
    truth = json.loads(b.tests_files["reference/truth.json"])
    sids = [s["id"] for s in truth["samples"]]
    sheet_types = [s["type"] for s in truth["samples"]]
    S = len(sids)
    head, *rows = [line.split("\t") for line in b.input_files["samples.tsv"].splitlines()]
    sheet = [dict(zip(head, r)) for r in rows]
    if [r["sample_id"] for r in sheet] != sids or [r["sample_type"] for r in sheet] != sheet_types:
        raise RuntimeError("level 7: sample sheet and planted truth of level 3 disagree")
    copies = np.array([float(r["qpcr_16s_copies_per_ul"]) for r in sheet])
    planted = truth["planted"]
    C3 = np.array([p["counts"] for p in planted], dtype=int)
    con = np.array([p["contaminant"] for p in planted], dtype=bool)
    scored = np.array([p["scored"] for p in planted], dtype=bool)
    cfrac = truth["info"]["contaminant_fraction"]

    # ---- which blank and which specimen: seeded draw among the exchanges that satisfy the four conditions ---
    drng = np.random.default_rng(77100 + seed)       # own stream: the level-3 batch of the seed is unchanged
    cands = [(jb, jx) for jb in range(S) if sheet_types[jb] == "negative_control"
             for jx in range(S) if sheet_types[jx] == "study" and cfrac[jx] <= MAX_REAGENT_SHARE]
    got, rejected = None, []
    for k in drng.permutation(len(cands)):
        jb, jx = cands[int(k)]
        got = _assess(C3, truth["total_pairs"], jb, jx, sheet_types, copies, con, scored)
        if got["usable"]:
            break
        rejected.append({"labelled_blank": sids[jb], "labelled_specimen": sids[jx],
                         "failed": [c for c, ok in got["checks"].items() if not ok]})
        got = None
    if got is None:
        raise RuntimeError(f"level 7, seed {seed}: no exchange of a blank and a specimen satisfies the conditions: "
                           f"{rejected}")
    content, C7, real_types, ev = got["content"], got["counts"], got["real_types"], got["evidence"]
    B, A = sids[jb], sids[jx]                        # B: labelled blank, holds the specimen; A: labelled specimen

    # ---- the two libraries filed under each other's name ---------------------------------------------------
    level3_file = {(s, m): b.input_files[f"fastq/{s}_{m}.fastq.gz"] for s in (A, B) for m in ("R1", "R2")}
    for src_id, dst_id in ((A, B), (B, A)):
        for mate in ("R1", "R2"):
            src = level3_file[(src_id, mate)]
            dst = workdir / "fastq_l7" / f"{dst_id}_{mate}.fastq.gz"
            if not isinstance(src, Path) or src.parent == dst.parent:
                raise RuntimeError("level 7: unexpected FASTQ source of the level-3 build")
            _refile(src, dst)                        # the level-3 files stay where they are
            b.input_files[f"fastq/{dst_id}_{mate}.fastq.gz"] = dst

    # ---- planted truth as filed ------------------------------------------------------------------------------
    for i, p in enumerate(planted):
        p["counts"] = [int(v) for v in C7[i]]
    truth["total_pairs"] = got["pairs"]
    for key in ("chimeric_pairs", "contaminant_fraction"):
        truth["info"][key] = [truth["info"][key][k] for k in content]
    truth["samples"] = [{"id": s, "type": real_types[j], "sheet_type": sheet_types[j]} for j, s in enumerate(sids)]
    truth["plants"] = {"control_swap": {
        "labelled_blank": B, "labelled_specimen": A, "holds": {B: f"library of specimen {A}", A: f"library of blank {B}"},
        "candidates": len(cands), "rejected_before": rejected, "checks": got["checks"], "evidence_on_planted_counts": ev}}
    b.tests_files["reference/truth.json"] = json.dumps(truth)
    # the level-3 verifier, unchanged, as a module; score.py runs it for both filings of the two columns
    b.tests_files["score_l3.py"] = b.tests_files["score.py"]
    b.tests_files["score.py"] = ASSETS / "score_l7.py"

    # ---- solutions ---------------------------------------------------------------------------------------------
    b.solution_files = {"solve.sh": ASSETS / "solve_l7.sh", "amplicon_pipeline_l7.py": ASSETS / "amplicon_pipeline_l7.py",
                        "solve_alt_l7.py": ASSETS / "solve_alt_l7.py"}

    # ---- battery: the level-3 battery on the level-7 data, plus the level-7 probes and attacks ----------------
    run3, run7 = "python3 /alegen_script/amplicon_pipeline.py ", "python3 /alegen_script/amplicon_pipeline_l7.py "
    l3_py = {"amplicon_pipeline.py": ASSETS / "amplicon_pipeline.py"}
    l7_py = {"amplicon_pipeline_l7.py": ASSETS / "amplicon_pipeline_l7.py"}
    ceil = PASS - 0.15

    def refiled(files: dict) -> dict:                # truth-derived submissions follow the files
        return dict(files, **{"asv_table.tsv": _swap_columns(files["asv_table.tsv"], A, B)})

    def rerun(script: str, name: str) -> str:        # the naive variants ignore exactly their own decision point
        if not script.startswith(run3):
            raise RuntimeError(f"level 7: unexpected script of the level-3 battery entry {name}")
        return run7 + script[len(run3):] + (" --roles sheet" if name == "naive_default_pipeline" else "")

    probes = []
    for p in b.probes:
        if p.name == "independent_solution":
            probes.append(Probe(p.name, p.description + "; the exchange is found from the mass balance against the "
                                "qPCR column (absolute abundance of constant-amount ASVs), without sample types and "
                                "without comparing profiles", p.min_score, script=p.script))
        elif p.script is None:
            probes.append(Probe(p.name, p.description, p.min_score, files=refiled(p.files)))
        else:
            probes.append(Probe(p.name, p.description + " (controls checked, the two libraries used in their real "
                                "roles)", p.min_score, script=rerun(p.script, p.name), script_files=l7_py))
    truth_files = next(p.files for p in probes if p.name == "planted_truth")
    table = {line.split("\t")[0]: [int(v) for v in line.split("\t")[1:]]
             for line in truth_files["asv_table.tsv"].splitlines()[1:]}
    fa = truth_files["asv.fasta"].split()
    id_of = {seq: name[1:] for name, seq in zip(fa[0::2], fa[1::2])}
    for p in planted:                                # the refiled probe is the planted truth of this level
        if (p["seq"] in id_of) != (sum(p["counts"]) > 0) or (p["seq"] in id_of and table[id_of[p["seq"]]] != p["counts"]):
            raise RuntimeError("level 7: the planted-truth probe does not match the planted counts as filed")
    sheet_listed = [id_of[p["seq"]] for i, p in enumerate(planted) if got["sheet_flags"][i] and p["seq"] in id_of]
    probes.append(Probe(
        "drop_the_odd_blank", "the odd blank is recognised, but the two libraries are not used in their real roles: "
        "both are left out of the contaminant identification (one blank against seven specimens and the mock); "
        "their counts are reported as for every file. A valid route: the count table is defined per file, no "
        "deliverable asks for sample types, and the remaining blank is enough to find the planted contaminants "
        "(a condition of the instance on the planted counts, re-measured here with the reference pipeline)", 0.99,
        script=run7 + "--roles drop", script_files=l7_py))
    probes.append(Probe(
        "counts_filed_under_the_other_sample_id", "the exchange is recognised and the count table is filed by the "
        "sample each library really is: planted truth with the right contaminants, the two columns the other way "
        "round (the level-3 table of the seed). Accepted like the table as the files are: the verifier scores "
        "both readings and takes the better one", 0.99,
        files=dict(truth_files, **{"asv_table.tsv": _swap_columns(truth_files["asv_table.tsv"], A, B)})))

    attacks = [
        Attack("l6_reference", "the reference pipeline of levels 3-6, unchanged, on the level-7 inputs: every decision "
               "point of the lower levels is handled, the roles of the libraries are taken from the sample sheet", ceil,
               script=run3.strip(), script_files=l3_py),
        Attack("l6_second_solution", "the independent second solution of levels 3-6, unchanged, on the level-7 inputs "
               "(absolute-abundance test with the sample types and qPCR values of the sheet)", ceil,
               script="python3 /alegen_script/solve_alt.py", script_files={"solve_alt.py": ASSETS / "solve_alt.py"}),
        Attack("sheet_roles_everything_else_right", "M_control_swap ignored and nothing else wrong: planted sequences, "
               "planted counts as filed and planted genus; contaminants as the reference test finds them on the "
               "planted counts with the roles of the sample sheet", ceil,
               files=dict(truth_files, **{"contaminants.txt": "asv_id\n" + "".join(i + "\n" for i in sheet_listed)})),
        Attack("drop_the_odd_blank_and_its_counts", "the odd blank is recognised and both libraries are removed from "
               "the batch: left out of the contaminant identification and their two columns of the count table set "
               "to zero. The profile of the specimen whose reads are in the blank's file is lost", ceil,
               script=run7 + "--roles drop_table", script_files=l7_py),
        Attack("one_column_moved_specimen_reads_lost", "inconsistent table: the exchange is recognised but only one "
               f"of the two columns is moved. Planted truth with the right contaminants, except that column {B} "
               f"carries the read pairs of the blank library (moved to the sample_id of the blank) while column {A} "
               "still carries them as filed; the read pairs of the specimen are in neither column. Neither of the "
               "two accepted filings", ceil, files=dict(truth_files, **{
                   "asv_table.tsv": _copy_column(truth_files["asv_table.tsv"], A, B)})),
    ]
    for a in b.attacks:
        if a.script is None:
            attacks.append(Attack(a.name, a.description, a.max_score, files=refiled(a.files)))
        else:
            note = (", roles from the sample sheet" if a.name == "naive_default_pipeline" else
                    " (controls checked, the two libraries used in their real roles)")
            attacks.append(Attack(a.name, a.description + note, a.max_score, script=rerun(a.script, a.name),
                                  script_files=l7_py))
    b.probes, b.attacks = probes, attacks

    # ---- task card -----------------------------------------------------------------------------------------------
    rp, ri, rs, ra, rb, sr = (ev[k] for k in ("read_pairs", "asvs_above_0.5pct", "reagent_share",
                                              "reagent_share_x_copies", "bray_curtis", "sheet_roles"))
    several = lambda v: "/".join(str(x) for x in v)          # one value per other blank             # noqa: E731
    span = lambda v: f"{v[0]:.2f}-{v[1]:.2f}"                # lowest-highest                         # noqa: E731
    b.traps = list(b.traps) + [Trap(
        "M_control_swap",
        f"The reads of extraction blank {B} and of study specimen {A} changed places before sequencing: the FASTQ "
        f"pair filed under {B} holds the library of the specimen ({rp['labelled_blank']} read pairs, a stool "
        f"community with {100 * rs['labelled_blank']:.1f} % reagent reads) and the pair filed under {A} holds the "
        f"library of the blank ({rp['labelled_specimen']} read pairs, {100 * rs['labelled_specimen']:.1f} % reagent "
        "reads). The sample sheet (sample_type, specimen, qPCR copies) is what the wet lab wrote down and is "
        "unchanged; no text mentions the exchange, sample identity or checking the controls. Contaminant "
        "identification with the roles of the sheet treats a community as a blank and a blank as a specimen. On the "
        f"planted counts the reference test then misses {sr['reagent_missed']} of the {int((con & scored).sum())} "
        f"reagent ASVs (their share in the labelled blanks is {span(sr['reagent_enrichment'])} times their share in "
        f"the other libraries; the test asks for 2; with the real roles "
        f"{span(ev['real_roles']['reagent_enrichment'])}) and lists {sr['specimen_flagged']} ASVs of the specimen; "
        "the absolute-abundance test of the second solution lists " + (
            "none at all" if not sr["second_method_flags"] else
            f"{sr['second_method_flags']} ASVs, {sr['second_method_reagent_flagged']} of them reagent ASVs")
        + ". Nothing fails. With everything else right such "
        f"a submission scores {sr['score']:.2f} (contaminant MCC {sr['mcc']:.2f}, lowest count agreement of a "
        f"specimen {sr['min_study_agreement']:.2f}).",
        "looking at the negative controls before they are used: read pairs per library, a bar plot or an ordination "
        "of the profiles, and the reagent background against the qPCR column. Measured on the planted counts as "
        f"filed. Read pairs: {B} (labelled blank) {rp['labelled_blank']}, the other blank "
        f"{several(rp['other_blanks'])}; {A} (labelled specimen) {rp['labelled_specimen']}, the other specimens "
        f"{rp['other_specimens'][0]}-{rp['other_specimens'][1]}. ASVs above 0.5 % of the library: {B} "
        f"{ri['labelled_blank']}, the other blank {several(ri['other_blanks'])}, the other specimens "
        f"{ri['other_specimens'][0]}-{ri['other_specimens'][1]}; {A} {ri['labelled_specimen']}. Share of the reads "
        f"in the ASVs that make up the other blank (the five reagent ASVs): {B} {rs['labelled_blank']:.3f}, the other "
        f"blank {several(f'{v:.3f}' for v in rs['other_blanks'])}; {A} {rs['labelled_specimen']:.3f}, the other specimens "
        f"{rs['other_specimens'][0]:.3f}-{rs['other_specimens'][1]:.3f}. Bray-Curtis dissimilarity of the relative "
        f"profiles: {B} to the other blank {rb['labelled_blank_to_other_blank']:.2f}, to the nearest other specimen "
        f"{rb['labelled_blank_to_other_specimens'][0]:.2f} (median {rb['labelled_blank_to_other_specimens'][1]:.2f}); "
        f"{A} to the other blank {rb['labelled_specimen_to_other_blank']:.2f}, while no other specimen comes nearer "
        f"to it than {rb['other_specimens_to_other_blank_min']:.2f}. Copies-versus-share relation (reagent share x "
        f"qPCR copies per microlitre of the row, which is the same amount for every extract): "
        f"{ra['other_libraries'][0]:.0f}-{ra['other_libraries'][2]:.0f} (median {ra['other_libraries'][1]:.0f}) in "
        f"the nine other libraries, {ra['labelled_specimen']:.0f} for {A} ({ra['labelled_specimen_factor_above']:.0f} "
        f"times above) and {ra['labelled_blank']:.1f} for {B} ({ra['labelled_blank_factor_below']:.0f} times below); "
        "with the two rows exchanged both lie in the range of the others.",
        "Verification / missing QC")]
    b.judge_explanation_md = b.judge_explanation_md.rstrip() + L7_JUDGE
    b.capabilities = list(b.capabilities) + ["checking negative controls against the sample sheet before using them"]
    b.oracle_steps = 11
    b.level = 7
    b.task_id = f"amplicon-16s-s{seed:04d}-l7"
    return b
