"""Family: tumour/normal pairs from FASTQ to filtered somatic SNV/indel calls.

Truth by construction: the simulator fixes the germline genotypes of three patients, the somatic
variants of their tumours (with allele fractions), library artefacts and the read pair every allele
is written to, and only then emits paired FASTQ from a real hg38 chr22 slice. The hidden reference
is the list of planted variants; no aligner or caller is run to obtain it.

Decision points (level 2 adds the first two, level 3 the third):
  sample_swap            the normal libraries of two patients are exchanged in the sample sheet
  duplicated_segment     the panel reference holds one exact intra-slice duplication; reads from
                         its interior cannot be placed, and a variant of one copy shows up at both
  orientation_artifacts  one or two tumour libraries carry damage-type substitutions that exist in
                         one read-pair orientation only (F1R2 or F2R1), at 8-30% allele fraction

Layout: `sim.py` simulator, `assets/` oracle (`somatic_pipeline.py`), independent solution
(`solve_alt.py`) and verifier, `dev/` development harness that is not shipped in tasks.
"""
from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np

from ... import dockerlib
from ...spec import Attack, Probe, TaskBuild, Trap
from . import sim

ASSETS = Path(__file__).parent / "assets"
CONDA = ["bwa=0.7.19", "minimap2=2.30", "samtools=1.25", "bcftools=1.25", "pysam=0.24.1"]
PASS = 0.9
DEPTH_T, DEPTH_N = 180, 90

INSTRUCTION = """\
You are the bioinformatician of a cancer-genomics core facility. Tumour and normal DNA of three
patients was sequenced on a small hybrid-capture panel. The clinical team needs, for every patient,
the tumour/normal pair and the somatic variants of the tumour.

Task directory:
- `/task`

Inputs (read-only) under `/task/input`:
- `fastq/<sample_id>_R1.fastq.gz`, `fastq/<sample_id>_R2.fastq.gz` - paired-end reads (2 x 100 nt),
  one pair of files per library (six libraries)
- `manifest.tsv` - sample sheet exported from the LIMS: `patient`, `sample_type` (`tumour` or
  `normal`), `sample_id`, `fastq_1`, `fastq_2`
- `reference/ref.fa` (+ `.fai`) - the panel reference, one contig (`{contig}`); no aligner index is
  provided
- `targets.bed` - the capture targets on the panel reference (0-based, half-open)
- `core_notes.md` - hand-over notes of the sequencing core

Available software (on `PATH`): `bwa`, `minimap2`, `samtools`, `bcftools` (with `bgzip`, `tabix`),
Python 3.12 with `pysam`, numpy, scipy and pandas. There is no internet access.

Your task:
1. Align the reads of every library to the panel reference.
2. Establish, for every patient, which normal library was prepared from the same individual as the
   patient's tumour library.
3. Call the somatic single-nucleotide variants and short insertions/deletions of every tumour
   inside the capture targets.

Required outputs under `/task/output`:
- `pairing.tsv` - tab-separated with the columns `patient`, `tumour_sample`, `normal_sample`; one
  row per patient of the manifest. `tumour_sample` is the `sample_id` of the patient's tumour
  library as listed in the manifest; `normal_sample` is the `sample_id` of the normal library that
  comes from the same individual as that tumour library.
- `somatic/<patient>.vcf.gz` - one gzip- or bgzip-compressed VCF per patient holding the somatic
  variants of that patient's tumour library.
  - A somatic variant is a single-nucleotide variant or short insertion/deletion that is present in
    tumour cells of the patient and absent from the patient's germline. Germline variants, and
    differences from the reference that arise during library preparation, sequencing or read
    alignment, are not somatic variants.
  - Reported variants have `FILTER` = `PASS`; records with any other FILTER value are not part of
    the deliverable.
  - `CHROM`, `POS`, `REF`, `ALT` follow the VCF specification relative to `reference/ref.fa`; one
    alternate allele per record; insertions and deletions left-aligned and parsimonious (the
    representation `bcftools norm` produces).
  - Every record carries `AF=<number>` in the `INFO` column: the fraction (0-1) of the tumour
    library's DNA fragments covering the site that carry the alternate allele.
  - The reporting range of the assay starts at an allele fraction of 5%; variants below it may be
    left out.
- `summary.json` - `{{"patients": {{"<patient>": {{"tumour_sample": "<sample_id>", "normal_sample":
  "<sample_id>", "n_snv": <int>, "n_indel": <int>}}, ...}}}}` where `n_snv` and `n_indel` are the
  numbers of `PASS` single-nucleotide and insertion/deletion records in that patient's VCF and the
  two sample ids repeat `pairing.tsv`.

Constraints:
- Do not modify anything under `/task/input`; write deliverables only under `/task/output`.
- Extra columns, INFO fields and sample columns are allowed.
- Do not ask for confirmation. Execute directly.
"""

JUDGE = """\
**Where the truth comes from.** Before any read exists the generator fixes (a) the germline
genotype of each patient at 110 population SNP sites, (b) the somatic SNVs and indels of each
tumour with their allele fractions, (c) artefact sites and the read-pair orientation that shows
them, (d) which library is written to which sample id, and then (e) which sequenced fragment
carries which allele. `reference/truth.json` holds the planted variants and, per variant, the
realised fraction of tumour fragments that carry it. No aligner or caller output is used.

**Scored set and gray zone.** A planted somatic variant is *scored* if the simulated tumour library
holds at least {min_alt} alt-carrying fragments at the site, a realised allele fraction of at least
{min_af}, and at least max({min_minor}, {min_frac} x alt fragments) alt fragments in the rarer
read-pair orientation. Planted variants below that floor, and planted variants in the interior of
the duplicated segment (every supporting read pair maps equally well to both copies), are *gray*:
calling them is neither rewarded nor penalised. Every other PASS call inside the targets is a false
positive; that includes the patient's germline variants, orientation artefacts and the mirror
position of a duplicated-segment variant in the other copy. Calls outside the targets are ignored.

**Hard gates (score 0).** Missing/unparseable file; a patient missing or duplicated in
`pairing.tsv`; a VCF record on another contig, outside the reference, or whose REF does not match
`ref.fa`; `summary.json` without a `patients` object.

**Components.** `recall` and `precision` are computed per patient and the **minimum** over the
three patients is used, so a patient whose pair was left wrong cannot be averaged away.

| weight | component | metric | full credit | zero |
|---|---|---|---|---|
| 0.20 | pairing | every patient's (tumour, normal) pair equals the planted one | all | any wrong |
| 0.30 | somatic_precision | min over patients of TP / (TP + FP) | >= 0.95 | <= 0.85 |
| 0.30 | somatic_recall | min over patients of TP / scored variants; multiplied by ramp(precision, 0.50, 0.85) | >= 0.95 | <= 0.85 |
| 0.10 | allele_fraction | worst patient's median abs(AF - realised fragment fraction) over true positives | <= 0.03 | >= 0.08 |
| 0.10 | summary | fraction of patients whose summary entry matches their own VCF counts and pairing row | all | none |

pairing, somatic_precision, allele_fraction and summary are multiplied by ramp(min recall, 0.2,
0.5): a deliverable without real calls earns nothing.

**Why these tolerances.** See `docs/families/somatic_tn.md` for the measured values of the oracle
(bwa + pysam fragment counting), the independent solution (minimap2 + bcftools mpileup on
orientation-split BAMs) and every trap-ignoring analysis on two seeds.

**What is deliberately not scored.** QUAL, genotype columns, variants below the scoring floor,
anything outside the targets, the numeric concordance values behind the pairing.
"""


def _vcf(contig: str, length: int, tumour: str, recs: list[tuple[int, str, str, float]]) -> bytes:
    lines = ["##fileformat=VCFv4.2", f"##contig=<ID={contig},length={length}>",
             '##INFO=<ID=AF,Number=1,Type=Float,Description="Tumour allele fraction">',
             f"##tumour_sample={tumour}",
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    for pos0, ref, alt, af in sorted(recs):
        lines.append(f"{contig}\t{pos0 + 1}\t.\t{ref}\t{alt}\t.\tPASS\tAF={af:.4f}")
    return gzip.compress(("\n".join(lines) + "\n").encode(), mtime=0)


def _notes(level: int, samples: dict, plan: dict, ref: dict, order: list[str], pid: dict) -> str:
    n_t = len(ref["targets"])
    bp = sum(e - s for s, e in ref["targets"])
    lines = [
        "# Sequencing core hand-over: run 214, panel CGP-22", "",
        "## Assay", "",
        f"- Hybrid-capture panel CGP-22: {n_t} target regions, {bp / 1000:.1f} kb in total (`targets.bed`). The panel",
        f"  reference `reference/ref.fa` (contig `{sim.CONTIG}`) is the hg38 chr22 slice the baits were designed on.",
        "- Libraries: enzymatic end repair / A-tailing, ligation of dual-indexed adapters, PCR amplification,",
        "  one pooled capture reaction for all six libraries.",
        "- Sequencing: NovaSeq 6000, flow cell HKF3WDSX7 lane 2, 2 x 100 cycles. Demultiplexed with bcl2fastq;",
        "  adapters were trimmed by the instrument software. FASTQ files are otherwise as produced.",
        f"- Requested mean depth on target: tumour {DEPTH_T}x, normal {DEPTH_N}x.", "",
        "## Libraries", "",
        "| sample_id | specimen | DNA input (ng) | shearing | PCR cycles | i7+i5 index |",
        "|---|---|---|---|---|---|",
    ]
    for sid in order:
        s = samples[sid]
        lines.append(f"| {sid} | {s['specimen']} | {s['input_ng']} | {s['shear']} | {s['pcr']} | {s['index']} |")
    lines += ["", "## Remarks", ""]
    lines.append("- Shearing: Covaris LE220, low-EDTA TE, standard setting, one batch.")
    if plan["damaged"]:
        lines.append("- FFPE sections were deparaffinised and extracted with the FFPE DNA kit of the pathology lab.")
    lines.append("- Tumour cell content estimated by the pathologist on the adjacent H&E section: "
                 + ", ".join(f"{pid[p]} {int(round(plan['patients'][p]['purity'] * 10)) * 10}%" for p in sorted(plan["patients"]))
                 + ".")
    lines.append("- All libraries passed the core's yield and fragment-size checks (mean insert 260-280 bp).")
    return "\n".join(lines) + "\n"


def generate(seed: int, level: int, workdir: Path) -> TaskBuild:
    workdir = Path(workdir)
    fq_dir = workdir / "fastq"
    fq_dir.mkdir(parents=True, exist_ok=True)
    for old in fq_dir.glob("*.fastq.gz"):
        old.unlink()
    ref = sim.build_reference(seed, with_dup=level >= 2)
    plan = sim.plant(ref, seed, level)
    seq = ref["seq"]
    rng = np.random.default_rng(72_000 + seed)
    patients = sorted(plan["patients"])

    # ---- sample ids, indices, sample sheet -------------------------------------------------
    nums = rng.permutation(np.arange(1100, 9900))[:6]
    case = rng.permutation(np.arange(203, 988))[:3]
    pid = {p: f"PT-{int(c):04d}" for p, c in zip(patients, sorted(case))}
    samples: dict[str, dict] = {}
    tum, nor = {}, {}
    k = 0
    for p in patients:
        for typ in ("tumour", "normal"):
            sid = f"CG24-{int(nums[k]):04d}"
            k += 1
            idx = "".join("ACGT"[i] for i in rng.integers(0, 4, 8)) + "+" + "".join("ACGT"[i] for i in rng.integers(0, 4, 8))
            damaged = typ == "tumour" and p in plan["damaged"]
            if typ == "normal":
                specimen = "peripheral blood (EDTA)"
            elif damaged:
                specimen = f"FFPE section (block from {int(rng.integers(2017, 2022))})"
            else:
                specimen = "fresh-frozen biopsy"
            samples[sid] = {"patient": p, "type": typ, "index": idx, "specimen": specimen,
                            "input_ng": int(rng.choice([50, 80, 100, 150, 200])) if typ == "tumour" else 200,
                            "shear": "batch A",
                            "pcr": 10 if "FFPE" in specimen else 8}
            (tum if typ == "tumour" else nor)[p] = sid
    swap = sorted(rng.permutation(patients)[:2].tolist()) if level >= 2 else []
    sheet_nor = dict(nor)
    if swap:
        a, b = swap
        sheet_nor[a], sheet_nor[b] = nor[b], nor[a]
    order = [str(s) for s in rng.permutation(sorted(samples))]
    manifest = "patient\tsample_type\tsample_id\tfastq_1\tfastq_2\n"
    for p in patients:
        for typ, sid in (("tumour", tum[p]), ("normal", sheet_nor[p])):
            manifest += f"{pid[p]}\t{typ}\t{sid}\tfastq/{sid}_R1.fastq.gz\tfastq/{sid}_R2.fastq.gz\n"

    # ---- reads -----------------------------------------------------------------------------
    eff = [float(x) for x in rng.uniform(0.85, 1.15, size=len(ref["targets"]))]
    truth_p = {}
    fastq_files = {}
    n_reads = {}
    for p in patients:
        for typ, sid in (("tumour", tum[p]), ("normal", nor[p])):
            depth = (DEPTH_T if typ == "tumour" else DEPTH_N) * float(rng.uniform(0.95, 1.05))
            r1, r2, st = sim.simulate_sample(ref, plan["patients"][p], typ == "tumour", depth, eff, rng)
            f1, f2 = fq_dir / f"{sid}_R1.fastq.gz", fq_dir / f"{sid}_R2.fastq.gz"
            n_reads[sid] = sim.write_fastq_pair(f1, f2, r1, r2, samples[sid]["index"], rng)
            fastq_files[f"fastq/{f1.name}"] = f1
            fastq_files[f"fastq/{f2.name}"] = f2
            if typ != "tumour":
                continue
            scored, gray, mirror = [], [], []
            for j, v in enumerate(plan["patients"][p]["somatic"]):
                alt = int(st["som_alt"][j].sum())
                minor = int(st["som_alt"][j].min())
                cov = int(st["som_cov"][j])
                af = alt / cov if cov else 0.0
                rec = {"pos": v["pos0"] + 1, "ref": v["ref"], "alt": v["alt"], "af": round(af, 5), "kind": v["kind"],
                       "cls": v["cls"], "planted_vaf": round(v["vaf"], 4), "alt_fragments": alt,
                       "minor_orientation_fragments": minor, "fragments": cov}
                ok = (v["cls"] != "dup" and alt >= sim.MIN_ALT_FRAG and af >= sim.MIN_AF
                      and minor >= max(sim.MIN_MINOR_FRAG, sim.MIN_MINOR_FRAC * alt))
                (scored if ok else gray).append(rec)
                if v["cls"] == "dup":
                    mirror.append({"pos": v["mirror0"] + 1, "ref": seq[v["mirror0"]], "alt": v["alt"]})
            arts = [{"pos": a["pos0"] + 1, "ref": a["ref"], "alt": a["alt"], "orientation": ["F1R2", "F2R1"][a["orient"]],
                     "chemistry": a["chemistry"],
                     "alt_fragments": int(st["art_alt"][j]), "fragments": int(st["art_cov"][j])}
                    for j, a in enumerate(plan["patients"][p]["artefacts"])]
            germ = [{"pos": x + 1, "ref": seq[x], "alt": a, "gt": f"{h[0]}|{h[1]}"}
                    for x, (a, h) in sorted(plan["patients"][p]["germline"].items())]
            truth_p[pid[p]] = {"tumour_sample": tum[p], "normal_sample": nor[p], "sheet_normal": sheet_nor[p],
                               "purity": round(plan["patients"][p]["purity"], 4),
                               "scored": scored, "gray": gray, "mirror": mirror, "artefacts": arts, "germline": germ}

    # the battery relies on these properties; fail loudly at generation instead of shipping a weak instance
    for q, t in truth_p.items():
        assert len(t["scored"]) >= 20 and sum(1 for r in t["scored"] if r["kind"] != "snv") >= 3, (q, len(t["scored"]))
        if level >= 2:
            assert len(t["mirror"]) >= 9 and len(t["mirror"]) / (len(t["mirror"]) + len(t["scored"])) >= 0.19, q
        strong = sum(1 for r in t["artefacts"] if r["alt_fragments"] >= 10 and r["alt_fragments"] / r["fragments"] >= 0.06)
        assert not t["artefacts"] or strong / (strong + len(t["scored"])) >= 0.22, (q, strong)
    if swap:
        # germline variants a patient does not share with the owner of the normal the sheet assigns to it
        a, b = swap
        priv = [len(set(plan["patients"][x]["germline"]) - set(plan["patients"][y]["germline"])) for x, y in ((a, b), (b, a))]
        assert max(priv) >= 12, priv

    truth = {"contig": sim.CONTIG, "length": len(seq), "targets": [list(t) for t in ref["targets"]],
             "patients": truth_p, "duplication": ref["dup"], "swapped_patients": [pid[p] for p in swap],
             "damage_mode": plan["damage_mode"], "damaged_patients": [pid[p] for p in plan["damaged"]],
             "floor": {"min_alt_fragments": sim.MIN_ALT_FRAG, "min_af": sim.MIN_AF,
                       "min_minor_fragments": sim.MIN_MINOR_FRAG, "min_minor_fraction": sim.MIN_MINOR_FRAC},
             "reads": n_reads, "reference_source": ref["source"]}
    fa, fai = sim.fasta(sim.CONTIG, seq)
    bed = "".join(f"{sim.CONTIG}\t{s}\t{e}\tCGP22_T{i + 1:02d}\n" for i, (s, e) in enumerate(ref["targets"]))
    notes = _notes(level, samples, plan, ref, order, pid)

    # ---- truth-derived probes / attacks ------------------------------------------------------
    prng = np.random.default_rng(seed + 991)
    L = len(seq)

    def decoys(pat: str, n: int) -> list[tuple[int, str, str, float]]:
        """n false calls inside targets, away from every planted site of the patient."""
        t = truth_p[pat]
        used = {r["pos"] - 1 for key in ("scored", "gray", "mirror", "artefacts", "germline") for r in t[key]}
        out = []
        while len(out) < n:
            s, e = ref["targets"][int(prng.integers(0, len(ref["targets"])))]
            x = int(prng.integers(s + 50, e - 50))
            if any(abs(x - u) < 5 for u in used):
                continue
            used.add(x)
            out.append((x, seq[x], "ACGT"[("ACGT".index(seq[x]) + 1 + int(prng.integers(0, 3))) % 4], 0.2))
        return out

    def submission(pairing: dict, calls: dict) -> dict:
        files = {"pairing.tsv": "patient\ttumour_sample\tnormal_sample\n" + "".join(
            f"{q}\t{truth_p[q]['tumour_sample']}\t{pairing[q]}\n" for q in truth_p)}
        summ = {}
        for q in truth_p:
            files[f"somatic/{q}.vcf.gz"] = _vcf(sim.CONTIG, L, truth_p[q]["tumour_sample"], calls[q])
            summ[q] = {"tumour_sample": truth_p[q]["tumour_sample"], "normal_sample": pairing[q],
                       "n_snv": sum(1 for c in calls[q] if len(c[1]) == len(c[2])),
                       "n_indel": sum(1 for c in calls[q] if len(c[1]) != len(c[2]))}
        files["summary.json"] = json.dumps({"patients": summ}, indent=1)
        return files

    true_pair = {q: truth_p[q]["normal_sample"] for q in truth_p}
    sheet_pair = {q: truth_p[q]["sheet_normal"] for q in truth_p}
    true_calls = {q: [(r["pos"] - 1, r["ref"], r["alt"], r["af"]) for r in truth_p[q]["scored"]] for q in truth_p}
    tol_calls, fp_calls, fn_calls = {}, {}, {}
    for q in truth_p:
        c = list(true_calls[q])
        drop = int(prng.integers(0, len(c)))
        tol_calls[q] = [(x, r, a, min(1.0, max(0.0, af + float(prng.uniform(-0.02, 0.02)))))
                        for i, (x, r, a, af) in enumerate(c) if i != drop] + decoys(q, 1)
        fp_calls[q] = c + decoys(q, int(np.ceil(len(c) * 0.30)))
        keep = set(prng.permutation(len(c))[:int(len(c) * 0.75)].tolist())
        fn_calls[q] = [r for i, r in enumerate(c) if i in keep]
    pipe = {"somatic_pipeline.py": ASSETS / "somatic_pipeline.py"}
    run = "python3 /alegen_script/somatic_pipeline.py "

    probes = [
        Probe("planted_truth", "submission written directly from the planted truth", 0.999,
              files=submission(true_pair, true_calls)),
        Probe("within_tolerance", "per patient one scored variant missed, one false call, AF noise +-0.02", 0.99,
              files=submission(true_pair, tol_calls)),
        Probe("independent_solution", "second tool chain: minimap2 (no bwa), bcftools mpileup allele depths on "
              "orientation-split BAMs (no pysam), genotype discordance from bcftools call for the pairing",
              0.99, script="SOLUTION=2 bash /solution/solve.sh"),
    ]
    attacks = [
        Attack("sheet_pairing_no_calls", "prior-only answer: pairing copied from the sample sheet, VCFs without records",
               0.10, files=submission(sheet_pair, {q: [] for q in truth_p})),
        Attack("false_calls_beyond_tolerance", "planted truth plus 30% false calls per patient", PASS - 0.15,
               files=submission(true_pair, fp_calls)),
        Attack("missed_calls_beyond_tolerance", "planted truth with a quarter of the variants missing", PASS - 0.15,
               files=submission(true_pair, fn_calls)),
        Attack("naive_tumour_only", "correct alignment and filters, but no matched normal is used: germline variants "
               "are reported as somatic", 0.50, script=run + "--no-normal", script_files=pipe),
    ]
    traps = []
    if level == 1:
        probes.append(Probe("default_pipeline_is_correct_at_level1", "sample sheet pairing, no mapping-quality and no "
                            "orientation filter: without planted obstacles the default analysis is right", 0.99,
                            script=run + "--pairing sheet --min-mapq 0 --no-orientation", script_files=pipe))
    if level >= 2:
        a, b = [pid[p] for p in swap]
        attacks += [
            Attack("naive_sheet_pairing", "all filters applied, but the pairs of the sample sheet are used as given",
                   PASS - 0.15, script=run + "--pairing sheet", script_files=pipe),
            Attack("naive_no_mapq_filter", "verified pairs and orientation filter, but reads of any mapping quality "
                   "are counted", PASS - 0.15, script=run + "--min-mapq 0", script_files=pipe),
            Attack("naive_default_pipeline", "sample sheet pairs, no mapping-quality filter, no orientation filter",
                   0.40, script=run + "--pairing sheet --min-mapq 0 --no-orientation", script_files=pipe),
        ]
        traps += [
            Trap("sample_swap",
                 f"The normal libraries of {a} and {b} are exchanged in manifest.tsv. With the sheet pairs, every "
                 "germline variant private to the patient is called somatic in two of three patients.",
                 "pairwise germline genotype concordance between the six libraries: each tumour agrees with exactly "
                 "one normal at ~100% of variable sites and with the others at well under 60%; for two patients the "
                 "matching normal is not the one in the sheet. Also visible downstream: dozens of 'somatic' calls at "
                 "allele fraction ~0.5 and ~1.0 in two patients",
                 "Understanding / unverified assumption about inputs"),
            Trap("duplicated_segment",
                 f"A {ref['dup']['len']} bp segment occurs twice in the panel reference (two targets). Read pairs from "
                 "its interior align equally well to both copies; a tumour variant of one copy appears at both, with "
                 "half the allele fraction, supported by mapping-quality-0 reads only.",
                 "mapping-quality distribution per target (two targets consist almost entirely of MAPQ 0 pairs); the "
                 "candidate calls there come in pairs at identical offsets in the two targets; aligning the reference "
                 "against itself shows the duplication",
                 "Approach / missing QC of alignment evidence"),
        ]
    if level >= 3:
        attacks.append(Attack("naive_no_orientation_filter", "verified pairs and mapping-quality filter, but no check of "
                              "the read-pair orientation of the supporting reads", PASS - 0.15,
                              script=run + "--no-orientation", script_files=pipe))
        attacks.append(Attack("naive_strand_bias_only", "verified pairs and mapping-quality filter; supporting reads "
                              "are required on both alignment strands (classic strand-bias filter) instead of in both "
                              "read-pair orientations", PASS - 0.15,
                              script=run + "--strand-instead-of-orientation", script_files=pipe))
        attacks.append(Attack("shortcut_spectrum_filter", "verified pairs and mapping-quality filter; instead of "
                              "looking at the supporting reads, every C>T, G>A, G>T and C>A call below 30% allele "
                              "fraction is dropped in the tumours the notes list as FFPE (removes the artefacts and "
                              "the true variants of those classes alike)",
                              PASS - 0.15, script=run + "--no-orientation --spectrum-filter 0.30", script_files=pipe))
        sub = "C>T / G>A (deamination) and G>T / C>A (oxidation)"
        traps.append(Trap(
            "orientation_artifacts",
            f"The tumour libraries of {', '.join(pid[p] for p in plan['damaged'])} carry 14-18 {sub} artefact sites "
            "each at 8-30% allele fraction. Every supporting fragment has the same read-pair orientation (F1R2 for "
            "one substitution, F2R1 for its complement); forward/reverse strand counts are balanced, so a plain "
            "strand-bias filter does not see them.",
            "the unfiltered calls of the affected tumours contain an excess of two substitution types that the "
            "other tumours do not show; tabulating alt-supporting reads by read number and strand "
            "shows 100% F1R2 or 100% F2R1 at those sites while true variants are balanced; core_notes.md lists the "
            "specimen type of each library without comment",
            "Understanding / domain knowledge gap"))

    in_bytes = sum(f.stat().st_size for f in fastq_files.values())
    return TaskBuild(
        task_id=f"somatic-tn-s{seed:04d}-l{level}", family="somatic_tn", seed=seed, level=level,
        title="Tumour/normal pairs: FASTQ to filtered somatic variants",
        summary="Align six panel libraries of three cancer patients, establish which normal belongs to which tumour, "
                "and deliver filtered somatic SNV/indel calls with allele fractions per patient.",
        subdomain="3.2 Genomics & Sequence Analysis",
        software=["bwa", "minimap2", "samtools", "bcftools", "pysam", "Python"],
        instruction_md=INSTRUCTION.format(contig=sim.CONTIG),
        agent_dockerfile=dockerlib.agent_dockerfile(conda=CONDA),
        verifier_dockerfile=dockerlib.verifier_dockerfile(pip=[]),
        input_files={**fastq_files, "manifest.tsv": manifest, "reference/ref.fa": fa, "reference/ref.fa.fai": fai,
                     "targets.bed": bed, "core_notes.md": notes},
        solution_files={"solve.sh": ASSETS / "solve.sh", "somatic_pipeline.py": ASSETS / "somatic_pipeline.py",
                        "solve_alt.py": ASSETS / "solve_alt.py"},
        tests_files={"score.py": ASSETS / "score.py", "reference/truth.json": json.dumps(truth),
                     "reference/ref.fa": fa},
        traps=traps, attacks=attacks, probes=probes,
        capabilities=["sample identity QC before analysis", "read-level artefact diagnosis (orientation, mapping quality)",
                      "long tool chain from raw reads to filtered calls", "tumour/normal comparison without a turnkey caller",
                      "exact output contracts (VCF normalisation)"],
        judge_explanation_md=JUDGE.format(min_alt=sim.MIN_ALT_FRAG, min_af=sim.MIN_AF, min_minor=sim.MIN_MINOR_FRAG,
                                          min_frac=sim.MIN_MINOR_FRAC),
        pass_threshold=PASS, oracle_steps=33, oracle_runtime_min=0.4, cpus=4, memory_mb=8000, storage_mb=10000,
        extra_metadata={"reference_source": sim.REF_URL, "reference_sha256": sim.REF_SHA256,
                        "fastq_bytes": in_bytes},
    )
