"""Family: pooled CRISPR-Cas9 dropout screen, FASTQ -> guide counts -> gene fitness calls.

Truth-by-construction: the simulator decides which genes have a fitness effect and how large it
is, then emits reads. Gene identities are blinded (TGT0001...), so prior knowledge about real
essential genes cannot substitute for analysing the data.

Decision points (level 2 adds the first two, level 3 adds the third):
  staggered_amplicon     spacer position varies read to read (P5 stagger primers)
  read2_only_sample      one sample was delivered as Read 2 -> reads are reverse complements
  hit_enriched_library   validation library where most genes are true hits: median normalisation
                         and rank-based (RRA) testing both assume the opposite
"""
from __future__ import annotations

import gzip
import io
import json
from pathlib import Path

import numpy as np

from ... import dockerlib
from ...spec import Attack, Probe, TaskBuild, Trap

ASSETS = Path(__file__).parent / "assets"
PREFIX = "TCTTGTGGAAAGGACGAAACACCG"
SCAFFOLD = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"
READLEN = 75
COMP = str.maketrans("ACGT", "TGCA")


def _rc(s: str) -> str:
    return s.translate(COMP)[::-1]


def _guides(rng: np.random.Generator, n: int) -> list[str]:
    """Random 20-mers: GC 35-70%, no homopolymer > 4, pairwise Hamming distance >= 4."""
    out: list[str] = []
    arr = np.zeros((0, 20), dtype=np.int8)
    bases = np.array(list("ACGT"))
    while len(out) < n:
        cand = rng.integers(0, 4, size=20).astype(np.int8)
        s = "".join(bases[cand])
        gc = (s.count("G") + s.count("C")) / 20
        if not 0.35 <= gc <= 0.70 or any(b * 5 in s for b in "ACGT"):
            continue
        if s[:6] in PREFIX or s[-6:] in SCAFFOLD[:30] or "GTTTTAGA" in s:
            continue
        if len(arr) and (arr != cand).sum(axis=1).min() < 4:
            continue
        out.append(s)
        arr = np.vstack([arr, cand])
    return out


def _write_fastq(path: Path, origin: np.ndarray, fwd: list[str], rev: list[str], reverse: bool,
                 stagger: bool, rng: np.random.Generator, run: str) -> None:
    bases = "ACGT"
    n = len(origin)
    stag = rng.integers(0, 8, size=n) if stagger else np.zeros(n, dtype=int)
    nerr = rng.binomial(READLEN, 0.002, size=n)
    tiles = rng.integers(1101, 1140, size=n)
    xs = rng.integers(1000, 30000, size=n)
    ys = rng.integers(1000, 30000, size=n)
    rand = rng.integers(0, 4, size=(n, 8))
    hdr = "@NB551608:214:HV2KJBGXK:1" if run == "A" else "@A01789:96:HFG7TDSX5:3"
    mate = "2" if reverse else "1"
    qual = "I" * (READLEN - 6) + "FFFFFF"
    # mtime=0 and no stored file name: the same seed must give byte-identical archives
    with open(path, "wb") as raw, gzip.GzipFile(filename="", fileobj=raw, mode="wb", compresslevel=6, mtime=0) as gz, \
            io.TextIOWrapper(gz, encoding="ascii", newline="\n") as fh:
        for i in range(n):
            g = origin[i]
            if g < 0:   # primer-dimer / unrelated read
                seq = "".join(bases[b] for b in rng.integers(0, 4, size=READLEN))
            else:
                s = int(stag[i])
                body = rev[g] if reverse else fwd[g]
                seq = "".join(bases[b] for b in rand[i, :s]) + body[:READLEN - s]
                if nerr[i]:
                    sl = list(seq)
                    for pos in rng.integers(0, READLEN, size=int(nerr[i])):
                        sl[pos] = bases[(bases.index(sl[pos]) + int(rng.integers(1, 4))) % 4]
                    seq = "".join(sl)
            fh.write(f"{hdr}:{tiles[i]}:{xs[i]}:{ys[i]} {mate}:N:0:1\n{seq}\n+\n{qual}\n")


INSTRUCTION = """\
You are the computational biologist of a functional-genomics core facility. A client ran a pooled
CRISPR-Cas9 knockout dropout screen with a secondary validation library and needs guide-level
counts and gene-level fitness results.

Task directory:
- `/task`

Inputs (read-only) under `/task/input`:
- `fastq/<sample_id>.fastq.gz` - one single-end FASTQ per sample ({n_samples} samples)
- `samples.tsv` - sample sheet (sample_id, timepoint, replicate, fastq_file, sequencing_run)
- `library.tsv` - sgRNA library with columns `sgRNA`, `sequence`, `gene`. Non-targeting control
  sgRNAs have `gene = NON_TARGETING`. Gene identities were blinded by the client.
- `screen_design.md` - wet-lab notes on the library, PCR and sequencing

Available software (on `PATH`): `cutadapt`, `mageck` (MAGeCK), Python 3.12 with numpy, scipy and
pandas. There is no internet access.

Your task:
1. Count how many reads carry each library sgRNA in every sample.
2. For every targeted gene, estimate the fitness effect of the knockout between T0 and T14 and
   decide whether the gene is a fitness gene.

Required outputs under `/task/output`:
- `counts.tsv` - tab-separated, header `sgRNA`, `gene`, then one column per `sample_id` in the
  order of `samples.tsv`. One row per library sgRNA (non-targeting sgRNAs included). Values are raw
  integer read counts.
- `gene_summary.tsv` - tab-separated with at least the columns `gene`, `lfc`, `is_hit`; one row per
  targeted gene (no row for `NON_TARGETING`). Extra columns are allowed.
  - `lfc`: log2 fold change of the gene, T14 versus T0, defined as the median over the gene's
    sgRNAs of the sgRNA-level log2 ratio of depth-normalised abundance (mean of the T14 replicates
    over mean of the T0 replicates). The scale must be such that a knockout without any fitness
    effect has `lfc` = 0.
  - `is_hit`: `1` if the gene's sgRNAs are significantly depleted compared with sgRNAs that have no
    fitness effect, at a false-discovery rate of 5%; otherwise `0`.
- `qc.json` - `{{"samples": {{"<sample_id>": {{"total_reads": <int>, "mapped_reads": <int>}}, ...}}}}`
  where `total_reads` is the number of reads in the FASTQ and `mapped_reads` the number of reads
  assigned to a library sgRNA.

Constraints:
- Do not modify anything under `/task/input`; write deliverables only under `/task/output`.
- Every sample in the sample sheet belongs to the experiment and must be reported.
- Do not ask for confirmation. Execute directly.
"""

JUDGE = """\
**Where the truth comes from.** The generator fixes, before any read exists, (a) which blinded
genes have a fitness effect and its size, (b) each sgRNA's efficiency, and (c) how many reads each
sgRNA contributes to each sample. `reference/truth.json` stores exactly these planted values. No
tool output is used as reference, so a mistake in the reference pipeline cannot become "truth".

**Hard gates (score 0).** Missing/unparseable file; `counts.tsv` not listing every library sgRNA
exactly once with non-negative integers; `gene_summary.tsv` not listing every targeted gene exactly
once; `is_hit` not in {0,1}; non-finite numbers.

**Components.**

| weight | component | metric | full credit | zero |
|---|---|---|---|---|
| 0.35 | guide_counts | per sample `sum(min(obs,true)) / max(sum obs, sum true)`, then the **minimum** over samples | >= 0.90 | <= 0.50 |
| 0.30 | gene_lfc | median absolute error of `lfc` against the planted latent median-guide LFC | <= 0.20 | >= 0.60 |
| 0.30 | hit_calls | Matthews correlation of `is_hit` with the planted fitness genes | >= 0.90 | <= 0.50 |
| 0.05 | qc_report | fraction of samples with exact `total_reads` and `mapped_reads` within [0.90, 1.005] x true | 1.0 | 0 |

**Why these tolerances.** Reads carry 0.2% substitution errors, so an exact-match counter recovers
about 96% of the true counts and a 1-mismatch counter about 100%; both must pass, hence 0.90. The
minimum (not the mean) over samples is used because a sample that was silently lost is a wrong
deliverable, not a 20% smaller one. Sampling noise on the gene LFC at this depth is about 0.05-0.10
log2 units and pseudocount choices move strongly depleted genes by < 0.1, while a wrong
normalisation shifts every gene by more than 1 log2 unit, so the 0.20/0.60 band separates the two
cleanly. MCC is used instead of F1 because in a hit-enriched library "call every gene a hit"
reaches F1 > 0.9 but MCC = 0. BH at 5% with ~40 neutral genes is expected to produce 0-2 false
hits, which still gives MCC >= 0.95.

**What is deliberately not scored.** P-values/FDR columns (method dependent) and any plot.
"""


def generate(seed: int, level: int, workdir: Path) -> TaskBuild:
    if level >= 4:   # level 4 has its own generator (own random stream; levels 1-3 stay byte-identical)
        from . import l4
        return l4.generate(seed, level, workdir)
    rng = np.random.default_rng(1000 + seed)
    hit_rich = level >= 3
    stagger = level >= 2
    n_fit, n_neu = (190, 40) if hit_rich else (60, 170)
    n_ntc, per_gene = 80, 4
    n_genes = n_fit + n_neu

    gene_ids = [f"TGT{i:04d}" for i in rng.permutation(np.arange(1, 9000))[:n_genes]]
    gene_ids.sort()
    fit_mask = np.zeros(n_genes, dtype=bool)
    fit_mask[rng.permutation(n_genes)[:n_fit]] = True
    gene_lfc = np.where(fit_mask, -rng.uniform(1.5, 4.0, size=n_genes), 0.0)

    seqs = _guides(rng, n_genes * per_gene + n_ntc)
    order = rng.permutation(len(seqs))
    guides = []   # (id, seq, gene, expected_lfc)
    k = 0
    for gi, g in enumerate(gene_ids):
        inactive = int(rng.integers(0, per_gene)) if rng.random() < 0.45 else -1   # at most one dud
        for j in range(per_gene):
            eff = rng.uniform(0.0, 0.08) if j == inactive else rng.uniform(0.75, 1.0)
            guides.append((f"{g}_sg{j + 1}", seqs[order[k]], g, float(eff * gene_lfc[gi])))
            k += 1
    for j in range(n_ntc):
        guides.append((f"NTC_{j + 1:03d}", seqs[order[k]], "NON_TARGETING", 0.0))
        k += 1
    lib_order = rng.permutation(len(guides))
    guides = [guides[i] for i in lib_order]
    G = len(guides)
    exp_lfc = np.array([g[3] for g in guides])

    samples = [("T0_rep1", "T0", 1), ("T0_rep2", "T0", 2),
               ("T14_rep1", "T14", 1), ("T14_rep2", "T14", 2), ("T14_rep3", "T14", 3)]
    rc_sample = samples[int(rng.integers(0, len(samples)))][0] if level >= 2 else None

    base = rng.lognormal(0.0, 0.45, size=G)
    fq_dir = workdir / "fastq"
    fq_dir.mkdir(parents=True, exist_ok=True)
    fwd = [(PREFIX + g[1] + SCAFFOLD) for g in guides]
    rev = [_rc(PREFIX + g[1] + SCAFFOLD[:40]) for g in guides]
    true_counts = np.zeros((G, len(samples)), dtype=int)
    total_reads, guide_reads = [], []
    for j, (sid, tp, rep) in enumerate(samples):
        ab = base * rng.lognormal(0.0, 0.12, size=G)
        if tp == "T14":
            ab = ab * 2.0 ** (exp_lfc + rng.normal(0.0, 0.10, size=G))
        depth = int(G * rng.uniform(150, 190))
        c = rng.multinomial(depth, ab / ab.sum())
        true_counts[:, j] = c
        n_junk = int(depth * rng.uniform(0.03, 0.06))
        origin = np.concatenate([np.repeat(np.arange(G), c), -np.ones(n_junk, dtype=int)])
        origin = origin[rng.permutation(len(origin))]
        reverse = sid == rc_sample
        _write_fastq(fq_dir / f"{sid}.fastq.gz", origin, fwd, rev, reverse, stagger, rng,
                     "B" if reverse else "A")
        total_reads.append(int(len(origin)))
        guide_reads.append(int(depth))

    # latent gene truth: median of the expected guide LFCs
    genes_truth = {}
    for g in gene_ids:
        vals = [x[3] for x in guides if x[2] == g]
        genes_truth[g] = {"lfc": float(np.median(vals)), "fitness": int(bool(fit_mask[gene_ids.index(g)]))}
    sids = [s[0] for s in samples]
    truth = {"samples": sids, "guides": {g[0]: g[2] for g in guides},
             "true_counts": {g[0]: [int(v) for v in true_counts[i]] for i, g in enumerate(guides)},
             "genes": genes_truth, "total_reads": total_reads, "guide_reads": guide_reads}

    library_tsv = "sgRNA\tsequence\tgene\n" + "".join(f"{g[0]}\t{g[1]}\t{g[2]}\n" for g in guides)
    samples_tsv = "sample_id\ttimepoint\treplicate\tfastq_file\tsequencing_run\n" + "".join(
        f"{sid}\t{tp}\t{rep}\t{sid}.fastq.gz\t{'B' if sid == rc_sample else 'A'}\n" for sid, tp, rep in samples)

    notes = ["# Screen design notes (wet lab)", "",
             "- Cells: Cas9-stable clone (identity blinded by the client), lentiviral infection at MOI ~0.3, "
             "500x library coverage maintained at every passage.",
             "- Library: VAL-2 secondary validation library in lentiGuide-Puro. "
             + (f"{n_fit} genes carried forward as candidate fitness genes from the primary genome-wide dropout "
                f"screen and {n_neu} genes picked as putative neutral" if hit_rich else
                f"{n_genes} genes from the client's pathway panel")
             + f"; {per_gene} sgRNAs per gene, plus {n_ntc} non-targeting sgRNAs.",
             "- Timepoints: T0 = day 3 after puromycin selection (2 replicates); T14 = 14 days later (3 replicates).",
             ("- gDNA PCR: one-step PCR with the P5 stagger primer mix (8 forward primers, 0-7 nt stagger) and one "
              "P7 index primer per sample." if stagger else
              "- gDNA PCR: one-step PCR with a single P5 primer (no stagger) and one P7 index primer per sample."),
             f"- Amplicon, sense strand 5'->3': ...{PREFIX} [20-nt spacer] {SCAFFOLD[:30]}... (sgRNA scaffold).",
             "- Sequencing run A: NextSeq 550, 75-cycle single read."]
    if rc_sample:
        notes.append(f"- {rc_sample} failed QC on run A. The core re-sequenced the same PCR product on run B "
                     "(NovaSeq, paired-end lane shared with other users); only the Read 2 file was delivered "
                     "and that is the file provided here.")
    screen_md = "\n".join(notes) + "\n"

    # ---- truth-derived probes / attacks ---------------------------------------------------
    def tsv_counts(mat):
        return "sgRNA\tgene\t" + "\t".join(sids) + "\n" + "".join(
            f"{g[0]}\t{g[2]}\t" + "\t".join(str(int(v)) for v in mat[i]) + "\n" for i, g in enumerate(guides))

    def tsv_genes(lfc_fn, hit_fn):
        return "gene\tlfc\tis_hit\n" + "".join(
            f"{g}\t{lfc_fn(g):.4f}\t{hit_fn(g)}\n" for g in gene_ids)

    def qc_json(mapped):
        return json.dumps({"samples": {s: {"total_reads": total_reads[j], "mapped_reads": int(mapped[j])}
                                       for j, s in enumerate(sids)}})

    prng = np.random.default_rng(seed + 77)
    thinned = prng.binomial(true_counts, 0.955)
    noisy_lfc = {g: genes_truth[g]["lfc"] + float(prng.normal(0, 0.08)) for g in gene_ids}
    neutral = [g for g in gene_ids if not genes_truth[g]["fitness"]]
    two_fp = set(neutral[:2])
    pipeline_py = {"crispr_pipeline.py": ASSETS / "crispr_pipeline.py"}

    probes = [
        Probe("planted_truth", "submission written directly from the planted truth", 0.999,
              files={"counts.tsv": tsv_counts(true_counts),
                     "gene_summary.tsv": tsv_genes(lambda g: genes_truth[g]["lfc"], lambda g: genes_truth[g]["fitness"]),
                     "qc.json": qc_json(guide_reads)}),
        Probe("within_tolerance", "4.5% of reads lost to sequencing errors, LFC noise sd 0.08, two false hits",
              0.99,
              files={"counts.tsv": tsv_counts(thinned),
                     "gene_summary.tsv": tsv_genes(lambda g: noisy_lfc[g],
                                                   lambda g: 1 if (genes_truth[g]["fitness"] or g in two_fp) else 0),
                     "qc.json": qc_json(thinned.sum(axis=0))}),
        Probe("independent_solution", "second implementation: no cutadapt/MAGeCK, bootstrap null from control guides",
              0.99, script="SOLUTION=2 bash /solution/solve.sh"),
    ]
    attacks = [
        Attack("all_genes_hit", "prior-only answer: every gene a hit with lfc -3, counts all 1", 0.05,
               files={"counts.tsv": tsv_counts(np.ones_like(true_counts)),
                      "gene_summary.tsv": tsv_genes(lambda g: -3.0, lambda g: 1),
                      "qc.json": qc_json([1] * len(sids))}),
        Attack("lfc_shifted_beyond_tolerance", "perfect counts and calls but every lfc shifted by +0.7", 0.75,
               files={"counts.tsv": tsv_counts(true_counts),
                      "gene_summary.tsv": tsv_genes(lambda g: genes_truth[g]["lfc"] + 0.7,
                                                    lambda g: genes_truth[g]["fitness"]),
                      "qc.json": qc_json(guide_reads)}),
    ]
    if level >= 2:
        attacks += [
            Attack("naive_fixed_offset", "spacer read at a fixed position, forward strand only", 0.60,
                   script="python3 /alegen_script/crispr_pipeline.py --extract fixed", script_files=pipeline_py),
            Attack("naive_drop_failed_sample", "correct trimming, but the Read-2 sample is not rescued: reported as zeros and left out of the statistics", 0.75,
                   script="python3 /alegen_script/crispr_pipeline.py --no-revcomp", script_files=pipeline_py),
        ]
    if level >= 3:
        attacks += [
            Attack("naive_library_defaults", "correct counts, but median normalisation and MAGeCK RRA defaults", 0.75,
                   script="python3 /alegen_script/crispr_pipeline.py --norm median --test rra", script_files=pipeline_py),
            Attack("naive_mageck_end_to_end", "mageck count on raw FASTQ with the stagger offsets, forward strand only, then mageck test defaults", 0.60,
                   script="python3 /alegen_script/crispr_pipeline.py --extract mageck_auto --norm median --test rra",
                   script_files=pipeline_py),
        ]

    traps = []
    if level >= 2:
        traps += [
            Trap("staggered_amplicon", "The spacer starts 24-31 nt into the read; a fixed-offset parser assigns ~1/8 of reads.",
                 "screen_design.md names the P5 stagger mix; per-sample mapping rate collapses with a fixed offset",
                 "Approach / wrong strategy"),
            Trap("read2_only_sample", f"{rc_sample} is the reverse complement of the amplicon; forward search maps ~0% of it.",
                 "screen_design.md says only Read 2 was delivered; read names end in ' 2:N:0'; mapping rate ~0 for one sample",
                 "Understanding / domain knowledge gap"),
        ]
    if level >= 3:
        traps.append(
            Trap("hit_enriched_library", "83% of genes are true fitness genes. Median/total normalisation shifts every "
                 "LFC by > 1 log2 unit and rank-based RRA misses most hits.",
                 "screen_design.md describes a validation library of carried-forward hits; non-targeting guides show a "
                 "large positive LFC after default normalisation",
                 "Understanding / domain knowledge gap"))

    return TaskBuild(
        task_id=f"crispr-screen-s{seed:04d}-l{level}", family="crispr_screen", seed=seed, level=level,
        title="Pooled CRISPR dropout screen: FASTQ to fitness genes",
        summary="Quantify sgRNA counts from staggered amplicon reads of a pooled CRISPR knockout screen and "
                "call fitness genes in a validation library with blinded gene identities.",
        subdomain="3.2 Genomics & Sequence Analysis",
        software=["cutadapt", "MAGeCK", "Python", "numpy/scipy/pandas"],
        instruction_md=INSTRUCTION.format(n_samples=len(samples)),
        agent_dockerfile=dockerlib.agent_dockerfile(conda=["mageck=0.5.9.5", "cutadapt=5.2"]),
        verifier_dockerfile=dockerlib.verifier_dockerfile(pip=[]),
        input_files={**{f"fastq/{s}.fastq.gz": fq_dir / f"{s}.fastq.gz" for s in sids},
                     "library.tsv": library_tsv, "samples.tsv": samples_tsv, "screen_design.md": screen_md},
        solution_files={"solve.sh": ASSETS / "solve.sh", "crispr_pipeline.py": ASSETS / "crispr_pipeline.py",
                        "solve_alt.py": ASSETS / "solve_alt.py"},
        tests_files={"score.py": ASSETS / "score.py", "reference/truth.json": json.dumps(truth)},
        traps=traps, attacks=attacks, probes=probes,
        capabilities=["read-structure reasoning from lab notes", "per-sample QC before analysis",
                      "choosing statistics that match the experimental design", "exact output contracts"],
        judge_explanation_md=JUDGE, pass_threshold=0.9,
        oracle_steps=9, oracle_runtime_min=2.0, cpus=4, memory_mb=8000, storage_mb=10000,
    )
