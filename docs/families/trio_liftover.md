# Family `trio_liftover` - rare-disease trio, cross-build variant prioritisation

ALE subdomain 3.2 (Genomics & Sequence Analysis). Code: `pipeline/alegen/families/trio_liftover/`.
Spec: `docs/specs/trio_liftover.md`. Not a copy of an ALE public task.

## Scenario

A clinical genetics lab receives a joint trio VCF (affected child + both parents, chr22 only) from an
external provider. The lab's transcript set and population-frequency table are on GRCh38. The analyst
must (1) deliver the call set on GRCh38 in canonical representation (contig `chr22`, one ALT per record,
left-aligned, minimal), (2) list the input alleles that have no GRCh38 position under the provided
chain, and (3) apply the lab SOP (usable genotypes, allele balance, rarity, qualifying consequence,
panel, three inheritance models) and report the candidate causal variants.

Agent sees: `trio.vcf.gz(.tbi)`, `pedigree.ped`, `panel_genes.txt`, `lab_sop.md`, and `resources/` with
N-masked GRCh37/GRCh38 chr22 FASTA (bgzip + fai + gzi), both chr22 liftOver chains, a GFF3 with one
canonical transcript per gene and a GRCh38 frequency table. Tools: bcftools 1.21 (+liftover plugin),
samtools/htslib 1.21, picard-slim 3.3.0, Python 3.12 with pysam 0.23.0, numpy, scipy, pandas.
Input size 2.5 MB, agent image 1.96 GB, generation 5-10 s per attempt (79 of 80 seed/level combinations
succeed on the first attempt), oracle 0.07-0.35 min in `make` on the shared build host.

## Reference material (cached under `.cache/trio_liftover/dl/`)

| file | URL | sha256 |
|---|---|---|
| hg19.chr22.fa.gz | https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr22.fa.gz | 4a5a0049b23232b163e575bf506b8af86e7d563fa9990311fb9fa3346d6865ec |
| hg38.chr22.fa.gz | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz | 05f9d97d6fbfd08a44ca45b50837ca2ae9c471f35ba79dffec04d2cb5eaaf695 |
| hg19ToHg38.over.chain.gz | https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz | 5c0598e500ceb5a78c73086929e8ef993aec309bcafb595139b53d440b125a1d |
| hg38ToHg19.over.chain.gz | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz | 14a712e8e147d9fc8e9d87d51977b46f6f8ddb93efbe5d0843d86b6205f587b1 |
| Homo_sapiens.GRCh38.113.chromosome.22.gff3.gz | https://ftp.ensembl.org/pub/release-113/gff3/homo_sapiens/Homo_sapiens.GRCh38.113.chromosome.22.gff3.gz | 20b5eff9aa90c936d5cc9c6fa514c6a9678cfebe0aa3b3882b84337284ecc102 |

Derived, seed-independent resources are cached in `.cache/trio_liftover/derived_v1/resources/` (built once
inside the tools image). Five 600-kb GRCh38 windows keep real sequence (18.90-19.50, 21.40-22.00,
23.60-24.20, 41.55-42.15, 50.05-50.65 Mb); the matching hg19 windows are the hg19 spans of the main
chain blocks that map into them; everything else is N. The windows hold 88 protein-coding genes with an
Ensembl-canonical mRNA (87 with a clean CDS, which are the ones used for planting) and different
hg19->GRCh38 offsets per window (+12477, about -354289, -342187/-396448, -396004, -438429), so no single
shift replaces the lift-over. Chains are restricted to chr22->chr22 (11 chains forward). With that restriction
the 54-kb hg19 segment chr22:24,343,276-24,397,527 (GSTT region, on an ALT contig in the full chain) has
no GRCh38 image.

## Truth model

Everything is planted on GRCh38 first; the provider VCF is derived from it.

1. Every allele is created as a left-aligned minimal GRCh38 allele (own normaliser; `bcftools norm -c e`
   on the planted site list must return it unchanged, otherwise generation aborts).
2. An allele is only planted where (a) the position lies >= 100 bp inside an ungapped plus-strand block
   of the main chain, (b) the hg19 and GRCh38 sequences are identical for 50 bp to the left and
   50 bp + allele span + maximal right-shift to the right, and (c) nothing else is planted within 12 bp
   (deliberate shared-record pairs excepted). Under these conditions the hg19 allele is the same
   allele at `pos - offset`, and any correct lift-over returns the planted key.
3. Coding variants are designed from the CDS of the single shipped transcript: codon substitutions
   (missense / stop_gained / synonymous, >= 6 bp from CDS-segment ends, not in the first/last 8 codons),
   1-2 bp frameshift indels (>= 8 bp inside the segment over the whole shift range, no stop within the
   next 3 codons), canonical GT/AG splice-site SNVs, 3-bp in-frame deletions. Only genes whose position
   is covered by exactly one gene span are used. Non-coding background stays >= 30 bp away from any exon.
4. Trio genotypes are planted (HWE parents + Mendelian transmission for common background; explicit
   patterns for rare variants); `GT:AD:DP:GQ` is simulated from the planted genotype.
5. Causal set C = 1 de novo + 1 compound-het pair (frameshift from one parent, missense from the other)
   + (with probability 0.6) 1 homozygous-recessive, each in a different panel gene. One causal SNV
   (seed-dependent: the de novo or the compound-het missense) shares its input record with a common
   in-frame deletion and is the 2nd ALT; the compound-het frameshift is delivered right-shifted.
6. Everything else fails at least one SOP rule by construction: non-coding; common; synonymous;
   non-panel gene; a single inherited heterozygous hit per spare panel gene. The generator runs the SOP
   as code over the planted data and requires the result to be exactly C.
7. hg19-only alleles: 26-33 SNVs at hg19 positions that no shipped chain block covers (all but 8 in the
   54-kb region, 8 in 1-60 bp chain gaps; 21-28 of them end up scored).
8. Input derivation: back-map to hg19, contig `22` (level >= 2), 65% of indels right-shifted or padded,
   2% of background SNVs padded, ~125 two-SNV sites and ~110 deletion+SNV pairs merged into
   multi-allelic records (never both ALTs in one sample), no assembly tag in the header (level >= 2).

Guards with domain tools (tools image, generation time only):
- `bcftools csq -l` on the planted sites: the consequence of every causal variant and of every decoy
  that is meant to be damaging must be exactly `[(gene, designed class)]`; synonymous / wrong-build
  decoys must have no qualifying consequence; the SOP run on the csq consequences of *all* planted
  alleles must return exactly C. Any failure -> the whole instance is redrawn with the next sub-seed
  (max 6).
- Lift-over consensus: the derived input is renamed, split/left-aligned on hg19 and lifted with Picard
  LiftoverVcf 3.3.0 and with bcftools +liftover (2024-09-27), each followed by `bcftools norm`. A planted
  allele is scored only if both return the planted key *and* the planted genotypes; an hg19-only allele
  is scored only if both reject it. Everything else is gray (not rewarded, not penalised, including any
  key either tool produced for it). A causal or decoy variant in the gray zone triggers a redraw.
  Measured gray zone: 4 / 5 / 5 alleles (s1-l3 / s2-l3 / s1-l1) and 4-7 in 16 further generated
  instances, all of them hg19-only alleles in tiny chain gaps that bcftools +liftover places on GRCh38
  (reference allele added / alleles swapped) while Picard rejects them. No mapped allele was ever gray,
  and no instance needed a redraw because of a tool disagreement.

Why the reference is right even if a tool is wrong: keys, genotypes and the candidate list come from
the planting code; the tools can only *remove* items from scoring or force a redraw.

## Traps and their discoverability

| trap (level) | what goes wrong | how an analyst finds it |
|---|---|---|
| build_mismatch (2+) | input is hg19, nothing says so; annotated as GRCh38 everything is misplaced, and one planted record (hg19 coordinate = a missense position of a panel gene on GRCh38, REF identical in both assemblies) becomes a "de novo missense" | `##contig` length 51304566 = hg19 chr22 (GRCh38: 50818468); REF concordance 100% vs GRCh37 FASTA, 10.3-10.6% vs GRCh38 (measured on s1-l3, s2-l3); `bcftools csq` aborts with "the fasta reference does not match the VCF REF allele"; SOP states lab resources are GRCh38 |
| contig_naming (2+) | VCF says `22`, all resources say `chr22` | Picard rejects 100% as `NoTarget`; `bcftools norm` stops with `sequence "22" not found` |
| representation (2+) | 25% of alleles are right-shifted, padded or in multi-allelic records; a causal SNV is the 2nd ALT with a padded representation; a common frameshift (AF 1.2-6%, hom-alt in the child, both parents het) is right-shifted so a literal frequency lookup misses it | output contract and SOP rule 1 state the canonical form; `bcftools norm` reports split/realigned counts; multi-allelic ALT columns are visible |
| genotype_quality_decoys (3) | 2 "de novo" with an unusable parental genotype (DP 3-6 or GQ 5-15), 2 with proband allele balance 9-19% | SOP rules 2 and 3 + FORMAT fields |
| same_parent_pairs (3) | 2 panel genes with two rare damaging hets inherited from the same parent | SOP rule 7 + parental genotypes |

Not declared as a trap: hg19-only alleles (`unmapped.tsv`). Both lift-over tools hand the rejects over
and the instruction asks for them explicitly, so they are an ordinary deliverable. Also planted at all
levels (plain SOP filtering, not traps): a synonymous de novo in a panel gene, a damaging de novo and
up to three recessive-looking damaging variants in non-panel genes.

## Output contract

- `trio.GRCh38.vcf.gz` + `.tbi`/`.csi`: GRCh38, contig `chr22`, one ALT per record, left-aligned minimal,
  REF = GRCh38, the three original sample names with GT.
- `unmapped.tsv`: `chrom_in, pos_in, ref, alt` as written in the input.
- `candidates.tsv`: `chrom, pos, ref, alt, gene, consequence, model`.

## Scoring (pass threshold 0.90)

Gates (score 0): files missing/unparseable, VCF samples != the three input samples, no index file,
VCF > 10 MB compressed or > 50 MB uncompressed (streaming capped decompression; a real submission is
< 1 MB), > 200 candidate rows, `model` outside the vocabulary, non-nucleotide candidate alleles.

| weight | component | metric | full | zero |
|---|---|---|---|---|
| 0.15 | lifted_callset | F1 of literal (chrom,pos,ref,alt) vs scored planted alleles; gray keys ignored | >= 0.985 | <= 0.90 |
| 0.10 | genotypes | min over samples of fraction of scored alleles present with the planted GT (phasing/order ignored) | >= 0.985 | <= 0.90 |
| 0.10 | unmapped | F1 of (pos_in,ref,alt), `22`/`chr22` both accepted | >= 0.95 | <= 0.50 |
| 0.45 | candidates | F1 of candidate keys vs planted causal set | 1.0 | <= 0.70 |
| 0.20 | candidate_annotation | fraction of causal variants with correct gene + consequence class + model (`*_variant` suffix and `&`-joined terms accepted) | 1.0 | 0 |

Weights differ from the spec sketch (0.20/0.10/0.10/0.40/0.20): the three "earlier artefact" components
now carry 0.35 as AUTHORING.md asks.

### Measured values

Each cell: lifted F1 / min GT concordance / unmapped F1 / candidate TP-FP-FN (F1) / annotated -> **total score**.
Values are read from `details.json` of the final `make` runs (`results/work/<task>/`).

| submission | s1-l3 (3 causal) | s2-l3 (3 causal) | s3-l3 (4 causal) | s1-l1 (4 causal) |
|---|---|---|---|---|
| oracle (Picard + bcftools csq) | 1.000 / 1.000 / 1.000 / 3-0-0 (1.000) / 3 -> **1.000** | 1.000 / 1.000 / 1.000 / 3-0-0 (1.000) / 3 -> **1.000** | 1.000 / 1.000 / 1.000 / 4-0-0 (1.000) / 4 -> **1.000** | 1.000 / 1.000 / 1.000 / 4-0-0 (1.000) / 4 -> **1.000** |
| independent solution (pure Python) | 1.000 / 1.000 / 1.000 / 3-0-0 (1.000) / 3 -> **1.000** | 1.000 / 1.000 / 1.000 / 3-0-0 (1.000) / 3 -> **1.000** | 1.000 / 1.000 / 1.000 / 4-0-0 (1.000) / 4 -> **1.000** | 1.000 / 1.000 / 1.000 / 4-0-0 (1.000) / 4 -> **1.000** |
| within_tolerance probe | 0.997 / 0.994 / 0.980 / 3-0-0 (1.000) / 3 -> **1.000** | 0.997 / 0.994 / 0.981 / 3-0-0 (1.000) / 3 -> **1.000** | 0.997 / 0.994 / 0.981 / 4-0-0 (1.000) / 4 -> **1.000** | 0.997 / 0.993 / 0.976 / 4-0-0 (1.000) / 4 -> **1.000** |
| naive: no lift-over | 0.000 / 0.000 / 0.000 / 0-1-3 (0.000) / 0 -> **0.000** | 0.000 / 0.000 / 0.000 / 0-1-3 (0.000) / 0 -> **0.000** | 0.000 / 0.000 / 0.000 / 0-1-4 (0.000) / 0 -> **0.000** | n/a (not in battery) |
| naive: contig name as delivered | 0.000 / 0.000 / 0.017 / 0-0-3 (0.000) / 0 -> **0.000** | 0.000 / 0.000 / 0.018 / 0-0-3 (0.000) / 0 -> **0.000** | 0.000 / 0.000 / 0.018 / 0-0-4 (0.000) / 0 -> **0.000** | n/a (not in battery) |
| naive: no split / left-align | 0.774 / 0.744 / 1.000 / 1-1-2 (0.400) / 1 -> **0.167** | 0.780 / 0.749 / 1.000 / 1-1-2 (0.400) / 1 -> **0.167** | 0.771 / 0.740 / 1.000 / 2-2-2 (0.500) / 2 -> **0.200** | n/a (not in battery) |
| naive: quality rules ignored | 1.000 / 1.000 / 1.000 / 3-4-0 (0.600) / 3 -> **0.550** | 1.000 / 1.000 / 1.000 / 3-4-0 (0.600) / 3 -> **0.550** | 1.000 / 1.000 / 1.000 / 4-4-0 (0.667) / 4 -> **0.550** | n/a (not in battery) |
| naive: phase-blind compound het | 1.000 / 1.000 / 1.000 / 3-4-0 (0.600) / 3 -> **0.550** | 1.000 / 1.000 / 1.000 / 3-4-0 (0.600) / 3 -> **0.550** | 1.000 / 1.000 / 1.000 / 4-4-0 (0.667) / 4 -> **0.550** | n/a (not in battery) |
| naive: all of the above | 0.000 / 0.000 / 0.000 / 0-1-3 (0.000) / 0 -> **0.000** | 0.000 / 0.000 / 0.000 / 0-1-3 (0.000) / 0 -> **0.000** | 0.000 / 0.000 / 0.000 / 0-1-4 (0.000) / 0 -> **0.000** | n/a (not in battery) |
| truth + 1 false candidate | 1.000 / 1.000 / 1.000 / 3-1-0 (0.857) / 3 -> **0.786** | 1.000 / 1.000 / 1.000 / 3-1-0 (0.857) / 3 -> **0.786** | 1.000 / 1.000 / 1.000 / 4-1-0 (0.889) / 4 -> **0.833** | 1.000 / 1.000 / 1.000 / 4-1-0 (0.889) / 4 -> **0.833** |
| truth - 1 candidate | 1.000 / 1.000 / 1.000 / 2-0-1 (0.800) / 2 -> **0.633** | 1.000 / 1.000 / 1.000 / 2-0-1 (0.800) / 2 -> **0.633** | 1.000 / 1.000 / 1.000 / 3-0-1 (0.857) / 3 -> **0.736** | 1.000 / 1.000 / 1.000 / 3-0-1 (0.857) / 3 -> **0.736** |
| perfect call set, 190 bogus candidates | 1.000 / 1.000 / 1.000 / 0-190-3 (0.000) / 0 -> **0.350** | 1.000 / 1.000 / 1.000 / 0-190-3 (0.000) / 0 -> **0.350** | 1.000 / 1.000 / 1.000 / 1-189-3 (0.010) / 0 -> **0.350** | 1.000 / 1.000 / 1.000 / 0-190-4 (0.000) / 0 -> **0.350** |
| constant empty answer | 0.000 / 0.000 / 0.000 / 0-0-3 (0.000) / 0 -> **0.000** | 0.000 / 0.000 / 0.000 / 0-0-3 (0.000) / 0 -> **0.000** | 0.000 / 0.000 / 0.000 / 0-0-4 (0.000) / 0 -> **0.000** | 0.000 / 0.000 / 0.000 / 0-0-4 (0.000) / 0 -> **0.000** |

Generic attacks (empty output, inputs copied, zero-byte files, header only) score 0.000 on all four tasks.
Harbor (s1-l3): `oracle` reward 1.0, `nop` reward 0.0.

Reading of the bands:
- lifted_callset / genotypes: both correct solutions 1.000 on both seeds; three further correct routes
  measured outside the battery on s1-l3 (bcftools +liftover on raw records; +liftover after
  pre-normalisation; Picard on un-split records; each followed by `bcftools norm -m-any` on GRCh38) also
  give F1 1.000 / GT 1.000 / unmapped F1 1.000. The no-normalise variant gives F1 0.774-0.780 and GT
  0.744-0.749. The band [0.90, 0.985] sits between with 0.015 margin above and 0.12 below; the
  within-tolerance probe (0.7% of records dropped) sits at 0.997 / 0.993.
- candidates: with 3-4 causal variants one extra row gives F1 0.857-0.889 and one missing row 0.80-0.857,
  i.e. component value 0.33-0.63; measured totals 0.633-0.833 < 0.90, so a wrong short list never
  passes. Every single-trap naive variant ends at F1 <= 0.667 and gets nothing for the component;
  the highest naive total is 0.55 (ceiling 0.75 = pass threshold - 0.15).
- unmapped: the tools hand this over; the band only has to separate "reported" from "not reported".

## Deviations from the spec

1. **CrossMap replaced by bcftools +liftover.** CrossMap 0.7.3 (the current bioconda build for
   Python 3.12 that I tested) lifts only the first base of REF and writes deletions as `REF=A ALT=`
   (empty) and padded substitutions as bogus insertions (verified on a hand-made VCF and in
   `cmmodule/mapvcf.py`). It cannot serve as a consensus tool. The agent image therefore ships Picard and
   the bcftools liftover plugin instead of CrossMap; older CrossMap versions were not tested.
2. `liftover_rejects` is a deliverable, not a declared trap (see above); hence no naive variant for it.
3. Weights 0.15/0.10/0.10/0.45/0.20 instead of 0.20/0.10/0.10/0.40/0.20.
4. Windows are fixed across seeds (panel, causal genes, every variant and all decoys vary).
5. Decoys were multiplied (2 same-parent genes, 4 quality decoys) so that every single-trap naive
   variant loses the whole candidate component, not just part of it.
6. The "common, non-normalised" decoy is a right-shifted common frameshift in a recessive pattern; the
   "two hets from one parent" decoy exists twice; "rare damaging variant in a non-panel gene" exists as a
   de novo and as recessive-looking variants.
7. GFF3 is shipped with `chr22` as seqid (Ensembl ships `22`) so that the contig-name trap is confined
   to the provider VCF.

## Known limitations / things I am not sure about

- **Scored set is conditioned on two tools.** By construction this only removed hg19-only alleles in
  1-60 bp chain gaps (5-7 per task). Whether such an allele "has a GRCh38 position" is a genuine
  definitional question, which is why it is gray rather than decided by me.
- **All mapped alleles are easy lift-over cases by construction** (plus strand, identical flanks, far
  from block edges). Minus-strand blocks, REF/ALT swaps between assemblies and block-edge indels - the
  really nasty part of real lift-overs - are deliberately absent, because the two tools disagree there
  and no planted truth exists. The task tests build detection, naming, normalisation and reject
  handling, not lift-over corner cases.
- **Candidate F1 on 3-4 items is coarse**: the component is effectively all-or-nothing plus one
  intermediate step. That is intended (a short list with a wrong entry fails) but gives little gradient.
  The annotation component (0.20) and the call-set components (0.35) are the partial-credit signal.
- **`naive_no_liftover` relies on a "drop REF-inconsistent records and carry on" reaction**, because
  bcftools csq refuses to annotate hg19 coordinates on GRCh38. A different careless reaction would give
  a different but equally wrong output; all of them score 0 on the call set.
- **`naive_no_rename` is degenerate** (0 records lifted, everything in unmapped.tsv); it is the literal
  result of running Picard on the file as delivered, but a real analyst would notice.
- **Consequence vocabulary**: only five classes qualify and planted coding variants avoid start/stop
  codons, splice regions and exon-spanning indels, so the class is unambiguous for the scored variants.
  Background variants are not guaranteed to get the same class from every annotator; the guard only
  guarantees that with bcftools csq *and* with the designed classes none of them becomes a candidate.
  An annotator with different transcript handling could in principle promote a background variant; the
  independent solution (own annotator) did not on either seed.
- **Haplotype-aware csq**: without `-l`, `bcftools csq -p a` treats the unphased `0/1` calls of the
  compound-het pair as one haplotype and may report the downstream variant as `@pos`/`*missense`. The SOP
  now says consequences are per variant; an agent still has to notice.
- Sample columns: Picard reorders samples alphabetically; the verifier matches by name, so both orders
  score the same.
- Level 2 (seed 1) passed every gate with `--no-harbor` in the round before the verifier's decompression
  was hardened (log kept as `results/trio_s1_l2.pre_hardening.log`); it was not repeated afterwards because
  the build host ran short of disk, and its task directory was removed so that no stale task is shipped.
- The index gate only checks that a non-empty `.tbi`/`.csi` file exists; it is not validated.
- Both required level-3 seeds happened to draw the same configuration (3 causal variants, compound-het
  missense as 2nd ALT). Seed 3 (4 causal variants, de novo as 2nd ALT) was therefore validated in
  addition; see the table.
- Real gene symbols are kept (not blinded): the panel and the causal genes are drawn at random per seed,
  so prior knowledge about chr22 disease genes does not help.
- The frequency table contains the normalised GRCh38 key of ~75% of the planted alleles mixed with
  2500-3200 unrelated rows. It cannot replace the lift-over (no genotypes, ~50% decoy rows), but it is
  a weak hint about canonical representation.

---

# Level 4 (2026-11-12)

Spec: `docs/specs/l4/trio_liftover.md`. Code: the `lvl4` branches of `pipeline/alegen/families/trio_liftover/__init__.py`
and five new assets (`score_l4.py`, `trio_pipeline_l4.py`, `solve_alt_l4.py`, `gen_check_l4.sh`, `solve_l4.sh`).
Levels 1-3 are untouched: level 3 seed 1 regenerates byte-identical to `tasks/trio-liftover-s0001-l3` (`diff -r`,
nothing ignored), and level 3 seeds 2-3 and level 1 seed 1 regenerate the `truth.json`, `trio.vcf` and `stats.json`
of the earlier runs kept under `results/gen/`.

## Why (pilot finding -> element)

The pilot model scored 1.0 on level 3 in 3 min 23 s with 18 tool calls. The level-3 limitations above name the
reasons: every mapped allele is an easy lift-over case ("run the tool" is the whole call-set stage), consequences
are per variant, and three or four candidates give no gradient. Level 4 keeps the level-3 task with all its five
traps and adds the lift-over cases that were left out, with exact truth: the family's alleles are planted on GRCh38
and the provider's file is derived from them, so no tool has to be trusted for the new cases.

| element | pilot finding it answers | what is planted | why the truth is exact |
|---|---|---|---|
| **P_reference_differences** | (b) a plausible shortcut - the lift-over tool's output taken as the deliverable - that the generator closes | 178 aligned positions where the provider's reference base differs from GRCh38: the 173 that exist between the real hg19 and GRCh38 sequences inside the shipped windows plus 5 planted into the shipped GRCh37 slice inside coding sequence. On seed 1: 99 have **no provider record** (all three members homozygous for the provider's base = homozygous ALT on GRCh38; 95-107 on seeds 2-8), 43 carry a record whose ALT is the GRCh38 base (complementary genotypes on GRCh38), 21 carry a third allele, 12 are "everybody carries the GRCh38 base" (the record disappears), 3 are the planted candidate sites. On such sites sit a homozygous-recessive variant (provider shows the child as 0/0), one partner of a compound het (provider shows the "other" parent as 1/1) and a de novo third allele (provider shows parents 1/1, child 1/2). | each member's allele pair is written once against the GRCh38 base (truth) and once against the provider's base (input); no lift-over tool is involved |
| **P_minus_strand_block** | same | a segment around an exon of a panel gene (3.0-6.0 kb of flanks plus the exon by design, 4.6-6.5 kb on seeds 1-8) is replaced by its reverse complement in the shipped GRCh37 slice; both chain files carry it as an extra minus-strand chain. 38-49 scored alleles inside (13-21 indels) on seeds 1-8, including a compound-het pair (frameshift + missense). | the provider allele is the reverse complement of the planted allele, left-aligned on the provider's strand; scored only if Picard and bcftools +liftover both bring it back to the planted key and genotypes (they did for every allele, see below) |
| **P_haplotype_consequence** | (b) per-variant logic where the stated rule is per haplotype | SOP rule 5 assesses single-nucleotide variants that share a codon on one chromosome copy of the proband together. Planted: a *trans* pair in one codon (each missense, together silent: a true compound het, both reported), a *cis* pair (each qualifying, together silent) plus a damaging variant from the other parent (gene must not be reported; per-variant logic reports three rows), a de novo that is silent alone but a stop/missense together with a homozygous position without provider record (reported), and the reverse (missense alone, silent on the real haplotype: not reported). | the generator runs the SOP including rule 5 as code on the planted data and requires exactly the planted set; `bcftools csq -l` confirms every single-variant class |
| **C_dependent_stages** | (a)/(c) an early mistake that only shows in the final deliverable | the two de novo cases above depend on a position that no input record mentions: missing it in the lift-over stage changes the candidate list (one lost, one false) | as above |
| **Q_gradient** | three items gave no gradient | 13 scored candidates in 9 genes (level 3: 3-4 in 2-3). This is a scoring change, **not** a Q element in the sense of AUTHORING.md (no in-data control). | - |

The instruction defines the call set by meaning (how the three family members differ from GRCh38 at every aligned
position, given the provider's statement that the file is complete against the provider's reference), defines `GT`
as the number of copies of the record's ALT allele, and defines `unmapped.tsv` by "position not aligned by the
chain". `lab_sop.md` (rev. 5) restates every rule in terms of ALT copies against GRCh38, says that genotypes that
follow from the completeness statement are usable, and states rule 5 completely (which variants are assessed
together, how a copy is assigned from the trio genotypes, what happens when it cannot be assigned, that every
carried variant counts as a partner, and what a homozygous variant needs). No method, tool option or threshold of a
method is named.

## Reference solutions (level 4)

- Oracle `trio_pipeline_l4.py` (24 external tool invocations, listed in its docstring): rename -> `bcftools +liftover`
  on the records **as delivered** -> sort -> `bcftools norm -m-any` on GRCh38 -> drop records nobody carries ->
  walk the chain blocks (both strands), compare the assemblies and add the positions without a record as `1/1` ->
  `bcftools csq -l` for the single-variant classes -> assign each proband variant to a copy -> one `bcftools csq`
  run per copy on a VCF that holds the SNVs of that copy as homozygous -> SOP in Python.
- Independent solution `solve_alt_l4.py`: no Picard, no bcftools. Own chain parser with strand, own re-expression
  at differing bases (records and positions without a record), own left-alignment, own codon logic per copy.
- A third implementation of rule 5 lives in the generator (`haplotype_cons`) and is what the truth is checked
  against.

## Scoring (pass 0.90)

| weight | component | metric | full | zero |
|---|---|---|---|---|
| 0.08 | callset_standard | F1 over scored alleles at plus-strand positions where both assemblies agree; a TP needs the literal key and the planted ALT copies of all three samples | >= 0.985 | <= 0.90 |
| 0.12 | callset_reference_differences | the same over alleles at positions where the provider's base differs from GRCh38 | >= 0.97 | <= 0.80 |
| 0.07 | callset_inverted_block | the same over alleles inside the minus-strand block | >= 0.95 | <= 0.70 |
| 0.08 | unmapped | F1 of (pos_in, ref, alt) | >= 0.95 | <= 0.50 |
| 0.45 | candidates | F1 of candidate keys vs the 13 planted | 1.0 | <= 0.87 |
| 0.20 | candidate_annotation | fraction of planted candidates with right gene + class + model | 1.0 | 0 |

Call set + unmapped = 0.35, candidate list = 0.65. Records in which no sample carries the ALT allele are ignored.
Differences from the level-3 verifier: genotypes are folded into the call-set F1 (a key with wrong copies is a
wrong record), false positives are attributed to the category of their position, and the candidate ramp starts at
0.87 instead of 0.70 because with 13 candidates one wrong row gives F1 0.96 instead of 0.86.

## Measured values

`make --keep-work`, values read from `results/work/<task>/**/details.json`. Each cell: F1 standard / F1 reference
differences / F1 inverted block / F1 unmapped / candidate TP-FP-FN (F1) / annotated -> **total**.

| submission | s1-l4 (`make`, with Harbor) | s2-l4 (`make --no-harbor`) |
|---|---|---|
| oracle | 1.000 / 1.000 / 1.000 / 1.000 / 13-0-0 (1.000) / 13 -> **1.000** | 1.000 / 1.000 / 1.000 / 1.000 / 13-0-0 (1.000) / 13 -> **1.000** |
| probe:independent_solution | 1.000 / 1.000 / 1.000 / 1.000 / 13-0-0 (1.000) / 13 -> **1.000** | 1.000 / 1.000 / 1.000 / 1.000 / 13-0-0 (1.000) / 13 -> **1.000** |
| probe:planted_truth | 1.000 / 1.000 / 1.000 / 1.000 / 13-0-0 (1.000) / 13 -> **1.000** | 1.000 / 1.000 / 1.000 / 1.000 / 13-0-0 (1.000) / 13 -> **1.000** |
| probe:within_tolerance | 0.997 / 0.997 / 1.000 / 0.980 / 13-0-0 (1.000) / 13 -> **1.000** | 0.997 / 1.000 / 1.000 / 0.981 / 13-0-0 (1.000) / 13 -> **1.000** |
| attack:added_sites_missing_from_truth | 1.000 / 0.639 / 1.000 / 1.000 / 13-0-0 (1.000) / 13 -> **0.880** | 1.000 / 0.652 / 1.000 / 1.000 / 13-0-0 (1.000) / 13 -> **0.880** |
| attack:candidates_at_provider_coordinates | 1.000 / 1.000 / 1.000 / 1.000 / 0-13-13 (0.000) / 0 -> **0.350** | 1.000 / 1.000 / 1.000 / 1.000 / 0-13-13 (0.000) / 0 -> **0.350** |
| attack:codon_merge_ignoring_phase | 1.000 / 1.000 / 1.000 / 1.000 / 11-0-2 (0.917) / 11 -> **0.681** | 1.000 / 1.000 / 1.000 / 1.000 / 11-0-2 (0.917) / 11 -> **0.681** |
| attack:constant_answer | 0.000 / 0.000 / 0.000 / 0.000 / 0-0-13 (0.000) / 0 -> **0.000** | 0.000 / 0.000 / 0.000 / 0.000 / 0-0-13 (0.000) / 0 -> **0.000** |
| attack:empty_output | gate: required output file missing: trio.GRCh3 -> **0.000** | gate: required output file missing: trio.GRCh3 -> **0.000** |
| attack:every_panel_variant | 1.000 / 1.000 / 1.000 / 1.000 / 0-190-13 (0.000) / 0 -> **0.350** | 1.000 / 1.000 / 1.000 / 1.000 / 0-190-13 (0.000) / 0 -> **0.350** |
| attack:header_only | gate: trio.GRCh38.vcf.gz is empty -> **0.000** | gate: trio.GRCh38.vcf.gz is empty -> **0.000** |
| attack:ignore_reference_differences | 1.000 / 0.000 / 1.000 / 0.349 / 8-1-5 (0.727) / 8 -> **0.273** | 1.000 / 0.000 / 1.000 / 0.348 / 8-1-5 (0.727) / 8 -> **0.273** |
| attack:inputs_copied | gate: required output file missing: trio.GRCh3 -> **0.000** | gate: required output file missing: trio.GRCh3 -> **0.000** |
| attack:inverted_block_lost | 1.000 / 1.000 / 0.000 / 1.000 / 11-0-2 (0.917) / 11 -> **0.611** | 1.000 / 1.000 / 0.000 / 1.000 / 11-0-2 (0.917) / 11 -> **0.611** |
| attack:inverted_block_unmapped | 1.000 / 1.000 / 0.000 / 0.531 / 11-0-2 (0.917) / 11 -> **0.536** | 1.000 / 1.000 / 0.000 / 0.551 / 11-0-2 (0.917) / 11 -> **0.540** |
| attack:l3_reference | 1.000 / 0.000 / 1.000 / 0.349 / 8-4-5 (0.640) / 8 -> **0.273** | 1.000 / 0.000 / 1.000 / 0.348 / 8-4-5 (0.640) / 8 -> **0.273** |
| attack:naive_default_pipeline | 0.000 / 0.000 / 0.000 / 0.000 / 0-1-13 (0.000) / 0 -> **0.000** | 0.000 / 0.000 / 0.000 / 0.000 / 0-1-13 (0.000) / 0 -> **0.000** |
| attack:naive_ignore_quality | 1.000 / 1.000 / 1.000 / 1.000 / 13-4-0 (0.867) / 13 -> **0.550** | 1.000 / 1.000 / 1.000 / 1.000 / 13-4-0 (0.867) / 13 -> **0.550** |
| attack:naive_no_liftover | 0.000 / 0.000 / 0.000 / 0.000 / 0-1-13 (0.000) / 0 -> **0.000** | 0.000 / 0.000 / 0.000 / 0.000 / 0-1-13 (0.000) / 0 -> **0.000** |
| attack:naive_no_normalise | 0.866 / 0.821 / 1.000 / 1.000 / 12-1-1 (0.923) / 12 -> **0.533** | 0.872 / 0.825 / 1.000 / 1.000 / 12-1-1 (0.923) / 12 -> **0.536** |
| attack:naive_no_rename | 0.000 / 0.000 / 0.000 / 0.016 / 0-0-13 (0.000) / 0 -> **0.000** | 0.000 / 0.000 / 0.000 / 0.018 / 0-0-13 (0.000) / 0 -> **0.000** |
| attack:naive_phase_blind | 1.000 / 1.000 / 1.000 / 1.000 / 13-4-0 (0.867) / 13 -> **0.550** | 1.000 / 1.000 / 1.000 / 1.000 / 13-4-0 (0.867) / 13 -> **0.550** |
| attack:one_false_candidate | 1.000 / 1.000 / 1.000 / 1.000 / 13-1-0 (0.963) / 13 -> **0.872** | 1.000 / 1.000 / 1.000 / 1.000 / 13-1-0 (0.963) / 13 -> **0.872** |
| attack:one_missing_candidate | 1.000 / 1.000 / 1.000 / 1.000 / 12-0-1 (0.960) / 12 -> **0.846** | 1.000 / 1.000 / 1.000 / 1.000 / 12-0-1 (0.960) / 12 -> **0.846** |
| attack:per_variant_consequence | 1.000 / 1.000 / 1.000 / 1.000 / 12-4-1 (0.828) / 12 -> **0.535** | 1.000 / 1.000 / 1.000 / 1.000 / 12-4-1 (0.828) / 12 -> **0.535** |
| attack:provider_genotypes | 1.000 / 1.000 / 1.000 / 1.000 / 10-0-3 (0.870) / 10 -> **0.504** | 1.000 / 1.000 / 1.000 / 1.000 / 10-0-3 (0.870) / 10 -> **0.504** |
| attack:records_only_no_added_sites | 1.000 / 0.639 / 1.000 / 1.000 / 12-1-1 (0.923) / 12 -> **0.598** | 1.000 / 0.652 / 1.000 / 1.000 / 12-1-1 (0.923) / 12 -> **0.598** |
| attack:split_before_liftover | 0.997 / 0.911 / 1.000 / 1.000 / 12-0-1 (0.960) / 12 -> **0.804** | 0.997 / 0.909 / 1.000 / 1.000 / 12-0-1 (0.960) / 12 -> **0.803** |
| attack:zero_byte_files | gate: trio.GRCh38.vcf.gz is empty -> **0.000** | gate: trio.GRCh38.vcf.gz is empty -> **0.000** |

Seeds 3-8 were generated into scratch space and scored with the same verifier inside the tools image (no `make`,
no agent image): generation on attempt 0 for five seeds and on attempt 1 for
one (seed 8: "too few alleles fit into the inverted block"); oracle 1.000 and second solution 1.000 on all six;
l3_reference 0.273 and ignore_reference_differences 0.273 on all six; records_only_no_added_sites 0.598
(reference-difference F1 0.575-0.664); split_before_liftover 0.798-0.839 (F1 0.902-0.960);
inverted_block_unmapped 0.531-0.546; inverted_block_lost 0.611; per_variant_consequence 0.535;
codon_merge_ignoring_phase 0.681; provider_genotypes 0.488-0.504; naive_no_normalise 0.527-0.589;
naive_ignore_quality 0.550; naive_phase_blind 0.550. Inverted blocks of these seeds: 4561-6471 bp, 38-49 scored
alleles, 13-20 indels, six different genes. An earlier robustness run on seeds 3-10 with the version before the
inverted block (11 candidates) also gave 8 of 8 on attempt 0 and 1.000 for both solutions.

Instance numbers (seed 1 / seed 2): input records 2898 / 2793; alleles on GRCh38 3228 / 3095; positions
where the assemblies differ 178 / 178, of which without a record 99 / 95; scored reference-difference alleles
183 / 180; inverted block 4929 bp in PPM1F with 46 scored alleles (21 indels) / 4584 bp in CIMAP1B with 44 scored alleles (14 indels); gray zone 5 / 4
alleles (as at level 3: provider alleles in 1-60 bp chain gaps that +liftover places and Picard rejects). Input
size 2.5 MB, agent image 1.96 GB (same tool list as level 3), generation 18 s for the final version on a loaded
host (14 s before the inverted block was added), oracle 0.09 / 0.11 min in `make`, whole validation 5.0 min with
Harbor (seed 1) and 6.1 min without (seed 2, host under load).

Measured outside the battery (scratch copy of the oracle, seeds 1 / 2): Picard LiftoverVcf with
`RECOVER_SWAPPED_REF_ALT=true` after splitting, plus the added sites and rule 5, scores 0.702 / 0.708
(reference-difference F1 0.828 / 0.832, unmapped F1 0.703 / 0.720, candidates 12-0-1): the two-allele swaps are
recovered, the three-allele records are rejected and the three-allele de novo is lost.

Verifier behaviour on equivalent encodings of the seed-1 oracle output (run by hand, all 1.000): phased genotypes
with the alleles in the other order, `GT`-only FORMAT, extra records that nobody carries, the two gray positions
dropped or written as one multi-base record, `*_variant` consequence names and upper-case models. A truncated
BGZF file scores 0 (gate); every genotype set to `0/1` with the planted candidate list scores 0.73; candidate rows
with chromosome `22` instead of `chr22` score 0.35 (the key is literal, as at level 3).

Reading of the bands:
- Both correct solutions give 1.000 in every call-set category on eight seeds. Reference differences: "records
  re-expressed by the tool, no site added" gives 0.575-0.664, "multi-allelic records split before the lift-over"
  0.902-0.960, "records rejected" 0; the band 0.80-0.97 gives nothing to the first and part of the credit to the
  second. Inverted block: 1.000 or 0.000; nothing in between was observed.
- Candidates: one false row 0.872, one missing row 0.846 (both below 0.90); every "ignore one element" variant loses
  the candidate component completely or keeps at most 41% of it (`split_before_liftover`: 69%).
- The highest careless variant is `split_before_liftover` (0.80-0.84): it gets the three-allele de novo wrong and
  loses a third to two thirds of the reference-difference credit. Its ceiling is 0.89, not 0.75, because it is a
  partial handling of P_reference_differences, not "ignore the element".

## Deviations from the level-4 spec

1. **Nothing from the spec was dropped.** All four elements are in and validated.
2. **13 candidates instead of the 8-12 target**: making the inverted block decisive needed two candidates inside it
   (with one candidate there, the scoring formula gives 0.76 for "block silently lost", above the 0.75 ceiling;
   computed, not measured; with two it is 0.611 measured), and the three-allele site got
   its own candidate so that mishandling it changes the list and not only a call-set fraction.
3. **"Two harmless changes that give a stop" does not exist** in the standard genetic code (brute force over all
   codons: 0 pairs of individually synonymous single-base changes that together give a stop). Used instead: silent
   alone + missense partner -> stop or a third amino acid (18 stop combinations), and for "rescue" the only four
   pairs of individually non-synonymous changes that are silent together (Ser TCT/TCC <-> AGT/AGC).
4. **Phase "follows from the parents' genotypes"** is written as an explicit rule (exactly one parent carries the
   allele; homozygous = both copies; otherwise unassigned and assessed only together with homozygous variants). The
   stricter Mendelian deduction (mother 1/1, father 0/1 -> maternal) is deliberately not part of the rule; no planted
   pair depends on it.
5. **No de novo reference allele** (fairness note of the spec): the three reference-difference candidates are a
   homozygous-recessive ALT, an inherited ALT and a de novo third allele that is non-reference in both assemblies.
6. **Reference differences are mostly real**, not planted: 173 real + 5 planted. The shipped GRCh37 slice therefore
   differs from the public hg19 by five substitutions and one inversion per instance (seed dependent). The
   instruction calls the two FASTA files "exactly as they are in use at the provider and at the laboratory".
7. The level-3 trap attacks `naive_no_liftover`, `naive_no_rename`, `naive_default_pipeline` still run the level-3
   script (they score 0 anyway); `naive_no_normalise`, `naive_ignore_quality`, `naive_phase_blind` run the level-4
   reference with exactly that trap ignored.
8. "Candidates taken from the provider's coordinates" is covered in both readings: `provider_genotypes`
   (inheritance patterns read from the genotypes as written against the provider's reference) and
   `candidates_at_provider_coordinates` (right list, provider positions/alleles).

## Known limitations and open doubts (level 4)

- **The tools handle most of the new call-set cases when used a certain way.** bcftools +liftover run on the raw
  records swaps REF/ALT, adds a new reference allele at three-allele sites and lifts the inverted block; Picard
  (after splitting) lifts the inverted block too, which is why `l3_reference` keeps full credit for it (measured
  F1 1.000). What no tool produces are the positions without a record (99 of the 183 scored reference-difference
  alleles on seed 1, and the two de novo cases that hang on them). So the inverted block is a realism element that only bites a home-made
  lift-over, and the hardness of P_reference_differences rests mainly on the added sites.
- **`split_before_liftover` is a sharp edge.** `bcftools norm -m-any` before the lift-over is a common habit and is
  what the level-3 reference did. At three-allele sites it writes "the other ALT" as the provider's reference
  allele, which becomes a wrong ALT copy after the swap. The instruction defines `GT` by copies and the provider
  record shows `1/2`, so I consider it fair, but it is the most likely way for a strong agent to end at 0.80-0.84.
- **Positions without a record depend on one sentence of the instruction** (the provider's completeness statement)
  and on the "how the family differs from GRCh38" definition. If the pilot shows agents reading the deliverable as
  "the provider's records, lifted", that sentence is where to look first.
- **`added_sites_missing_from_truth` scores 0.88**, only 0.02 below the pass mark: a perfect candidate list with a
  call set that lacks every added site would almost pass. It cannot be produced by an analysis (the candidate list
  depends on two of those sites), but the margin is thin by construction (weight 0.12).
- **Rescue pattern has little diversity**: the trans pair, the cis pair and the "hidden rescue" decoy are always a
  serine codon with changes at codon positions 1 and 2 (adjacent bases), because the genetic code offers nothing
  else. Whether the pair is reported depends on phase, so the pattern alone does not give the answer away, but it
  marks where to look.
- **Allele balance at three-allele sites.** After the lift-over the GRCh38 reference allele of an "added reference"
  record has no read count (`.` in +liftover output); the reference solutions treat it as 0 reads, which is what
  the data say. A proband with two different non-reference alleles then has balance 1.0 for each. To keep this
  convention away from the candidate list, real assembly differences within 30 bp of an exon are always "no record"
  sites, and the planted three-allele candidate has a GRCh38-reference copy in the child.
- **Homozygous variant with a partner on one copy.** Rule 5 says it must qualify on both copies and that differing
  classes are joined by `&`; no planted candidate exercises this (the homozygous partners are common alleles), so
  the `&` branch of the two solutions is untested against a scored case.
- **Unassigned variants.** The oracle puts unassigned SNVs into extra per-copy runs so that two of them within 3 bp
  never share a run; a codon split by an intron is not covered by that distance rule (no planted codon is split).
- **Duplicate records are counted as standard false positives whatever their position** (carried over from the
  level-3 verifier). `split_before_liftover` produces 18 duplicates at three-allele sites on seeds 1 and 2; they lower the
  standard F1 to 0.997 (still full credit) instead of the reference-difference F1. No effect on a correct
  submission, but the per-category diagnosis is slightly off in that case.
- **Indels that span a differing base are not planted** (nothing is planted within 50 bp of a real difference or
  60 bp of a planted one). The second solution skips such a record and counts it; the oracle would rely on
  +liftover. Untested territory by design.
- **Gray zone additions.** Real assembly differences within 2 bp of each other (one pair, GRCh38 23,981,469 and
  23,981,471) are not scored, because one record or two is a convention. Real differences 9 bp apart are scored as
  two single-base records.
- **The 173 real differences are the same positions in every seed** (which of them carry a record, and with which
  genotypes, is drawn per seed). They are computable from the shipped FASTA files and the chain, so nothing leaks,
  but an agent trained on many seeds sees the same hidden sites again.
- **Seed 2**: all 33 checks of `make --no-harbor` ok (`results/trio_s2_l4.log`). The build had to wait from 01:32
  to 01:47 UTC because the shared host had less than 8 GB free; in the meantime the same checks were run without
  building images (`results/trio_l4_emu_validate.py`, log `results/trio_s2_l4.no_image_build.log`), with identical
  scores. Following the level-3 convention only seed 1 is kept under `tasks/`; the seed-2 task directory was
  removed after validation (the generator is byte-deterministic: seed 1 regenerated into a scratch directory is
  `diff -r`-identical to `tasks/trio-liftover-s0001-l4`).
- **Not piloted** (by instruction). The level-4 difficulty claims rest on the attack battery, not on agent runs.
- Level-3 limitations that still apply: consequence vocabulary of background variants across annotators, the
  unvalidated index file, real gene symbols.
