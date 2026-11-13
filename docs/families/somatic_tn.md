# Family `somatic_tn` - tumour/normal pairs from FASTQ to filtered somatic calls

ALE subdomain 3.2 (Genomics & Sequence Analysis). Code: `pipeline/alegen/families/somatic_tn/`
(`__init__.py` task assembly, `sim.py` simulator, `assets/` oracle, second solution and verifier,
`dev/` development harness that is not shipped in tasks).

## Scenario

A cancer-genomics core sequenced tumour and normal DNA of three patients on a small hybrid-capture
panel (13 targets, about 46 kb, 2 x 100 nt). The analyst receives twelve FASTQ files, the LIMS
sample sheet, the panel reference, the target BED and the core's hand-over notes, and delivers the
tumour/normal pair of every patient, one VCF of PASS somatic SNVs/indels per patient with tumour
allele fractions, and a summary. The work is long (index, six alignments with duplicate marking,
identity check across six libraries, three tumour-vs-normal comparisons, artefact review, three
normalised VCFs) and there is no turnkey somatic caller in the image: the comparison and the
filters have to be built from `samtools`, `bcftools` and `pysam`.

## Truth model (what is planted, and why it is correct)

Everything scored is decided in `sim.py` (host side, numpy only) before a read exists.

1. **Panel reference.** A 120 kb window of hg38 chr22 (seeded start between 18 and 48 Mb, no N,
   GC 36-60%), upper-cased. Twelve capture targets of 2.6-4.4 kb are placed where the UCSC soft-mask
   marks at most 35% repeat. From level 2 one 3.4-4.0 kb segment of a target is copied **exactly**
   into an inter-target gap (the reference grows by that length) and the copy, with 350 bp unique
   flanks, becomes a 13th target. The patients' genomes carry the same duplication.
2. **Eligible positions.** A planted site (outside the duplicated segment) must lie >= 200 bp inside
   a target (full depth), have +-60 bp free of soft-masked repeat and of any 35-mer that occurs twice
   in the reference (either strand), keep 150 bp from the duplicated segment and not touch a
   homopolymer of >= 5. A window with fewer than 15,000 such positions is rejected. Sites keep
   25 bp (germline), 30 bp (somatic SNV, artefact) or 60 bp (indel) from every planted site.
3. **Germline.** 110 population SNP sites with an allele frequency (Beta(0.7, 1.3), clipped to
   0.06-0.92); each patient's two haplotypes are drawn independently. Re-drawn until every pair of
   patients differs at >= 30 sites and every patient has >= 30 non-reference sites. Seed 1: 56/54/59
   non-reference sites, seed 2: 69/60/66.
4. **Somatic.** Per tumour: purity uniform 0.58-0.90; 13-16 clonal SNVs (allele fraction purity/2),
   6-8 SNVs of a subclone with cancer-cell fraction 0.45-0.65, 5-7 SNVs of a subclone with
   0.20-0.34, three SNVs at 2-4% (below the assay's range by design), 3-4 deletions (1-8 bp) and 2-3
   insertions (1-5 bp), clonal or first subclone. Substitutions follow a C>T-rich spectrum (36% C>T,
   20% T>C, 16% C>A, ...). Indels are placed where their representation cannot be shifted, and the
   planted form is checked to be the left-normalised one.
5. **Duplicated-segment variants (level >= 2).** 9-11 clonal SNVs per tumour in the interior of the
   duplicated segment, each in one randomly chosen copy, at least 470 bp from the copy edges (the
   largest insert is 450 bp), so no read pair that covers the site reaches unique sequence.
6. **Orientation artefacts (level 3).** One or two tumours (seeded) are "damaged". A damaged tumour
   gets 14-18 artefact sites of two chemistries (5-11 each): deamination (ref C: C>T in F1R2 pairs
   only; ref G: G>A in F2R1 pairs only) and oxidation (ref G: G>T in F1R2; ref C: C>A in F2R1). At
   each site 16-56% of the fragments *of the affected orientation* carry the substitution: nominal
   allele fraction 8-28%, realised 0.06-0.33 on seeds 1-3. Normals carry no artefact. In a damaged
   tumour five true SNVs of the first subclone (allele fraction 0.13-0.29) are planted in the same
   substitution classes; true variants are spread over both orientations like any other.
7. **Fragments and reads.** Per library and target, fragment midpoints are uniform over the target
   +-50 bp, insert length N(270, 45) clipped to 150-450, haplotype and pair orientation (F1R2 /
   F2R1) are fair coins, capture efficiency per target 0.85-1.15, depth 180x (tumour) / 90x (normal)
   +-5%. A fragment carries a germline allele if its haplotype does, a somatic allele with
   probability equal to the variant's allele fraction, an artefact allele with the site's rate if
   its orientation matches. 6.5% of fragments are sequenced twice and 1.5% three times (PCR
   duplicates with independent sequencing errors; `samtools markdup` flags about 8% of reads).
   Base qualities are F/:/, (95.5/3.5/1%) with substitution error 0.04/0.4/4%. Reads are shuffled.
8. **Sample sheet.** Library ids are random accession numbers. From level 2 the normal rows of two
   patients (seeded) are exchanged in `manifest.tsv`; tumour rows are always right.

`tests/reference/truth.json` holds, per patient, the true tumour and normal library, the planted
somatic variants with the **realised** fragment counts (alt fragments, alt fragments in the rarer
orientation, covering fragments), the mirror positions, artefact sites and germline genotypes.
Generation is deterministic: `truth.json` and all twelve FASTQ files of seed 1 level 3 have the same
md5 in two independent generator processes (the `make` run and a `dev/score_files.py` run).

### Scored set and gray zone (the rule for "what the depth supports")

A planted somatic variant is **scored** iff, in the simulated tumour library,
`alt fragments >= 8` and `alt fragments / covering fragments >= 0.07` and
`alt fragments in the rarer pair orientation >= max(3, 0.25 x alt fragments)`, and it is not a
duplicated-segment variant. Everything else that was planted as somatic is **gray**: a call is
neither a true nor a false positive, a miss is not a false negative.

- The instruction states a reporting range from 5% allele fraction; the floor is 7% so that a
  variant the analyst estimates one or two points low is still inside the stated range.
- The orientation condition removes true variants that happen to look like an artefact (for 8
  alt fragments P(rarer orientation <= 1) = 7%); no analyst could tell those from artefacts.
- Duplicated-segment variants: every supporting pair has two equally good placements, so the copy
  cannot be determined from the data.
- Measured (level 3): 31/33/30 scored variants on seed 1, 34/31/29 on seed 2, 31/32/35 on seed 3;
  lowest scored allele fraction 0.072 / 0.075 / 0.076, lowest alt-fragment count 12 / 11 / 12; the
  oracle needs 5 reads and 4%. Gray besides the 9-11 duplicated-segment variants: 3-6 per patient
  (the 2-4% variants and the lowest subclone when purity is low).

Every other PASS call inside the targets is a false positive: the patient's germline variants,
artefact sites, mirror positions, sequencing or alignment noise. Calls outside targets are ignored.

### Reference material

| file | source | sha256 |
|---|---|---|
| `.cache/somatic_tn/chr22.hg38.fa.gz` (12.3 MB) | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz (copied from `.cache/variant_review/`, checksum verified on every run) | `05f9d97d6fbfd08a44ca45b50837ca2ae9c471f35ba79dffec04d2cb5eaaf695` |

## Inputs, software, outputs

Inputs: `fastq/<sample_id>_R{1,2}.fastq.gz` x 6 libraries, `manifest.tsv` (`patient`, `sample_type`,
`sample_id`, `fastq_1`, `fastq_2`), `reference/ref.fa` + `.fai` (one contig `chr22_panel`, 123-124 kb;
no aligner index), `targets.bed`, `core_notes.md`. Measured input size: 21.6 MB (seed 1 level 3),
21.1 MB (seed 2 level 3), 20.3 MB (seed 3 level 3), 21.5 MB (seed 1 level 2), 18.8 MB (seed 1
level 1); 205,000-217,000 read pairs per level-3 instance, 49 bytes of FASTQ per read.

Software in the agent image: bwa 0.7.19, minimap2 2.30, samtools 1.25, bcftools 1.25 (with bgzip and
tabix), pysam 0.24.1, numpy/scipy/pandas, Python 3.12. `docker images` lists the agent image with
2.0 GB disk usage / 490 MB content.

Outputs (exact definitions are in the instruction): `pairing.tsv` (`patient`, `tumour_sample`,
`normal_sample`), `somatic/<patient>.vcf.gz` (PASS records, normalised indels, `INFO/AF` = fraction
of tumour fragments carrying ALT), `summary.json` (pair and PASS SNV / indel counts per patient).
The patient id is anchored to the tumour library in the output definition, because "normals
exchanged" and "tumours exchanged" give the same genotype evidence.

`core_notes.md` lists facts only: panel, library chemistry, run, requested depth, a table of the six
libraries (specimen type, DNA input, shearing batch, PCR cycles, index) and the pathologist's tumour
cell content. It does not mention the duplication, the pairing, or what a specimen type means for
the analysis. Damaged tumours are listed as FFPE sections, the others as fresh-frozen biopsies.

## Traps and how an analyst finds them

Numbers are from the oracle's alignments of the shipped instances (`dev/evidence.py`, host pysam on
BAMs copied out of the tools image).

### `sample_swap` (level >= 2)
What goes wrong: with the sheet's pairs, every germline variant the patient does not share with the
other patient is "tumour-only". Seed 1: 23 and 21 false positives in the two affected patients;
seed 2: 16 and 7; seed 3: 15 and 16.
Evidence: genotype concordance of each tumour with each normal at sites that vary among the normals.

| seed | sites | same individual | other individuals |
|---|---|---|---|
| 1 | 86 | 1.00, 1.00, 1.00 | 0.34 - 0.48 |
| 2 | 84 | 1.00, 1.00, 1.00 | 0.42 - 0.52 |

The same is visible after a naive run: many calls at allele fraction ~0.5 and ~1.0 in two patients
and none in the third.

### `duplicated_segment` (level >= 2)
What goes wrong: read pairs from the interior of the segment get mapping quality 0 and a random
copy; a clonal tumour variant of one copy shows up at both copies with half its allele fraction.
Counting reads of any mapping quality adds 9-11 false positives per patient (the mirror positions).
Evidence: fraction of reads with MAPQ < 20 per target.

| seed | the two targets that hold the copies | the other eleven targets |
|---|---|---|
| 1 | 0.78, 0.77 | 0.000 each |
| 2 | 0.74, 0.80 | 0.000 each |

The candidate calls there come in pairs at the same offset in the two targets.

### `orientation_artifacts` (level 3)
What goes wrong: 14-18 extra calls per damaged tumour at 6-33% allele fraction, not in the normal,
with good base and mapping quality.
Evidence: alt-supporting reads of a site tabulated by read-pair orientation, and for comparison by
alignment strand (F1R2/F2R1 | forward/reverse).

| seed, patient | artefact sites (first six shown) | scored true SNVs of the same tumour |
|---|---|---|
| 1, PT-0760 | 0/24 \| 15/9, 0/19 \| 10/9, 0/52 \| 30/22, 20/0 \| 13/7, 0/34 \| 15/19, 42/0 \| 25/17; all 14 sites: rarer orientation 0 | rarer-orientation fraction min 0.25, median 0.45 |
| 2, PT-0204 | 40/0 \| 18/22, 17/0 \| 8/9, 0/34 \| 12/22, 13/0 \| 8/5, 45/0 \| 22/23, 0/16 \| 7/9; 15 sites: 0, one site 1/35 | min 0.30, median 0.44 |

Alignment strands are balanced at artefact sites, so a forward/reverse strand-bias filter keeps
every one of them (`naive_strand_bias_only` scores exactly like "no orientation filter"). The
artefacts are of two pyrimidine-reference classes (seed 1: 9 C>T + 5 C>A; seed 2: 8 + 8) and the
same tumour has 11 (seed 1) and 17 (seed 2) scored true SNVs of those classes, so the spectrum is
a cue that something is wrong but not a way to tell the calls apart. The specimen column of the
notes (FFPE vs fresh-frozen) is a fact a knowledgeable analyst connects; the notes do not.

Not a planted trap, but real in the data: bwa draws a 1-bp deletion close to a read end as a
mismatch. On the first instance generated during development this produced one false call (6 reads,
AF 0.05, next to a planted deletion). It is a false positive by the instruction's definition; the
precision band tolerates one per patient at full credit.

## Scoring

Hard gates (score 0): missing/unparseable file; a patient missing or duplicated in `pairing.tsv`;
a VCF record on another contig, outside the reference, or whose REF differs from `ref.fa`;
`summary.json` without a `patients` object. PASS and `.` count as PASS; other FILTER values are
ignored. Multi-allelic records are split; every allele is left-normalised by the verifier before
matching (stdlib implementation, `tests/reference/ref.fa`).

| weight | component | metric | full credit | zero |
|---|---|---|---|---|
| 0.20 | pairing | every patient's (tumour, normal) equals the planted pair | all | any wrong |
| 0.30 | somatic_precision | **min over patients** of TP / (TP + FP) | >= 0.95 | <= 0.85 |
| 0.30 | somatic_recall | **min over patients** of TP / scored; x ramp(min precision, 0.50, 0.85) | >= 0.95 | <= 0.85 |
| 0.10 | allele_fraction | worst patient's median abs(AF - realised fragment fraction) over TPs | <= 0.03 | >= 0.08 |
| 0.10 | summary | fraction of patients whose entry equals their own VCF counts and pairing row | all | none |

pairing, precision, allele_fraction and summary are multiplied by ramp(min recall, 0.2, 0.5).
Pass threshold 0.9.

Why precision and recall instead of the spec's single F1: with 30 true variants, 9 mirror calls
lower F1 only to 0.87. Scored separately, the same 9 calls give precision 0.77, below the band.

### Measured values (in-container, `results/work/<task>/*/logs/details.json`, final code)

Columns: min precision / min recall / worst AF error; TP,FP,FN per patient in sheet order; score.

| submission | seed 1 level 3 | seed 2 level 3 |
|---|---|---|
| oracle (bwa, markdup, pysam fragment counts, MAPQ >= 20, both orientations) | 1.000 / 1.000 / 0.0114; 31,0,0 33,0,0 30,0,0; **1.0** | 1.000 / 1.000 / 0.0089; 34,0,0 31,0,0 29,0,0; **1.0** |
| independent solution (minimap2, no dedup, bcftools mpileup AD on orientation-split BAMs, bcftools call genotypes for the pairing) | 1.000 / 1.000 / 0.0073; 31,0,0 33,0,0 30,0,0; **1.0** | 1.000 / 1.000 / 0.0109; 34,0,0 31,0,0 29,0,0; **1.0** |
| planted truth | 1.000 / 1.000 / 0; **1.0** | 1.000 / 1.000 / 0; **1.0** |
| within tolerance (1 miss + 1 false call per patient, AF +-0.02) | 0.967 / 0.967 / 0.0102; **1.0** | 0.966 / 0.966 / 0.0125; **1.0** |
| naive: sheet pairing (2 pairs wrong) | 0.574 / 1.000 / 0.0114; 31,23,0 33,21,0 30,0,0; **0.263** | 0.680 / 1.000 / 0.0089; 34,16,0 31,7,0 29,0,0; **0.354** |
| naive: no mapping-quality filter | 0.756 / 1.000 / 0.0114; 31,10,0 33,10,0 30,9,0; **0.620** | 0.725 / 1.000 / 0.0089; 34,10,0 31,10,0 29,11,0; **0.593** |
| naive: no orientation filter | 0.682 / 1.000 / 0.0114; 31,0,0 33,0,0 30,14,0; **0.556** | 0.680 / 1.000 / 0.0089; 34,16,0 31,0,0 29,0,0; **0.554** |
| naive: strand-bias filter instead of orientation | 0.682 / 1.000 / 0.0114; 31,0,0 33,0,0 30,14,0; **0.556** | 0.680 / 1.000 / 0.0089; 34,16,0 31,0,0 29,0,0; **0.554** |
| shortcut: blanket spectrum filter (FFPE tumours, damage classes below AF 0.30 dropped) | 1.000 / 0.733 / 0.0115; 31,0,0 33,0,0 22,0,8; **0.700** | 0.963 / 0.765 / 0.0089; 26,1,8 31,0,0 29,0,0; **0.700** |
| naive: all three ignored (2 pairs wrong) | 0.484 / 1.000 / 0.0114; 31,33,0 33,31,0 30,23,0; **0.200** | 0.447 / 1.000 / 0.0089; 34,42,0 31,17,0 29,11,0; **0.200** |
| naive: tumour-only | 0.337 / 1.000 / 0.0114; 31,56,0 33,54,0 30,59,0; **0.400** | 0.305 / 1.000 / 0.0089; 34,69,0 31,60,0 29,66,0; **0.400** |
| truth + 30% false calls | 0.756 / 1.000 / 0; **0.620** | 0.756 / 1.000 / 0; **0.619** |
| truth - 25% of calls | 1.000 / 0.727 / 0; **0.700** | 1.000 / 0.724 / 0; **0.700** |
| sheet pairing, VCFs without records | 0 / 0 / -; **0.000** | 0 / 0 / -; **0.000** |

Extra instances (validated without Harbor): seed 3 level 3 (two damaged tumours): oracle 1.0,
independent solution 1.0, sheet pairing 0.343, no MAPQ filter 0.620, no orientation filter 0.514,
strand-bias only 0.514, spectrum filter 0.589, all ignored 0.200, tumour-only 0.400. Seed 1 level 2:
oracle 1.0, independent 1.0, sheet pairing 0.270, no MAPQ filter 0.624, all ignored 0.200.
Level 1 (seed 1): oracle 1.0 (AF error 0.0107), independent solution 1.0 (0.0106), the default
pipeline (sheet pairs, no mapping-quality filter, no orientation filter) 1.0, tumour-only 0.400.

### Tolerances and where they sit

- **Precision 0.85-0.95.** Correct routes: 1.000 on every instance, both routes, five analyst
  variants (below). One false call in 30 is 0.968, two are 0.938. Trap-ignoring routes: 0.31-0.76.
  The closest one is the mapping-quality trap (0.725-0.762), 0.09 below the zero bound. The
  generator asserts at least 9 mirror sites and a mirror share of at least 19% of (mirror + scored)
  per patient, so the worst case by construction is 0.81.
- **Recall 0.85-0.95.** Correct routes: 1.000 on every instance and variant. One miss in 29 is
  0.966 (full credit), two are 0.931 (credit 0.81, score 0.94), three are 0.897 (score 0.84). The
  band was 0.80-0.92 until the spectrum-filter attack below; it was tightened, not loosened.
- **Allele fraction 0.03-0.08.** Correct routes: 0.0073-0.0114 (0.0034-0.0040 without the read-end
  margin). Binomial noise of the planted, not realised, fraction at this depth would be about
  0.035, which is why the realised fraction is the reference.
- **Analyst-variation runs** (dev harness, tools image, final seed 1 and seed 2 level 3, all 1.0
  with 0 FP / 0 FN): oracle without read-end margin; oracle with 8 reads / 6%; oracle with 4 reads /
  3%; independent solution without its indel-proximity rule; the same additionally without BAQ.

### The spectrum-filter shortcut (found by attacking the verifier, then closed)

With one damage chemistry per instance and the 0.80-0.92 recall band, "drop every call of the
damage class below 20-25% allele fraction in the FFPE tumours" scored **0.91-0.97 on both seeds**
without looking at a single read (`dev/spectrum_filter.py`). After planting both chemistries and
artefacts up to 30% it still scored 0.94 on seed 1 (0.87 on seed 2). Three changes closed it: both
chemistries per damaged tumour, five true first-subclone SNVs in the damage classes, recall band
0.85-0.95. Now the filter loses 7-13 true variants in a damaged tumour at any cut-off that removes
the artefacts: 0.60 / 0.70 / 0.70 for cut-offs 0.25 / 0.30 / 0.35 on seed 1, 0.64 / 0.70 / 0.70 on
seed 2 (dev harness), and 0.700 / 0.700 / 0.589 in the battery (`shortcut_spectrum_filter`, ceiling
0.75, seeds 1-3). The score is exactly 0.70 whenever recall is at or below 0.85 with full precision.
The five forced variants alone do not guarantee that: a damaged tumour with 38 scored variants that
lost only those five would have recall 0.868 and score 0.754, just above the ceiling. The measured
instances lose 7-13 (recall 0.62-0.77) because the natural spectrum adds true variants of the same
classes; a seed where it adds none would trip the ceiling in `make`, not ship silently.

## Levels

1: no trap (sheet correct, no duplication, no artefacts); the default pipeline scores 1.0 and is a
must-pass probe. 2: `sample_swap` + `duplicated_segment`. 3: + `orientation_artifacts`.

## Deviations from the spec

- **No GATK4 / Mutect2** (image size; the brief asked for <= about 2 GB). Tools: bwa, minimap2,
  samtools, bcftools, pysam. There is no dedicated somatic caller; both reference routes build the
  tumour-vs-normal comparison from allele counts.
- **Orientation artefacts are read-pair-orientation artefacts**, not strand bias. This is the real
  signature and it is not removed by the filter most pipelines have by default.
- **Exact duplication, no paralog-specific mismatches.** The spec's "paralog-driven mismatches" are
  realised as true variants of one copy appearing at the other copy. A diverged paralog would have
  made the recoverability of each site depend on the nearest paralog-specific base; with an exact
  copy the gray zone is defined by construction.
- **Precision and recall (0.30 + 0.30) instead of one F1 (0.60)**; reason above.
- `somatic/<patient>.vcf.gz` carries `INFO/AF` (the spec did not say where the allele fraction goes).
- `summary.json` carries the pair and the PASS counts per patient.

## Known limitations and open points

- **Not piloted with a model.** Whether the three obstacles are found from the data by a small or a
  large model is untested. My expectation: the pairing is asked for explicitly and will be checked
  by most agents; `-q 20` is a habit, so the duplication trap will often be avoided without being
  noticed; the orientation trap is the one that separates.
- **The pairing deliverable points at the identity check.** That is the spec's output contract; the
  instruction says what to deliver, not that the sheet is wrong, but an agent reading "establish
  which normal ..." has been told where to look.
- **The artefacts are stronger and cleaner than real damage**: 100% one orientation, allele
  fractions up to 0.33, and five true subclonal variants per damaged tumour deliberately placed in
  the same substitution classes. Real FFPE/oxidation artefacts are graded and mostly below 5-10%.
  This was the price of closing the spectrum shortcut; a filter tuned here is not validated for
  real data.
- **The artefact chemistry does not vary with the seed** (both chemistries in every damaged tumour).
  What varies: which tumours (1 or 2 of 3), number of sites, the mix, positions, rates.
- **Other shortcuts I did not attack**: a filter on allele-fraction clusters (e.g. keep only calls
  near the clonal or subclonal peak), and a "minimum 2 reads in the rarer orientation" rule without
  a fraction test, which I expect to work and regard as a legitimate solution.
- **The mapping-quality trap has the smallest margin among the naive variants** (score 0.59-0.62
  against a ceiling of 0.75).
- **Sheet-pairing false positives vary with the seed** (7 to 23 per affected patient, depending on
  how many germline variants the two patients share). The generator asserts that at least one of
  the two affected patients has >= 12 unshared germline sites.
- **A submission with perfect calls and no usable AF scores exactly 0.90**, i.e. it passes. AF and
  summary are 0.10 each as in the spec sketch.
- **Earlier-artefact weight is 0.40** (pairing + AF + summary for a recall-complete but imprecise
  call set; tumour-only calling scores 0.40), slightly above the ~0.35 guideline.
- **Recall band**: three misses in a 29-variant patient fail the task (0.84). Every correct route
  I ran has zero misses, but an analyst with a conservative allele-fraction cut-off near 10% would
  fail although the instruction's stated range starts at 5%.
- **No copy-number changes, no LOH, no tumour-in-normal contamination, no germline indels, no
  indel sequencing errors, no adapter read-through.** Normals show 0-1 alt reads at somatic sites.
- **Duplicates carry no trap**: the independent solution does not remove them and scores 1.0.
- **Seeds whose window is repeat-rich are skipped inside `build_reference`** (the seed then uses a
  later window); a seed for which planting or one of the generation-time assertions fails raises at
  generation (seen once, before the site spacing was reduced; seeds 1-8 generate at level 3).
- The duplication is exact and several kb long, which is rarer in real genomes than diverged
  segmental duplications; the reference is called a panel reference to keep that plausible.
- Harbor (G6) was run for seed 1 level 3 only, as the contract asks.
