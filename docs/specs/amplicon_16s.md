# Family spec: `amplicon_16s` - 16S rRNA amplicon reads to exact sequence variants and abundances

ALE subdomain: 3.2 Genomics & Sequence Analysis (microbiome). No counterpart among ALE public tasks.

## Scenario
A microbiome core receives paired-end 16S V4 amplicon reads for a study batch that includes
negative extraction controls and a mock community, and must deliver denoised sequence variants,
an abundance table, taxonomy and the list of reagent contaminants.

## Truth model
- Real 16S sequences (NCBI RefSeq targeted loci, bacteria 16S; cache the download). Extract the V4
  amplicon by primer match (515F/806R). Plant a community of ~40 strains, including several pairs
  whose amplicons differ by 1-2 nt, with per-sample abundances drawn from the seed. Ship a small
  reference database (sequence -> genus) that contains the planted genera among distractors.
- Read simulation: paired-end reads with primers (degenerate positions), 0-3 nt heterogeneity
  spacers, mixed read orientation, quality decay towards the 3' end, substitution errors,
  chimeras formed from abundant parents (5-8% of reads), low-level cross-talk between samples.
- Reagent contaminants: a few taxa present at low level in every sample and dominating the
  negative controls. A mock-community sample with known composition.

## Inputs
`fastq/<sample>_R1.fastq.gz`, `_R2`, `samples.tsv` (sample type: study / negative_control / mock),
`primers.tsv`, `reference_db.fasta` (+ taxonomy), `library_prep_notes.md`, `sop.md` (definition of
an ASV deliverable, taxonomy assignment rule, what counts as a contaminant).
Tools: vsearch, cutadapt, python.

## Outputs
`asv.fasta` (primer-free, forward orientation), `asv_table.tsv` (ASV x sample read counts),
`taxonomy.tsv` (ASV -> genus), `contaminants.txt` (ASV ids).

## Traps
primers and spacers; mixed orientation; chimeras; strain-level variants (97% OTU clustering merges
them); reagent contaminants identified from negative controls.

## Scoring sketch
gates: files parse and are mutually consistent. Components: sequence recovery (abundance-weighted
recall of exact planted sequences and precision against non-planted/chimeric sequences) (0.35);
per-sample composition distance to the planted abundances (0.30); contaminant flags (0.15);
genus assignment (0.10); mock community composition (0.10). Naive variants: no primer removal;
single orientation; no chimera removal; 97% OTU clustering; no contaminant handling.
