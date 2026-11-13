# Rare-disease trio: family call set against GRCh38 and haplotype-aware SOP candidate review

Express a provider's joint trio call set against GRCh38 (including the positions where the two reference assemblies differ), report the alleles in unaligned sequence, and return the candidate causal variants of the proband under a lab SOP with per-copy codon consequences.

- family: `trio_liftover`  seed: `2`  level: `4`
- ALE subdomain: 3.2 Genomics & Sequence Analysis
- software: bcftools, bcftools +liftover, bcftools csq, samtools/htslib, Picard, Python/pysam
- a-priori difficulty: 0.701 (last-exam)  |  pass threshold: 0.9
- capabilities trained: verifying the genome build of delivered data, cross-build lift-over with reject handling, re-expressing genotypes where two reference assemblies differ, variant normalisation before database lookup, trio phasing by transmission, haplotype-aware consequence prediction, trio inheritance reasoning under a written SOP, exact output contracts

## Undisclosed decision points (not shown to the agent)

- **build_mismatch** - The provider call set is on hg19/GRCh37 but nothing says so; annotating it against the GRCh38 resources puts every variant in the wrong place (and one record then looks like a de novo missense variant of a panel gene).
  - discoverable by: the ##contig length 51304566 is the hg19 length of chr22 (GRCh38: 50818468); REF alleles match GRCh37.chr22.fa.gz at 100% and GRCh38.chr22.fa.gz at ~10-12%; bcftools csq/norm stop with 'the fasta reference does not match the VCF REF allele'; lab_sop.md says all lab resources are on GRCh38
  - trains against ALE failure mode: Understanding / unverified assumption about inputs
- **contig_naming** - The VCF names the contig '22' while FASTA, chain, GFF3 and frequency table use 'chr22'; Picard rejects every record (NoTarget), bcftools norm stops with 'sequence 22 not found'.
  - discoverable by: tool error/100% reject rate on the first run; bcftools view -h | grep contig
  - trains against ALE failure mode: Execution / tool-interface mismatch
- **representation** - A quarter of the alleles are right-shifted, padded or share a multi-allelic record (a causal SNV is the 2nd ALT of a record whose 1st ALT is a common deletion; a common frameshift is delivered right-shifted so a literal lookup in the frequency table misses it).
  - discoverable by: bcftools norm --check-ref/-f reports realigned and split records; multi-allelic ALT columns are visible in the VCF; the frequency table and the SOP state the canonical representation
  - trains against ALE failure mode: Approach / missing normalisation step
- **genotype_quality_decoys** - Four apparent de novo variants in panel genes rest on an unusable parental genotype (low DP/GQ) or on a proband allele balance of 9-19%.
  - discoverable by: rules 2 and 3 of lab_sop.md; the FORMAT fields DP, GQ and AD of the three samples
  - trains against ALE failure mode: Verification / skipping stated acceptance criteria
- **same_parent_pairs** - Two panel genes carry two rare damaging heterozygous variants that both come from the same parent (in cis): not a compound heterozygote.
  - discoverable by: rule 7 of lab_sop.md; parental genotypes in the VCF and the pedigree
  - trains against ALE failure mode: Understanding / domain knowledge gap
- **P_reference_differences** - Inside aligned chain blocks the provider's reference and GRCh38 differ at 178 single bases. A provider record whose ALT is the GRCh38 base is the opposite allele with complementary genotypes on GRCh38 (a default Picard run rejects it); 95 of these positions have no provider record at all because the family is homozygous for the provider's base, which is a homozygous ALT allele on GRCh38 that no lift-over tool can produce from the input; some sites carry a third allele. A homozygous-recessive variant, a compound-het partner and a de novo variant (third allele) sit on such sites.
  - discoverable by: the instruction defines the call set by meaning (how the family differs from GRCh38, given that the provider's file is complete against the provider's reference); Picard reports MismatchedRefAllele rejects although the positions lift; bcftools +liftover reports swapped/added reference alleles; comparing the two FASTA files through the chain shows the differing bases
  - trains against ALE failure mode: Approach / prior decoupled from truth (tool output taken as the deliverable)
- **P_minus_strand_block** - A 4584-bp segment around an exon of a panel gene is inverted between the provider's reference and GRCh38: the chain files carry it as a separate minus-strand chain. It holds 44 scored alleles (14 indels, which have to be reverse complemented, re-anchored and left-aligned again) including a compound-het pair (frameshift + missense).
  - discoverable by: the header line of the extra chain has '-' in the query strand column; Picard LiftoverVcf and bcftools +liftover handle it (INFO/FLIP); a position-only lift-over yields REF alleles that contradict GRCh38 for every record of the block
  - trains against ALE failure mode: Execution / coordinate conversion without strand
- **P_haplotype_consequence** - SOP rule 5 assesses single-nucleotide variants that share a codon on one chromosome copy of the proband together. Planted: two missense variants in one codon on different copies (a true compound heterozygote; merged they would be silent), two qualifying variants in one codon on the same copy that are silent together (their gene must not be reported), and de novo variants next to a base where the provider's reference differs from GRCh38 (silent alone and damaging on the real haplotype, or the reverse).
  - discoverable by: rule 5 of lab_sop.md states the definition; neighbouring records in the same codon are visible in any annotation table; bcftools csq without -l reports compound consequences for phased genotypes
  - trains against ALE failure mode: Approach / per-variant logic where the stated rule is per haplotype
- **C_dependent_stages** - Lift-over -> call set against GRCh38 -> phase assignment -> consequence -> inheritance model: a position without a provider record that is missed in the first stage changes the consequence of a de novo variant in the same codon and thereby the candidate list (one causal variant lost, one false candidate gained).
  - discoverable by: the candidate list is derived from the delivered call set; rule 5 says that every variant the proband carries counts
  - trains against ALE failure mode: Planning / early error surfaces only in the final deliverable
- **Q_gradient** - 13 scored candidates in 9 genes instead of 3-4 in 2-3 genes, so that the candidate components give partial credit per element. This is a scoring change, not an in-data control in the sense of AUTHORING.md's Q element.
  - discoverable by: n/a (no decision point)
  - trains against ALE failure mode: n/a

## Judge

**Where the truth comes from.** The family's alleles are planted on GRCh38 (left-aligned minimal,
planted number of ALT copies per member); the provider VCF is *derived* by expressing these
haplotypes against the provider's reference. For alleles at positions where both assemblies agree
this is the level-3 construction (mapped back inside ungapped chain blocks with identical flanks,
then renamed/right-shifted/padded/merged; scored only if Picard LiftoverVcf and bcftools +liftover
both reproduce the planted key and genotypes). New at level 4:

- *Reference differences.* Positions inside aligned blocks where the provider's reference base
  differs from GRCh38: the 173 differences that exist between the two real assemblies inside the
  shipped windows plus five planted into the shipped GRCh37 slice inside coding sequence. There the
  truth needs no tool: the planted allele pair of each member is written once against GRCh38
  (truth) and once against the provider's base (input). A site where all three members are
  homozygous for the provider's base has no input record and is a homozygous ALT allele on GRCh38;
  a provider record whose ALT is the GRCh38 base turns into the opposite allele with complementary
  genotypes; a site where everybody carries the GRCh38 base disappears; some sites carry a third
  allele.
- *Inverted block.* A 3-7 kb segment of the shipped GRCh37 slice is replaced by its reverse
  complement and both chain files carry it as a separate minus-strand chain. Alleles inside it are
  written on the other strand and left-aligned there; they are scored only if both guard tools
  bring them back to the planted GRCh38 allele.
- *Haplotype consequence.* The generator runs the SOP, including the per-copy codon rule, as code
  on the planted data and requires exactly the planted causal set; `bcftools csq -l` confirms every
  designed single-variant class.

**Hard gates (score 0).** Missing/unparseable file; VCF without exactly the three input samples;
no `.tbi`/`.csi` index; VCF larger than 10 MB compressed / 50 MB uncompressed; more than 200
candidate rows; `model` outside the three allowed values; candidate alleles that are not plain
nucleotide strings.

**Components.** In the three call-set components a true positive needs the literal (chrom, pos,
ref, alt) key *and* the planted number of ALT copies in all three samples; a submitted record that
is not a scored allele is a false positive of the category its position belongs to.

| weight | component | metric | full credit | zero |
|---|---|---|---|---|
| 0.08 | callset_standard | F1 over scored alleles at plus-strand positions where both assemblies agree | >= 0.985 | <= 0.90 |
| 0.12 | callset_reference_differences | F1 over alleles at positions where the provider's reference base differs from GRCh38 | >= 0.97 | <= 0.80 |
| 0.07 | callset_inverted_block | F1 over scored alleles inside the minus-strand block | >= 0.95 | <= 0.70 |
| 0.08 | unmapped | F1 of (pos_in, ref, alt) against the planted alleles in unaligned provider sequence | >= 0.95 | <= 0.50 |
| 0.45 | candidates | F1 of candidate keys against the planted causal set (13 variants) | 1.0 | <= 0.87 |
| 0.20 | candidate_annotation | fraction of planted causal variants reported with the right gene, consequence class and model | 1.0 | 0 |

Records in which no sample carries the ALT allele are ignored (they claim nothing). Keys are compared
literally because the instruction defines the representation.

**Why these tolerances** (measured on seeds 1 and 2). Both reference solutions reproduce every
scored allele of all three categories (F1 = 1.000) and the exact candidate list. The standard band
is the level-3 band. Reference differences: a solution that re-expresses the provider's records but
adds no site reaches F1 0.64-0.65, one that splits multi-allelic records before the lift-over 0.91,
one that rejects the records 0; the band 0.80-0.97 gives nothing to the first, a third of the credit
to the second and full credit only above both. Inverted block: handled gives 1.000, not handled 0;
the band only has to leave room for one or two of about 45 alleles. With 13 planted candidates one
false row gives F1 0.963 and one missing row 0.960; the 0.87-1.0 ramp turns both into a total below
the pass mark (0.87 / 0.85), so a clinical short list with a wrong entry still does not pass, while
two errors keep partial credit and the annotation component is linear in the number of correct rows.

**What is deliberately not scored.** INFO/FORMAT fields other than GT, record order, extra columns,
records without a carrier, and the gray zone: alleles in 1-60 bp chain gaps on which the two guard
tools disagree, and assembly differences at neighbouring bases (one record or two is a matter of
convention).

## Verifier attack suite (must score low)

- `constant_answer` (max 0.05): valid empty call set, no unmapped alleles, no candidates
- `every_panel_variant` (max 0.4): perfect call set, but every allele seen in the proband is listed as de novo
- `one_false_candidate` (max 0.89): planted truth plus one extra candidate row (just beyond tolerance)
- `one_missing_candidate` (max 0.89): planted truth with one causal variant missing
- `naive_no_liftover` (max 0.3): contig renamed, coordinates left on the provider's assembly and annotated as if they were GRCh38
- `naive_no_rename` (max 0.3): lift-over run on the provider's contig name as delivered: every record is rejected and reported as unmapped
- `naive_default_pipeline` (max 0.3): no lift-over, no normalisation, no quality rules, phase-blind
- `l3_reference` (max 0.75): the level-3 reference method on the level-4 inputs: Picard LiftoverVcf (records whose REF is not the GRCh38 base are rejected and reported as unmapped), no sites added, one consequence per variant
- `ignore_reference_differences` (max 0.75): P_reference_differences ignored, everything else as in the reference solution: records at positions where the assemblies differ are rejected, no sites added
- `records_only_no_added_sites` (max 0.75): the lift-over tool re-expresses the provider's records at positions where the assemblies differ, but positions without a record are not considered
- `added_sites_missing_from_truth` (max 0.89): planted truth without the alleles at positions that have no record in the provider's file (candidates as planted)
- `split_before_liftover` (max 0.89): multi-allelic records are split against the provider's reference before the lift-over (the other ALT then counts as the provider's reference allele, which turns into a wrong ALT copy where the reference base changes); everything else as in the reference solution
- `inverted_block_unmapped` (max 0.75): P_minus_strand_block ignored: only plus-strand chains are used, the records of the inverted block are reported as unmapped; everything else as in the reference solution
- `inverted_block_lost` (max 0.75): positions of the inverted block are lifted without reverse complementing the alleles; the records then contradict the GRCh38 reference and are dropped silently
- `per_variant_consequence` (max 0.75): P_haplotype_consequence ignored, everything else as in the reference solution: every variant gets the consequence it has on its own
- `codon_merge_ignoring_phase` (max 0.75): variants of the proband in one codon are always assessed together, whichever chromosome copy carries them
- `provider_genotypes` (max 0.75): candidates taken from the provider's coordinates: inheritance patterns are read from the genotypes as written against the provider's reference
- `candidates_at_provider_coordinates` (max 0.4): planted call set and planted candidate list, but the candidate rows carry position and alleles as in the provider's file
- `naive_no_normalise` (max 0.75): level-4 reference solution without splitting / left-alignment
- `naive_ignore_quality` (max 0.75): level-4 reference solution without the SOP genotype-quality and allele-balance rules
- `naive_phase_blind` (max 0.75): level-4 reference solution, any two rare damaging heterozygous variants in a gene called compound_het

## Must-pass probes (guard against false negatives)

- `planted_truth` (min 0.999): submission written directly from the planted truth
- `within_tolerance` (min 0.99): 0.7% of the scored records missing, one hg19-only allele missing, phased GT, UCSC contig name in unmapped.tsv, *_variant consequence names, shuffled rows, extra column
- `independent_solution` (min 0.99): second implementation: own chain lift-over, own left-alignment and own consequence prediction in Python (no Picard, no bcftools csq)
