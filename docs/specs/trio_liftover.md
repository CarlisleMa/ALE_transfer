# Family spec: `trio_liftover` - rare-disease trio, cross-build variant prioritisation

ALE subdomain: 3.2 Genomics & Sequence Analysis. Not a copy of any ALE public task.

## Scenario
A clinical genetics lab receives a joint trio VCF (PROBAND, MOTHER, FATHER) from an external
sequencing provider. The lab's interpretation resources are on GRCh38. The analyst must bring the
call set onto GRCh38 in canonical representation and return the candidate causal variants under
the lab's SOP (de novo, compound heterozygous, homozygous recessive) in a gene panel.

## Truth model (planted in GRCh38 space, then the *input* is derived from it)
- Real reference material: UCSC hg19 + hg38 chr22 FASTA, both liftOver chains, Ensembl GRCh38 GFF3
  for chr22. Keep 4-6 windows (total <= ~3 Mb) containing >= 40 protein-coding genes; N-mask the rest
  of the chromosome so coordinates stay real and files stay small. Ship ONE canonical transcript per
  gene in the GFF3 so "consequence" has a single answer.
- Background: ~2-4k SNVs/indels with population AFs; parental genotypes from HWE, child by Mendelian
  transmission; FORMAT GT:AD:DP:GQ with realistic depths. A population-frequency table (GRCh38,
  left-aligned, biallelic) is shipped.
- Planted causal set C (varies with seed): 1 confident de novo protein-altering variant in a panel
  gene; 1 compound-het pair (one LoF from one parent, one missense from the other) in another panel
  gene; optionally 1 homozygous-recessive. Design each from the CDS so the consequence class is
  known by construction, then confirm it with `bcftools csq` at generation time (drop on mismatch).
- Decoys (each must fail exactly one SOP rule): two rare damaging hets in one gene inherited from
  the same parent; "de novo" with parental DP/GQ below SOP; proband alt-allele balance below SOP;
  rare damaging variant in a non-panel gene; common (AF above cutoff) damaging variant whose input
  representation is non-normalised (so a string lookup in the AF table misses it); synonymous de
  novo; a record whose hg19 coordinate, if (wrongly) read as GRCh38, lands in a panel-gene exon.
- Input derivation: back-lift every planted GRCh38 variant to hg19, rename contig `chr22` -> `22`,
  right-shift a subset of indels, merge some sites into multi-allelic records (causal allele as 2nd
  ALT in at least one case), add records in hg19-only sequence (no GRCh38 alignment).
- Tool consensus guard: at generation time lift hg19->GRCh38 with two independent tools (e.g.
  CrossMap and Picard LiftoverVcf, then `bcftools norm`). A planted variant is *scored* only if both
  reproduce the planted GRCh38 key; hg19-only records are scored only if both tools fail to lift
  them. Everything else stays in the input but is a gray zone (neither TP nor FP).

## Inputs shown to the agent (`/task/input`)
`trio.vcf.gz` (header without assembly tag; contig length line present), `pedigree.ped`,
`panel_genes.txt`, `lab_sop.md` (numeric SOP: min DP/GQ, allele balance, AF cutoff, qualifying
consequence classes, definitions of the three inheritance models; says deliverables are on GRCh38),
`resources/` with GRCh37 and GRCh38 chr22 FASTA (+index), both chain files, GFF3, AF table.
Tools: bcftools, samtools/htslib, CrossMap, picard, python (pysam).

## Outputs (`/task/output`)
- `trio.GRCh38.vcf.gz` (+ `.tbi` or `.csi`): GRCh38, contig `chr22`, left-aligned, multi-allelics split.
- `unmapped.tsv`: `chrom_in, pos_in, ref, alt` of input alleles that cannot be placed on GRCh38.
- `candidates.tsv`: `chrom, pos, ref, alt, gene, consequence, model` (GRCh38, normalised;
  `model` in {de_novo, compound_het, homozygous_recessive}; consequence from the SOP vocabulary).

## Traps
build_mismatch (contig length + REF mismatches reveal hg19); contig_naming; representation
(non-normalised indels, multi-allelics) feeding into AF lookup and compound-het logic;
liftover_rejects; inheritance decoys.

## Scoring sketch (calibrate with measurements)
gates: files parse, three samples present, candidates <= 200 rows. Components: lifted call set F1 on
scored keys (0.20) + genotype concordance on matched sites (0.10); unmapped set F1 (0.10);
candidate variants F1 on keys (0.40); gene+consequence+model correct on true candidates (0.20).
Naive variants: annotate without lifting; lift but do not normalise/split; ignore parental quality
rules; phase-blind compound het.
