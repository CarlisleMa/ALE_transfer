You are the clinical bioinformatician of a rare-disease genetics laboratory. An external sequencing
provider delivered a joint variant call set for a family trio (an affected child and both parents).
The laboratory needs to know how the family differs from the laboratory's own reference assembly,
and the candidate causal variants for the gene panel of the referral.

Task directory:
- `/task`

Inputs (read-only) under `/task/input`:
- `trio.vcf.gz` (+ `.tbi`) - the provider's joint call set of the three family members, chromosome 22
  only, FORMAT fields `GT:AD:DP:GQ`. The provider states that the file is complete with respect to
  the reference sequence the provider worked with: at every position of that reference which is not
  covered by a record of the file, all three family members are homozygous for the base of the
  provider's reference.
- `pedigree.ped` - six-column PED file (family, individual, father, mother, sex, phenotype; phenotype
  2 = affected)
- `panel_genes.txt` - gene symbols of the referral's panel, one per line
- `lab_sop.md` - the laboratory's SOP for trio candidate review
- `resources/GRCh37.chr22.fa.gz`, `resources/GRCh38.chr22.fa.gz` (+ `.fai`, `.gzi`) - chromosome 22 of
  both human assemblies, bgzip-compressed, exactly as they are in use at the provider and at the
  laboratory (sequence outside the regions the laboratory works on is masked with N; coordinates are
  those of the full chromosome)
- `resources/hg19ToHg38.chr22.over.chain`, `resources/hg38ToHg19.chr22.over.chain` - UCSC liftOver
  chain files restricted to chromosome 22
- `resources/GRCh38.chr22.canonical.gff3` - the laboratory's transcript set on GRCh38 (Ensembl GFF3
  conventions, one transcript per gene)
- `resources/popfreq.GRCh38.tsv.gz` - population allele frequencies, columns `chrom`, `pos`, `ref`,
  `alt`, `af` (GRCh38, 1-based, one ALT allele per row, left-aligned minimal representation)

Available software (on `PATH`): `bcftools` (with its plugins, including `+liftover`), `samtools`,
`bgzip`, `tabix`, `picard`, Python 3.12 with pysam, numpy, scipy and pandas. There is no internet
access.

Your task:
1. Deliver the call set of the family on GRCh38 in the laboratory's canonical representation.
2. Report which input alleles lie in sequence that cannot be placed on GRCh38.
3. Apply the laboratory SOP and report the candidate causal variants of the affected child.

Required outputs under `/task/output`:
- `trio.GRCh38.vcf.gz` plus its index (`.tbi` or `.csi`) - bgzip-compressed VCF on GRCh38 with the
  contig named `chr22`, sorted by position. It describes how the three family members differ from
  the GRCh38 reference sequence: it holds every allele that differs from the GRCh38 reference and is
  carried by at least one family member, for all GRCh38 positions that the provided chain file
  aligns to the provider's reference (the aligned blocks of the chain; sequence that the chain does
  not align and N-masked sequence are out of scope), and nothing else. One ALT allele per record;
  every record in left-aligned, minimal (parsimonious) representation with `REF` equal to the GRCh38
  reference sequence; differences of single bases are written as single-base records. The same
  three samples under their original names. `GT` of a sample states how many of its two chromosome
  copies carry the ALT allele of that record (`0/0`, `0/1` or `1/1`); a copy that carries any other
  allele counts as `0` in that record.
- `unmapped.tsv` - tab-separated with header `chrom_in`, `pos_in`, `ref`, `alt`: one row per ALT
  allele of the input that lies at a position of the provider's reference which the provided chain
  file does not align to GRCh38 chromosome 22, with chromosome, position and alleles as written in
  the input file. An input allele at an aligned position is never listed here, whatever it turns
  out to be on GRCh38. Header only if there is none.
- `candidates.tsv` - tab-separated with header `chrom`, `pos`, `ref`, `alt`, `gene`, `consequence`,
  `model`: one row per candidate variant under the SOP. `chrom`, `pos`, `ref`, `alt` exactly as in
  your `trio.GRCh38.vcf.gz`; `gene` is the gene symbol as written in `panel_genes.txt`;
  `consequence` is the qualifying class that the SOP assigns to the variant; `model` is one of
  `de_novo`, `compound_het`, `homozygous_recessive`. Header only if there is no candidate.

Constraints:
- Do not modify anything under `/task/input`; write deliverables only under `/task/output`.
- Do not ask for confirmation. Execute directly.
