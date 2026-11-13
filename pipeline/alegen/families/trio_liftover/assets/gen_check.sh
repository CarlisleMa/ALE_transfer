#!/bin/bash
# Generation-time guards, run inside the tools image (never shipped to the agent):
#   1. the planted GRCh38 alleles are in normal form (bcftools norm leaves them untouched)
#   2. bcftools csq consequences of every planted allele
#   3. the provider VCF, lifted by two independent tools (Picard LiftoverVcf, bcftools +liftover)
set -euo pipefail
mkdir -p /out /tmp/w && cd /tmp/w
RES=/res
bgzip -c /work/trio.vcf > /out/trio.vcf.gz && tabix -p vcf /out/trio.vcf.gz
bgzip -dc $RES/GRCh38.chr22.fa.gz > hg38.fa && samtools faidx hg38.fa && samtools dict hg38.fa -o hg38.dict 2>/dev/null
bgzip -dc $RES/GRCh37.chr22.fa.gz > hg19.fa && samtools faidx hg19.fa

bcftools norm -c e -f hg38.fa /work/truth.sites.vcf -Ou 2>norm_truth.log | bcftools query -f '%POS\t%REF\t%ALT\n' > /out/truth_norm.tsv
bcftools csq -l -f hg38.fa -g $RES/GRCh38.chr22.canonical.gff3 /work/truth.sites.vcf -Ou 2>csq.log \
  | bcftools query -f '%POS\t%REF\t%ALT\t%INFO/BCSQ\n' > /out/truth_csq.tsv

if [ "${RENAME:-1}" = "1" ]; then
  printf '22 chr22\n' > ren.txt
  bcftools annotate --rename-chrs ren.txt /out/trio.vcf.gz -Ou > renamed.bcf
else
  bcftools view -Ou /out/trio.vcf.gz > renamed.bcf
fi
bcftools norm -m-any -c e -f hg19.fa renamed.bcf -Ou 2>norm19.log \
  | bcftools annotate --set-id '%POS\_%REF\_%FIRST_ALT' -Oz -o split.vcf.gz
tabix -p vcf split.vcf.gz
bcftools query -f '%ID\n' split.vcf.gz > /out/split_ids.txt
bcftools query -l split.vcf.gz > /out/samples.txt

picard LiftoverVcf I=split.vcf.gz O=pic.vcf.gz CHAIN=$RES/hg19ToHg38.chr22.over.chain REJECT=pic_rej.vcf.gz R=hg38.fa \
  > picard.log 2>&1
bcftools view -S /out/samples.txt pic.vcf.gz -Ou | bcftools norm -m-any -c e -f hg38.fa -Ou 2>>picard.log | bcftools query -f '%ID\t%POS\t%REF\t%ALT[\t%GT]\n' > /out/picard.tsv

bcftools +liftover --no-version split.vcf.gz -Ou -- -s hg19.fa -f hg38.fa -c $RES/hg19ToHg38.chr22.over.chain \
    --reject plug_rej.vcf.gz --reject-type z 2>plugin.log \
  | bcftools sort -Ou 2>/dev/null | bcftools view -S /out/samples.txt -Ou | bcftools norm -m-any -c e -f hg38.fa -Ou 2>>plugin.log \
  | bcftools query -f '%ID\t%POS\t%REF\t%ALT[\t%GT]\n' > /out/plugin.tsv
cp norm19.log picard.log plugin.log csq.log /out/ 2>/dev/null || true
touch /out/done
