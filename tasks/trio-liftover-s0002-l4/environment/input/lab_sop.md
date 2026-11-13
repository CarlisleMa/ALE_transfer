# SOP-GEN-017 (rev. 5) - Review of trio call sets for candidate causal variants

Scope: germline joint call sets of an affected child ("proband") and both parents, restricted to the
gene panel of the referral. All interpretation resources of the laboratory (transcript set,
population frequencies) are kept on GRCh38 with UCSC-style contig names, and results are issued on
GRCh38. Every rule below refers to the call set as expressed against GRCh38: "ALT allele" is an
allele that differs from the GRCh38 reference, and "carries n copies" counts the chromosome copies
of a family member with that allele.

1. **Representation.** Variants are compared with laboratory resources and reported as one ALT
   allele per row, left-aligned and minimal.
2. **Usable genotypes.** A genotype is usable if `DP` >= 10 and `GQ` >= 20. A variant is reviewed
   only if the genotypes of all three family members are usable. Genotypes that follow from a
   provider's statement of completeness (positions without a record in the delivered file) are
   usable.
3. **Allele balance.** If the proband carries one copy of the ALT allele, the variant is accepted
   only if the fraction of reads supporting the ALT allele, `AD_alt / (AD_ref + AD_alt)` with the
   read counts of the ALT allele and of the GRCh38 reference allele, lies between 0.25 and 0.75
   (limits included).
4. **Population frequency.** A variant is rare if its allele frequency in the laboratory frequency
   table is below 0.001. A variant that is absent from the table is treated as having frequency 0.
5. **Consequence.** Consequences are taken relative to the laboratory transcript set and describe
   what the chromosome copy of the proband that carries the variant encodes.
   - A variant is assessed on its own, with one exception: if one chromosome copy of the proband
     carries two or more single-nucleotide variants in the same codon, these are assessed together.
     Each of them then has the consequence of the codon with all of these changes applied, compared
     with the reference codon (`synonymous` if the amino acid is unchanged, `missense` if it is
     changed, `stop_gained` if a stop codon results). Single-nucleotide variants in one codon that
     lie on different copies are assessed separately, and so are all other combinations of variants.
   - Which copy carries a variant follows from the trio genotypes. A variant of which the proband
     carries two copies is on both copies. A variant of which the proband carries one copy is on
     the copy inherited from the mother if the mother carries the allele and the father does not,
     and on the copy inherited from the father if the father carries it and the mother does not. In
     all other cases (neither parent or both parents carry the allele) the copy cannot be assigned;
     such a variant is assessed together only with variants of which the proband carries two copies.
   - All variants of the call set that the proband carries count for this, whatever their
     frequency, genotype quality or own consequence.
   - A variant of which the proband carries two copies qualifies only if its consequence qualifies
     on both copies (if the two classes differ, report both joined by `&`).

   Qualifying classes: `missense`, `stop_gained`, `frameshift`, `splice_acceptor`, `splice_donor`.
   Every other consequence (synonymous, in-frame insertion/deletion, UTR, intronic, intergenic, ...)
   does not qualify.
6. **Panel.** Only variants with a qualifying consequence in a gene of the referral's panel are
   reviewed.
7. **Inheritance models.** A rare variant that passes rules 2-6 is reported if it fits one of:
   - `de_novo` - the proband carries one copy of the ALT allele, neither parent carries it;
   - `homozygous_recessive` - the proband carries two copies, each parent carries exactly one;
   - `compound_het` - the proband carries one copy each of two different variants in the same gene,
     one of them inherited from the mother only (the mother carries one copy, the father none) and
     the other from the father only (the father carries one copy, the mother none). Both variants
     are reported.
   Variants that fit none of these models are not reported.
