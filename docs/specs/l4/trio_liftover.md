# Level-4 spec: `trio_liftover` - the lift-over cases the level-3 task left out (C + P + Q)

Pilot at level 3 (seed 1): Claude Haiku scored 1.0 in 3 min 23 s with 18 tool calls. Your own
limitations list says why: every mapped allele is an easy lift-over case, and consequences are
per variant.

## Idea
The family already chains lift-over -> canonical call set -> annotation -> inheritance models.
Level 4 makes the early stages matter for the final list by planting exactly the cases where
"run the tool" is not enough, with exact truth, because the truth is the family's haplotypes on
GRCh38 and the provider's file is *derived* from them.

## New elements
- **P_reference_differences.** Sites inside aligned blocks where the two assemblies carry
  different bases (plant them in the shipped FASTA slices if the real ones are too few or sit
  badly; record what you did). Derive the provider's records by expressing the true haplotypes
  against the provider's reference. Consequences, all exact by construction: a provider record
  whose ALT is the GRCh38 base (some genotypes turn into reference), a provider site with no
  record at all where the family carries the provider's base (a GRCh38 variant that the input
  never mentions), and a three-allele site. The instruction defines the deliverable by meaning
  ("every site where a family member differs from GRCh38 in the placed regions, given that the
  provider's file lists every site where one differs from the provider's reference"), not by tool.
- **P_minus_strand_block.** One inverted aligned block containing variants (alleles reverse
  complemented, indels re-anchored and re-normalised).
- **P_haplotype_consequence.** The SOP defines the consequence of a variant as the protein change
  of the child's haplotype. Plant pairs of variants in the same codon: in cis (phase follows from
  the parents' genotypes) the joint change differs from the two per-variant changes (e.g. a stop
  that is rescued, or two harmless changes that give a stop); in trans they are judged separately.
  At least one such pair decides whether a gene is on the candidate list.
- **Q_gradient.** More scored candidates (target 8-12 across the three models plus decoys) so the
  final list gives a gradient instead of three items.

## Fairness
Do not plant a causal variant whose status depends on a debatable convention (e.g. a de novo
reference allele). Every planted case needs one answer under the written definitions; if your two
solutions disagree on a case, fix the definition or drop the case.

## Scoring
Final candidate list >= 0.5 of the weight (variant identity, model, consequence); call set and
unmapped list <= 0.35; keep partial credit for the call set by category (easy alleles, reference
differences, inverted block) so the reward says where the lift-over went wrong. Pass 0.9.

## Battery
`l3_reference` (level-3 method: tool lift-over, normalise, per-variant consequence); one variant
per new element (reference differences ignored; inverted block dropped or not reverse
complemented; per-variant consequences); candidates taken from the provider's coordinates; beyond
tolerance. Probes as in levels 1-3, with a genuinely independent second lift-over.
