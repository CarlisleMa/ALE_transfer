# Family spec: `xlms_model_selection` - choosing among predicted complex models with XL-MS restraints

ALE subdomain: 3.1 Biomolecular Structure & Design. (ALE's ZDOCK task scores docking poses against
a *given* native with given formulas; here the native is hidden and experimental restraints in a
different numbering scheme are the only evidence.)

## Scenario
A structural biology group has a set of predicted models of a two-protein complex from a structure
prediction run, plus lysine-lysine cross-links from an XL-MS experiment. Evaluate every model
against the cross-links and select the model that best explains the data.

## Truth model
- Hidden native: a real heterodimer from RCSB with enough lysines on both chains (choose by seed
  from a short list; cache downloads). Never shipped.
- Models: 20-30 rigid-body variants of the native (ligand chain rotated/translated relative to the
  receptor) spanning near-native (interface RMSD < 2 A) to wrong, each written as a prediction
  output: chains relabelled/reordered per model, residues renumbered from 1, per-model confidence
  values that do NOT track quality (the top-confidence model is wrong).
- Cross-links: sampled from lysine pairs within range in the native (intra- and inter-protein)
  plus ~10-15% false positives, reported with protein accessions and positions in full-length
  UniProt numbering (offset from the model numbering, shifted again after an unresolved segment).
- Truth per model: exact satisfied counts under the correct residue mapping and the stated distance
  rule; interface RMSD / DockQ to the hidden native.

## Inputs
`models/model_XX.pdb` (+ `ranking.json` with confidences), `crosslinks.csv`, `sequences.fasta`
(full-length sequences with accessions), `xlms_notes.md` (cross-linker, distance rule, FDR of the
identification). Tools: python, biopython, numpy; optionally a sequence aligner.

## Outputs
`restraints.tsv` (`model, inter_satisfied, inter_total, intra_satisfied, intra_total`),
`selection.json` (`best_model`, `ranking` of all models best to worst), `mapping.tsv` (for every
cross-linked position: accession, position, model chain, model residue number in model_01).

## Traps
numbering (offset + gap); chain identity differs per model; misleading confidence; false-positive
cross-links (the best model does not satisfy all of them).

## Scoring sketch
gates: files parse, all models present. Components: exact per-model satisfied counts (0.35);
`best_model` is near-native by the hidden interface RMSD (0.35); ranking agreement with the hidden
quality order (0.20); residue mapping (0.10). Naive variants: positions used as model residue
numbers; chain letters assumed constant; best = highest confidence; intra+inter pooled.
