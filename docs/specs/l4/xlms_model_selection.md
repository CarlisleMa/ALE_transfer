# Level-4 spec: `xlms_model_selection` - cross-links that a straight line misjudges (P + Q)

Pilot at level 3 (seed 1): Claude Haiku scored 1.0 in 1 min 35 s with 15 tool calls. Once the
numbering is mapped, the task is a Euclidean distance count over rigid-body variants.

## New elements (implement at least two; show the measured gap for each on two seeds)
- **P_surface_distance.** The notes give the distance criterion as a solvent-accessible surface
  distance (the cross-linker cannot pass through the protein), with the definition written out
  (probe, grid, which atoms, cutoff). Plant links whose straight-line distance is below the cutoff
  while the surface path is far above it (and keep all scored links at least 3 A away from the
  cutoff under both of your implementations). Models that satisfy most links by straight line but
  not along the surface must exist, so that the Euclidean count selects a wrong model. If an
  established tool for this distance can be installed offline (e.g. Jwalk), ship it as the
  professional software; otherwise the analyst implements it. Your two solutions must be
  independent implementations and agree on every scored link.
- **P_homo_oligomer.** The complex contains two copies of one protein. A link between two
  positions of that protein is satisfied if any assignment to the copies satisfies it; a link
  between the same position of both copies can only be between copies. Classifying by "same
  protein = within one chain" misjudges these links and changes the selection. Write the
  definition into the instruction.
- **Q_two_states.** The cross-links come from a mixture of two arrangements and no single model
  explains them. The deliverable names the pair of models that jointly explain the most links and
  assigns each link (state A, state B, both, neither). Only if you can keep one answer.
- Remove the level-3 by-products that made shortcuts possible where you can: models that differ
  inside the chains (so pooling matters), wrong poses repacked so that physical plausibility alone
  does not single out the near-native models.

## Deliverables and scoring
As at level 3 with the definitions extended; per-link classification (if two states or copies are
used) scored with MCC. Exact counts stay strict only where the 3 A margin holds. Pass 0.9.

## Battery
`l3_reference` (Euclidean counts with the level-3 mapping logic); one variant per new element
(straight-line distance; same-protein = same chain; single best model for a two-state data set);
selection by the prediction confidence; beyond tolerance. Probes as in levels 1-3.
