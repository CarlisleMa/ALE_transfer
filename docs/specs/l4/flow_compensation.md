# Level-4 spec: `flow_compensation` - compensation that has to be right, not just done (P + Q)

Pilot at level 3 (seed 1): Claude Haiku scored 1.0 in 2 min 40 s with 20 tool calls. Your own
limitations list says why: populations are separated by 1.5 decades, so any reasonable gate lands
in band once the four decisions are made; the voltage change is on a scatter detector only.

## New elements (implement at least two; show the measured gap for each on two seeds)
- **P_embedded_matrix.** Every sample file carries a `$SPILLOVER` keyword written by the
  acquisition software from an older panel set-up. It is well formed and wrong for today's
  detectors by an amount that moves scored populations out of band. The notes state facts only
  (when the instrument was set up, when the controls were run).
- **P_voltage_change.** A fluorescence detector's voltage was raised part-way through the run and
  the affected single-stain control was re-acquired afterwards. Samples before and after need
  different matrices; one matrix for all samples fails for one group. Truth stays recoverable
  because the second control exists.
- **P_tandem_lot.** For a tandem dye the bead control was stained with another reagent lot than the
  samples and has a different spillover; a cell-based single stain with the sample lot is also
  provided. The notes record the lots as facts.
- **Q_dim_and_rare.** Add SOP populations without a clean valley: a dim activation marker whose
  positive gate the SOP defines on the matching FMO control (stated percentile, on compensated
  data), and a rare double-positive population (0.2-1% of parent) that under-compensation inflates
  through spreading. The FMO tubes are the in-data control: a correct analysis can see for itself
  that its gate holds.

## Deliverables
As at level 3. If matrices differ between sample groups, define in the instruction how
`spillover.csv` reports them (e.g. one file per acquisition block); the definition must be explicit.

## Scoring
Keep the level-3 structure; the new populations enter the frequency components with relative
tolerances set between your two solutions and the careless variants. Pass 0.9.

## Battery
`l3_reference`; one variant per new element (embedded matrix used; one matrix for all samples;
bead control used for the tandem dye; gates set on uncompensated or mis-compensated data, or
without the FMO); beyond tolerance. Probes as in levels 1-3.
