# Level-4 spec: `hcs_dose_response` - segmentation that has to be validated (Q + P)

Pilot at level 3 (seed 1): Claude Haiku scored 1.0 in 14 min 50 s with 35 tool calls. It found
Cellpose slow on CPU, switched to classical segmentation and still scored full, as your own
measurements predicted ("Cellpose and classical methods both count these images almost exactly").

## New elements (show the measured gap for each on two seeds)
- **Q_validated_counts.** Make the images hard enough that the counting method matters: dense
  touching and overlapping nuclei without an intensity dip, textured chromatin that a watershed
  over-splits, debris and apoptotic fragments that the SOP does not count, nuclei cut by the field
  border under a stated counting rule, nuclear size that changes with treatment (a fixed area
  filter miscounts treated wells). Use real nuclear texture if you can do it with exact truth
  (composite crops onto synthetic fields); say what you used. Ship the lab's manual point
  annotations for a handful of fields (`annotations/`), spanning sparse to dense and treated to
  untreated. The SOP states the acceptance criterion: agreement with the manual counts within a
  stated relative error on the annotated fields before the pipeline is used. A default threshold +
  watershed must fail that criterion visibly; at least one obtainable method (tuned Cellpose
  and/or a carefully built classical pipeline) must meet it. Meeting it on the annotated fields
  must be close to sufficient for the hidden wells.
- **P_potency_shape.** Curves for which the default four-parameter fit with free plateaus gives a
  wrong potency: partial responders whose potency the SOP defines on the absolute response scale,
  and a compound whose marker signal per cell rises only because cells are lost. Stay with one
  answer under the SOP's written definitions.
- Optional **I_time_budget**: fields large enough that running the model on everything with
  default settings does not fit the time limit on 4 CPU cores, so the analyst has to plan
  (resolution, batching). Only if the oracle stays within 20 minutes.

## Deliverables and scoring
As at level 3, plus `validation.tsv` (`field`, `manual_count`, `pipeline_count`) for the annotated
fields. Count accuracy on the hidden wells is scored on its own component with the band between
your two correct solutions and the default pipeline; potencies keep their weight. Pass 0.9.

## Battery
`l3_reference` (both level-3 routes: the classical pipeline and untuned Cellpose - report both
numbers); one variant per new element; counts tuned to the annotated fields only by a global
scale factor (must fail: the error is not a constant factor); beyond tolerance. Probes as in
levels 1-3.
