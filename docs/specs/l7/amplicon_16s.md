# Level-7 spec: `amplicon_16s` - a negative control and a specimen that changed places

Why: gemini-3.8-flash scored 1.0 on level 3 and twice on the pointer-free level 6 of this family. The shape
that stopped four models on crispr_screen (level 6) was: the sample sheet is wrong for two libraries, no text
says so, and the data show it to anyone who checks sample identity before trusting the sheet. This family
has the same lever in its controls.

## New decision point (silent)
**M_control_swap.** The reads of one extraction blank and of one study specimen changed places before
sequencing (the index pairs of the two libraries were exchanged): the FASTQ pair filed under the blank's
`sample_id` holds the specimen's library and vice versa. The sample sheet, including `sample_type` and
`qpcr_16s_copies_per_ul`, is what the wet lab wrote down and is left unchanged. Choose a specimen with
enough biomass that its library is clearly a community and not reagent background.
- Naive outcome: contaminant identification that takes `sample_type` from the sheet treats a real community
  as a blank and a blank as a specimen. Nothing fails; the contaminant list is merely wrong (specimen taxa
  flagged, reagent taxa missed), and so is the profile the investigator derives from it.
- Expert practice: before negative controls are used they are looked at. A library labelled blank that is
  as rich as a specimen, clusters with the specimens and is not dominated by the taxa of the other blank,
  and a "specimen" that looks like the other blank and contradicts the relation between qPCR copies and
  reagent-taxon share, are found by the first ordination or bar plot. The two libraries are then used in
  the roles their content supports.
- Truth: unchanged in meaning. `asv.fasta` and `taxonomy.tsv` are the same. `asv_table.tsv` is defined by
  the reads of each file, so the planted counts of the two columns change places. `contaminants.txt` is
  defined by where the template DNA comes from (SOP item 4) and is the same set as before.
- Nothing in the instruction, the SOP or the notes may mention the exchange, sample identity, or checking the
  controls. No deliverable asks for sample types.
- `Trap.discoverable_by`: name the QC and give measured numbers (richness and reagent-taxon share of the
  two libraries against the other blank and the other specimens; position in the copies-versus-share relation).

## What stays
Everything of level 3 (data model, the five decision points, deliverables, scorer, tolerances). The text of
level 7 is produced by the pipeline: the family returns a level-7 build (`generate(seed, 7, workdir)`,
`level == 7`, `task_id` ending in `-l7`) whose instruction, SOP and notes are the level-3 texts unchanged, and
`alegen/descaffold.py` then removes the pointer sentences exactly as for level 6. Do not edit those texts.

## Battery (in addition to the level-3 battery, on the level-7 data)
- `l6_reference` (required by the validator at level 7): the level-3 reference solution unchanged (roles from
  the sample sheet) on the level-7 inputs. Must score at most pass - 0.15. If the scorer's weights make that
  impossible without changing the scorer, do not change tolerances silently: report the measured score, what
  would have to change, and stop.
- `sheet_roles_everything_else_right`: planted truth for sequences, counts and taxonomy, contaminants as the
  reference method finds them with the sheet's roles - attack, must fail.
- `drop_the_odd_blank`: recognises the odd blank but simply discards both libraries instead of using them in
  their real roles - decide from measurements whether this is a valid route (then a probe) or loses too much
  information (then an attack), and say which and why.
- The independent second solution must detect the exchange by a different statistic than the oracle.

## Measurements to report (two seeds)
The numbers behind `discoverable_by`; the score of `l6_reference` and which components it loses; how many
contaminant calls change with the sheet's roles; that levels 1-3 and 6 regenerate byte-identically.
