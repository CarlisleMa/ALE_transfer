# Level-4 pilot: reading of every result (Claude Haiku, two attempts per task: batch D, then batches E and F)

Purpose: decide for every lost point whether it is difficulty (a method, understanding or execution error of the
agent) or a task defect (an ambiguous or missing statement). Sources: the verifier's `details.json`, the submitted
files, the agent's final message, and the agent's working files left in its container. The agents' session logs were
not read.

| task | reward | time | what happened |
|---|---|---|---|
| rnaseq-confound-l4 | 1.000 | 2 min 35 s | Re-quantified with the right library type, failed the mixed library on the SOP's identity criterion, excluded (rather than exchanged) the two mixed-up libraries, all 8 audit verdicts and all 40 erratum verdicts right. |
| flow-compensation-l4 | 1.000 | 3 min 44 s | Ignored the header matrix, took the cell-based control for the tandem dye because the bead lot differs from the panel lot, measured the gain ratio 1.891 from the PE-Cy7 pair and scaled the other fluorochromes' spill into that detector, gated on the matching FMO tube, ran the SOP's acceptance check (6.0 and 23.3 units). |
| crispr-screen-l4 | 0.910 | 4 min 14 s | Counts, both log fold changes, calls (reference genes 10/10 and 0/30), erratum all right; exchanged the two mislabelled libraries. Lost 0.09 on one audit verdict: see below. |
| variant-review-l4 | 1.000 | 9 min 23 s | Operated the IGV desktop (9 screenshots viewed), exported the alignments of eight windows per file with the track menu's "Export Alignments...", applied the SOP to the exported SAM; saved the session and 68 images from IGV. All 68 calls, tags and allele fractions right. |
| trio-liftover-l4 | 0.442 | 11 min 10 s | Call set perfect in every category (standard 2997/2997, reference-difference sites 183/183, inverted block 46/46, unmapped 26/26). Candidate list: 7 right, 1 false, 6 missed; one of the 7 with the wrong consequence. See below. |
| xlms-model-selection-l4 | 1.000 | 14 min 8 s | Wrote the lattice path search from the written definition; all 17 models' counts and all 884 link statuses identical to the planted truth. |
| hcs-dose-response-l4 | 1.000 | 16 min 45 s | Tried Cellpose, found it under-counting dense fields by 10-20, built a threshold + watershed + size-split pipeline and tuned it on the 8 annotated fields until every field met the SOP criterion; all hidden wells and all 16 potencies right. |

First attempts: six of seven pass (pass mark 0.9), mean 0.907. Second attempts: see the end of this file.

## crispr-screen-l4: the lost audit point was a task defect (fixed)
The agent rated `sample_columns` as `defective`; the planted verdict is `ok`. The inherited count table has one
column without any counts (the Read-2 library, checklist item `read2_library`, planted as defective). The item read
"Every sample column of the inherited count table holds the counts of the FASTQ file that the sample sheet lists for
that sample." An empty column does not hold the counts of its FASTQ file, so `defective` is a defensible reading.
The intent was "no column was taken from another sample's file".
Fix (generator, both seeds re-validated): the item now reads "No sample column ... was counted from a FASTQ file other
than the one the sample sheet lists for that sample. (A column without any counts is the subject of `read2_library`,
not of this item.)", and the SOP states that an effect of one item's defect is recorded under the item that caused it.
The pilot result above is on the wording before the fix.

## trio-liftover-l4: a method error, not a task defect
All six misses, the false candidate and the wrong consequence are in genes on the minus strand (PPM1F, MAPK11,
CHCHD10, DDT, SLC25A1); every candidate in a plus-strand gene is right. The agent did not use `bcftools csq` for the
consequences but wrote its own annotator, because the SOP's per-copy rule goes beyond the tool's default. In the
agent's transcript table (left in its container) the coding sequence of every minus-strand gene starts with `TAC`
and ends with `ACT`/`ATC`/`ATT`: the sequence was complemented but not reversed. Its codon lookup therefore gave
`synonymous` for the missense partners of both frameshift pairs and for both variants of the adjacent trans pair
(so three compound-het pairs disappeared), `missense` instead of `stop_gained` for the de novo that sits next to a
homozygous position without a record, and it kept the decoy de novo.
The SOP defines the classes and the transcript set is given; the reference solution (bcftools csq with phased
genotypes) and the independent solution (own annotator) agree on all 13. What level 4 added did bite, but indirectly:
the tool's default no longer answers the question, the agent re-implemented the tool, and the re-implementation had
an error that a comparison with the tool on any single variant would have shown.

## What the passes say about the level-4 elements
- H (inherited project): not a difficulty for this model. Both audits were done by re-running or re-deriving, and
  the sound-but-odd decoy items were not flagged (rnaseq 8/8; crispr 7/8, the one miss being the defect above).
- Q (quality bar with a control): the agents used the controls exactly as intended (FMO check, annotated fields,
  reference genes) and iterated until the stated criterion held. A stated criterion plus a control turns the
  quality bar into a checklist item.
- P (prior decoupled from truth): none of the planted defaults was taken (header matrix, bead lot, straight-line
  distance, same-chain rule, per-record lift-over).
- T (IGV is the only reader): solved through the GUI in nine minutes.
- The only failure came from an unforced re-implementation error inside a long chain with many scored items.

## Second attempts (batches E and F; fresh containers, same prompt and cap)

| task | first | second | time (2nd) | what happened in the second attempt |
|---|---|---|---|---|
| rnaseq-confound-l4 | 1.000 | 0.824 | 3 min 16 s | Audit 8/8 and identity control right again. This time it exchanged the *condition* of the two mixed-up libraries but left their library-prep session as in the sheet: 31 false discoveries next to 77 of 79 planted genes (false-discovery proportion 0.287, full credit up to 0.12). Two retracted genes wrongly confirmed as a consequence. |
| flow-compensation-l4 | 1.000 | 1.000 | 11 min 55 s | Same reasoning as the first attempt (gain ratio 1.894, cell control for the tandem dye, acceptance check 12.6 and 28.7 units; it also showed that the bead control fails the check with 859). |
| crispr-screen-l4 | 0.910 | 0.363 | 4 min 7 s | Counts and arms right. The reported drug effect of the modifiers is attenuated by a factor of 2-7 (planted -2.59 reported as -0.38), 12 modifiers missed, 8 of 10 reference modifiers right although the agent's message says 10 of 10: the check was not made on the table that was submitted. Audit: `sample_columns` right with the new wording, but `stagger_reads` rated defective because the inherited script loses 2.5% of the reads to scaffold mismatches, which is not what the item states. Erratum 11 of 40 wrong as a consequence of the calls. The container was removed before I looked for the cause of the attenuation, so that cause is not established. |
| variant-review-l4 | 1.000 | 1.000 | 7 min 28 s | GUI route again (three 20-kb exports per file); all 68 rows right. |
| trio-liftover-l4 | 0.442 | 1.000 | 9 min 22 s | Own annotator again, this time with the minus strand handled correctly; all 13 candidates with the right model and consequence. |
| xlms-model-selection-l4 | 1.000 | 1.000 | 82 min 9 s | All counts and statuses identical to the truth again, but the path search was much slower than in the first attempt (two of the 17 model runs died and were repeated); 82 of the 90 minutes allowed. |
| hcs-dose-response-l4 | 1.000 | 1.000 | 26 min 29 s | A different route from the first attempt: Cellpose (cyto3, diameter 12) after threshold + watershed and LoG had failed the criterion on the annotated fields; all hidden wells and potencies right. |

Second attempts: five of seven pass, mean 0.884. Over both attempts: 11 of 14 pass, mean 0.896; four tasks pass twice
(flow, variant-review, xlms, hcs), three pass once (rnaseq, crispr, trio), none fails twice.

Reading of the three failures:
- None is caused by an ambiguous statement. The one ambiguity that the first round exposed (crispr `sample_columns`)
  was read correctly after the fix.
- rnaseq (second attempt) fails on a *level-3* decision point, the mixed-up pair: the lab notes say that the
  sample-sheet fields follow the original records, so a library whose identity is established takes its whole record
  with it, session included, and the session shows in the data. The agent wrote that it was unsure about the session.
  The level-3 author had listed this half-fix as a point a reviewer might see as a hidden requirement; I read it as
  fair (stated in the notes, checkable in the data) but it is the closest case in the set.
- All three failures are slips in one link of a long chain, and they hit the three tasks with the most scored
  conclusions (crispr: five deliverables; rnaseq: five deliverables; trio: thirteen candidates in nine genes). The four
  tasks that passed twice have one to three deliverables. With two attempts per task this is an observation, not a
  measurement.
- The same model solved every one of the seven tasks at least once: none of them is beyond it.
