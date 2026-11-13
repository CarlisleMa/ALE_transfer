# Level-4 spec: `variant_review` - review on a locked-down workstation where IGV is the only reader (T + P)

Pilot at level 3 (seed 1): Claude Haiku scored 1.0 in 4 min 6 s with 13 tool calls, terminal only:
it read the BAM files with pysam and never opened IGV. The session-file check (0.05) can be
satisfied by hand-written XML.

## Goal of this level
The alignments can only be read through IGV. The instruction says so plainly. This is the one
family in which the professional GUI software is necessary, not optional.

## Environment
- Alignments are shipped as CRAM (reference-compressed) with the reference FASTA.
- The agent image contains IGV 2.19.8 (GUI, batch port, the bundled `igvtools`), the desktop layer
  and Python with the standard library only: no samtools, pysam, htslib, htsjdk outside IGV, no
  compiler toolchain that would make a home-made reader cheap. Check what the IGV bundle's Java
  runtime exposes (`jshell`, `javac`) and remove what is not needed to run IGV. Document every
  bypass you find and could not close; do not claim more than you verified.
- You own `pipeline/alegen/desktop/` this round. Add a small command-line front end to the same
  actions as the MCP server (`screenshot <png>`, `click`, `type`, `key`, `scroll`, `launch`,
  `windows`), installed in the image, so that a shell-only agent (and the pilot harness) can drive
  the desktop. Keep the MCP server working.

## Routes that must be shown to work (reference solutions)
1. GUI route: drive IGV with the desktop actions - navigate to each candidate, obtain the
   read-level evidence (e.g. the alignment track's "Export Alignments..." of the visible region,
   or the read pop-ups / "Copy read details"), save the session from IGV, take the snapshots.
2. IGV command-line route: `igvtools count` with per-base and per-strand output and its duplicate
   and mapping-quality switches, plus the batch port for snapshots and session.
The oracle uses one, the independent solution the other. If one route cannot be made reliable,
say so with the evidence and deliver the other one twice-checked; do not fake a route.

## Review content
Keep the SOP explicit (tags with precedence) and the level-3 classes. Add two or three classes for
which counting alleles is not enough (P), each planted with exact truth, for example: adjacent
substitutions carried by the same molecules that the SOP reports as one multi-nucleotide event;
support confined to molecules with an anomalous insert size or soft clipping; the variant present
at low level in the normal sample. Mixed rather than textbook-clean artefacts where the SOP
threshold still decides unambiguously (stay far from every threshold).

## Deliverables
`review.tsv` as at level 3 (plus what the new classes need), `review_session.xml` saved by IGV
(checked for the tracks, the reference and the SOP's display settings), `snapshots/<candidate>.png`
for every candidate (PNG, IGV-rendered size, not blank, distinct between candidates).

## Scoring
Calls, tags and allele fractions as at level 3 (>= 0.85 of the weight); session and snapshots
<= 0.15 together, and they cannot lift a failing review over the pass mark. Pass 0.9.

## Battery
`l3_reference` (the level-3 pysam solution has no reader: it must produce nothing useful);
counts-only review from `igvtools count` without the duplicate / quality switches; each new class
ignored; hand-written session and blank snapshots with a correct table (must lose exactly the
session/snapshot weight and no more); beyond tolerance. Probes as in levels 1-3.
