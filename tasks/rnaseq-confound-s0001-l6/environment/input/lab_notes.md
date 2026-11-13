# Wet-lab notes: compound vs vehicle, primary cells

- Cells: primary cells from 12 unrelated donors (one culture per donor), passage 2.
- Treatment: compound at 1 uM or vehicle (0.1% DMSO) for 24 h; one condition per donor.
- RNA: column purification with on-column DNase; RIN >= 8.9 for all samples.
- Gene panel: the transcripts in `reference/` (targeted panel plus sex-chromosome control genes).
- Library prep: poly-A selection, unstranded. Two prep sessions because the first kit lot ran out: 2024-03-05 (6 libraries, kit lot A, operator 1) and 2024-03-19 (6 libraries, kit lot B, operator 2, different PCR cycler). The session of each library is in the sample sheet.
- Sequencing: all libraries pooled equimolar on one lane, 50-cycle single read. The core delivered adapter-trimmed reads; read names and base qualities were stripped to placeholders by their export tool.
