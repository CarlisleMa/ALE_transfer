# Acquisition notes - experiment IMM-FC-012_PBMC_lineage (16-SEP-2025)

- Instrument: LSRFortessa X-20 (405/488/561/640 nm), FACSDiva 9.0, daily CS&T passed.
- Parameters recorded: FSC-A, FSC-H, SSC-A, eight fluorescence areas (see panel.tsv), Time. Detector voltages are stored in the files.
- Staining: surface antibody cocktail plus Zombie Aqua fixable viability dye, 1% PFA fixation. The reagent lots of the cocktail are listed in panel.tsv, the lot used for each control tube in controls.tsv.
- Material: cryopreserved PBMC, thawed on the day of staining; trypan blue viability after thaw between 80% and 92% depending on the donor.
- Tubes were vortexed briefly and run at medium flow rate (300-350 events/s) without straining through a filter cap.
- Controls recorded at the start of the session with the panel voltages: one unstained PBMC tube; antibody-capture beads (positive and negative beads in the same tube) for the seven antibody conjugates; Zombie Aqua on a 1:1 mix of live and heat-killed PBMC; an aliquot of the reference pool stained with CD8 PerCP-Cy5.5 alone (ctrl_PerCPCy55_cells.fcs); the HLA-DR FMO tube (reference pool, complete cocktail without HLA-DR PE-Cy7).
