# Assay notes - plate HCS-6608 (imaged 2025-02-21)

## Wet lab
- Cells: adherent reporter line (identity blinded), 96-well imaging plate with optical bottom, all 96 wells in use.
- Compounds: eight blinded test compounds, 8-point 3-fold series from 30 uM, one well per concentration; vehicle = DMSO (0.3% in every well); REF-INH = reference inhibitor of the marker pathway at 10 uM. Dispensed acoustically according to plate_layout.tsv, 48 h treatment.
- Incubation: 37 C / 5% CO2, standard polystyrene lid.
- Fixation 4% PFA, permeabilisation 0.1% Triton X-100.
- Staining: Alexa Fluor 488 conjugate against the marker protein, two PBS washes, wells left in 100 uL PBS.
- Nuclear counterstain: Hoechst 33342 in the last wash.

## Imaging
- Automated widefield microscope, 10x / NA 0.45 air objective, sCMOS camera with 2 x 2 binning (1.3 um per pixel), 16 bit.
- Channel 1: Hoechst (ex 385 nm). Channel 2: Alexa Fluor 488 (ex 475 nm). One field per well at the well centre, cropped to 320 x 320 pixels for transfer. Same exposure time for all wells of a channel.
- Focus: hardware autofocus on the plate bottom.
- Reference frames (flatfield/): average of 50 exposures of a uniform dye slide per channel and of 50 exposures with the shutter closed, recorded the same day.
- Files are raw camera frames: no illumination, background or any other correction was applied.

## Manual counts
- 8 fields of this plate were counted by hand by two technicians with the point tool of the image viewer, following SOP IA-07 section 1 (annotations/). Fields were picked to cover untreated wells of the inner and outer plate area and wells of the highest and of intermediate concentrations.
