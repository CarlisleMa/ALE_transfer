# SOP IA-07 (v5): image analysis of two-channel potency plates (nuclear stain / cytoplasmic marker)

This SOP applies to every plate of this assay format. It defines what is reported and which
criteria a result has to meet, not how a particular plate has to be processed.

## 1 Readouts (per well)
- **R1, cell number**: the number of cells in the imaged field. One nucleus is one cell; nuclei
  that touch or overlap each other are separate cells. Two rules decide what is counted:
  - *Fragments are not cells.* Nuclear fragments, apoptotic bodies and other stained debris are
    compact objects of a few micrometres (the largest are about 6.5 um across). Intact nuclei of
    this cell line measure at least 11 um in their shortest dimension, with or without treatment.
  - *Counting frame.* A nucleus that is cut by the top or the left edge of the image is counted.
    A nucleus that is cut by the bottom or the right edge of the image is not counted (also when
    it is cut by the top or the left edge as well).
- **R2, marker per cell**: for each counted cell, the mean marker-channel intensity of its
  cytoplasm that is due to the stain bound in the cell. Marker-channel signal that does not come
  from the cells is background and is not part of R2. The well value is the mean over the counted
  cells of the field.

## 2 Validation of the counting procedure
On every plate the lab counts a number of fields by hand under the rules of section 1 (point
annotations). Before a counting procedure is used on a plate, its counts are compared with these
manual counts. The procedure is fit for use when, for every annotated field,

    | procedure count - manual count |  <=  max(3 cells, 10 % of the manual count)

A procedure that does not meet this on every annotated field is not used: it is changed and
compared again. The comparison is reported for the procedure that was finally used. Annotated
fields are counted by the procedure like every other field; manual counts are not entered into
the well table.

## 3 Image QC
A field that was acquired out of focus cannot be quantified. Such a well fails QC (`qc_pass` = 0).
It stays in the well table but is not used for anything downstream (neither as a control nor as a
point of a curve). Every other well passes (`qc_pass` = 1). Few cells, a weak signal or any other
compound effect is a result, not a QC failure.

## 4 Response scale
Responses are relative to vehicle-treated cells:

    response = readout of the well / readout the same well would have given with vehicle only

so that vehicle-treated cells have response 1 for both readouts. The reference inhibitor (single
concentration) is a plate control; it does not define the scale.

## 5 Curve model, potency, `inactive`
For each compound and readout, the responses of the usable wells are described by a four-parameter
logistic function of the concentration c whose upper asymptote is fixed at the vehicle level:

    response(c) = bottom + (1 - bottom) / (1 + (c / EC50)^h),   0 <= bottom < 1,  h > 0

The reported potency is the concentration at which the fitted curve equals 0.5, i.e. the readout
is half of that of vehicle-treated cells (IC50 for R2, GI50 for R1). This equals EC50 only when
bottom = 0.

A compound is **inactive** on a readout when the fitted curve does not come down to 0.5 at or
below the highest concentration tested in a usable well (no effect, a plateau above 0.5, a
potency beyond the tested range, or a readout that does not fall at all). In that case the
literal `inactive` is reported, never a number.
