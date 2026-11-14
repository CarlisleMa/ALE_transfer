"""Level 7 of `hcs_dose_response`: the level-4 plate of the seed plus one silent decision point.

  T_plate_orientation   The plate lay on the microscope stage turned by 180 degrees. The acquisition software
                        names the images by stage position, so the image pair named for well (row r, column c)
                        shows plate well (row 9 - r, column 13 - c): `A01` shows the well the layout calls
                        `H12`. `plate_layout.tsv` is what the wet lab dispensed (plate coordinates) and is the
                        level-4 file; so are the instruction, the SOP and the notes (`python -m alegen make
                        --level 7` then removes the pointer sentences of the notes as for level 6).

What an analysis that joins images and layout by name gets wrong: compound wells serve as vehicle wells, the
reference wells are mid-series wells, every dose series is read backwards and credited to the compound of the
opposite row. Nothing fails: counting, its validation against the manual counts and the focus QC are per
image and are right; curves are fitted and potencies or `inactive` are reported.

Truth by construction, unchanged in meaning. The level-4 generator is run unchanged for the seed (same
cells, compounds, potencies, position effects, out-of-focus fields: all of these belong to the plate and turn
with it). This module then
  * files the two images of every plate well under the name of the opposite well. The pixel data are the
    level-4 files: a field is turned with its plate, but neither the simulated cell layer (placement of
    nuclei and fragments is symmetric under that turn) nor the optics (illumination falls off radially) has
    a direction, so the level-4 rendering of the well is a draw from the same distribution; flat-field and
    dark frames belong to the optical path and are the level-4 files;
  * keeps the level-4 choice of annotated fields: the lab counted the same eight physical fields (the plate
    wells chosen for coverage: dense and sparse untreated fields, fragment-rich fields, fields with enlarged
    nuclei). The acquisition software names images by stage position, so each annotation file carries the
    name of the image that shows its field, i.e. the opposite of the level-4 name. The files themselves are
    the level-4 files, byte for byte: same clicks, same manual counts (maintainer decision of 2026-11-14;
    the first version kept the level-4 names and annotated whatever image lay under them, which changed the
    coverage of the manual counts and with it the difficulty of the counting part);
  * stores the per-image truth (cell number, focus state, annotated or not) under the image name;
  * leaves `dose_response` truth alone: it is per compound.

Scorer: `assets/score_l7.py`, a thin wrapper around the level-4 verifier, which is shipped unchanged next to
it as `tests/score_l4.py`. The rows of the two per-well tables are accepted under the image name or under the
name of the plate well the image shows (see the task card).

The generator refuses a seed for which, on the planted values, the control wells do not tell the two
orientations apart (see `_checks`): the layout of the controls can be symmetric under the turn.

Levels 1-6 do not import this module.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from ...spec import Attack, Probe, Trap
from . import level4, sim, sim4

ASSETS = Path(__file__).parent / "assets"
OFF = 1.3                 # a vehicle well "looks different" beyond this factor (the reference solution's limit)
RHO = 0.7                 # a series "falls" / "rises" with concentration from this rank correlation

L7_JUDGE = """

**Level-7 decision point (the plate lay on the stage the other way round).** The image pair named for a well
shows the well at the opposite position of the plate (row 9 - r, column 13 - c); `plate_layout.tsv` is in
plate coordinates and unchanged. The planted truth is that of the level-4 plate of the same seed with nothing
redefined: potencies are per compound and are the level-4 values; cell number, focus state and manual
annotation belong to an image and are stored under the name of the image. Gates, components, weights and
tolerances are those of level 4: `tests/score_l4.py` is the level-4 verifier, byte for byte, and
`tests/score.py` only runs it. No deliverable asks for the orientation.

**The two per-well tables are accepted under either naming of their rows.** An analyst who found the turn can
file a measurement under the name of the image (and annotation file) it was made on, which is how the truth
is stored, or under the name of the plate well that image shows, which is what a reader who joins the table
with the layout needs. "One row for every well of `plate_layout.tsv`" and "`field` is the well id" can be
read either way once the turn is known, no text says which is meant, and none may be added at this level. The
verifier therefore runs the level-4 logic under the four combinations (each of `well_measurements.tsv` and
`validation.tsv` by image or by plate well) and the best total is the score; under every combination a
validation row is compared with the planted count of the image it refers to and with the `n_cells` reported
for that same image. `info.well_naming` in details.json records the combination used and all four totals.
Potencies (0.65 of the score) carry no well names, so the naming cannot rescue an analysis made with the
layout the wrong way round. Measured: the reference result scores 1.0 by image name and 1.0 when both tables
are filed by plate well; the by-name submissions below score 0.35 under the naming they use and 0.0 under
the three others.

**What the layout by name costs** (measured, seed 1 / seed 2). `l6_reference`, the level-4 reference
solution unchanged: 0.35 / 0.35. It keeps the three per-image components (every well without annotation
within the count tolerance, validation table 8 of 8, focus flags exact) and earns nothing on either potency
component: marker 4 / 3 of 8 compounds right, cell number 5 / 4 of 8, where a component needs a mean credit
above 0.8. On seed 1 it reports `inactive` for all 16 entries, so the entries it has right are the
compounds that are inactive anyway. On seed 2 it reports 12 times `inactive` and four numbers: marker IC50s
of 0.0013 and 0.0027 uM for compounds whose planted IC50s are 5.3 and 3.5 uM, 17 uM for a compound that is
inactive, and a GI50 of 0.0019 uM for a compound that is inactive on cell number; the three small values lie
below the lowest tested concentration (0.0137 uM). The curves behind these tables are not smooth inhibition
curves: a series holds the wells of the opposite row in reverse order, with that row's vehicle and reference
wells in between, so responses rise with the nominal concentration (seed 2, one series: 0.04, 0.08, 0.23,
0.54, 0.81, 0.95, 1.02, 0.90) or jump (seed 1, one series: 0.09, 1.47, 0.75, 0.08, 1.59, 1.53, 1.44, 0.94);
the SOP model, whose upper asymptote is fixed at the vehicle level, turns a rising series into `inactive`
or into a potency at the lower end of its search range. The table alone does not show any of this.
`layout_by_name_everything_else_right` (planted cell numbers, marker levels and focus flags; reference plate
model and fit with the layout by name): 0.35 / 0.35, the same potency tables as `l6_reference` on both
seeds. `controls_swapped_only` (layout by name, the roles of vehicle and reference wells exchanged): 0.35 /
0.35 (seed 1: 14 numeric potencies, marker 1 of 8 right; seed 2: marker 3 of 8, cell number 5 of 8).
`prior_only` shows the floor from the other side: `inactive` everywhere is right for four of eight marker
and five of eight count entries on every plate of this family, which earns no potency credit.

**Routes that reach the truth** (probes): the reference solution finds the orientation from the control wells
(reference wells inhibited, vehicle wells alike), the second solution from the direction of the dose series
without using the controls; both then apply the layout turned and report the per-well tables by image name.
`per_well_tables_by_plate_well` and `well_table_by_plate_well_validation_by_image` are the planted truth
filed the other way.

**The annotated fields are those of level 4.** The lab counted the same eight physical fields (the plate
wells chosen at level 4 to cover dense and sparse untreated fields, fragment-rich fields and fields with
enlarged nuclei). The annotation files are the level-4 files, byte for byte, under the names of the images
that show those fields; `validation.tsv` lists exactly these eight fields (by image name or by plate well,
see above). The counting part of the task is therefore that of levels 4 and 6. Measured on seeds 1, 2, 3 and
10: the reference solution selects the outline threshold it selects on the level-4 inputs (0.12; seed 1 was
compared directly, same worst annotated field at 0.76 of the SOP limit), and `naive_no_focus_qc` scores
0.544 / 0.544 / 0.669 / 0.581, under its level-4 ceiling of 0.75. (In the first version of this level the
annotation files kept their level-4 names and held whatever image lay under them; the tuning then chose
another threshold on some seeds, that attack reached 0.825 and its ceiling had been raised to 0.85. Both
are gone.)

**A side effect that is evidence, not text.** The annotation files now carry names under which the layout
lists something else than what was counted: the entry of the decision point above gives, for the eight
files of this seed, the manual count, what the layout calls the name and which plate well the field is.
Seed 1 / seed 2: 3 / 3 of the 8 files hold fewer than half the cells of an untreated field of their plate
area (about 120 cells in inner wells, 66 in the outer ring), 2 / 1 of them under a name that the layout
calls vehicle or reference (seed 1: `H01` 10 cells, layout vehicle; `E03` 31 cells, layout reference; seed
2: `H01` 8 cells, layout vehicle). An analyst who compares the manual counts with the layout sees it; the
assay notes say how the fields were chosen (untreated wells, highest and intermediate concentrations) and
are the level-4 notes.

**Seeds without a level 7.** A seed whose control wells sit symmetrically under the turn would show
well-behaved controls in both orientations. The generator refuses such seeds (41 of 400 surveyed): it
requires, on the planted values, that with the layout by name fewer than half of the wells called reference
are inhibited and at least four wells called vehicle are off by more than a factor 1.3, and that the
balance of falling over rising dose series separates the two orientations by at least five series. The
check does not use the annotated fields.
"""


def turned(well: str) -> str:
    """The well at the opposite position of the plate (rotation by 180 degrees)."""
    return f"{sim.ROWS[7 - sim.ROWS.index(well[0])]}{sim.NCOL + 1 - int(well[1:]):02d}"


def _spearman(x, y) -> float:
    def rank(v):
        v = np.asarray(v, dtype=float)
        order = np.argsort(v, kind="stable")
        r = np.empty(len(v))
        r[order] = np.arange(len(v), dtype=float)
        for val in np.unique(v):                      # average rank for ties
            m = v == val
            r[m] = r[m].mean()
        return r
    rx, ry = rank(x), rank(y)
    sx, sy = rx.std(), ry.std()
    return float(((rx - rx.mean()) * (ry - ry.mean())).mean() / (sx * sy)) if sx > 0 and sy > 0 else 0.0


def _reading(lay: list, content: list) -> dict:
    """Control and series statistics of one way to lay the layout over the images, on planted values.
    lay[i] = layout entry (role, compound, dose, row, col); content[i] = the planted well whose image that
    entry is applied to (n_cells, marker_mean, in focus, planted responses free of position effects)."""
    use = [i for i in range(len(lay)) if content[i]["focus"] and content[i]["n"] > 0]
    veh = [i for i in use if lay[i]["role"] == "vehicle"]
    ref = [i for i in use if lay[i]["role"] == "reference"]
    # reference wells against the vehicle wells of their own column (every column holds two)
    ratios = []
    for i in ref:
        col = [content[v]["m"] for v in veh if lay[v]["col"] == lay[i]["col"]]
        if col:
            ratios.append(content[i]["m"] / float(np.mean(col)))
    # the two vehicle wells of a column against each other (marker per cell)
    pair = []
    for c in range(sim.NCOL):
        col = [content[v]["m"] for v in veh if lay[v]["col"] == c]
        if len(col) == 2:
            pair.append(max(col) / min(col))
    # cell number of the vehicle wells against the median of their plate area (outer ring / inner wells)
    rel_n = []
    for area in (True, False):
        grp = [content[v]["n"] for v in veh if sim.is_edge(lay[v]["row"], lay[v]["col"]) == area]
        rel_n += [n / float(np.median(grp)) for n in grp]
    off_n = sum(1 for v in rel_n if v < 1 / OFF or v > OFF)
    # vehicle wells that an analyst would call different: planted response beyond the factor on either readout
    off = sum(1 for v in veh if not (1 / OFF <= content[v]["resp_m"] <= OFF and 1 / OFF <= content[v]["resp_n"] <= OFF))
    # direction of the dose series
    fall = {"marker": 0, "count": 0}
    rise = {"marker": 0, "count": 0}
    for name in sorted({e["compound"] for e in lay if e["role"] == "compound"}):
        s = [i for i in use if lay[i]["compound"] == name]
        if len(s) < 6:
            continue
        d = [lay[i]["dose"] for i in s]
        for key, k in (("marker", "resp_m"), ("count", "resp_n")):
            rho = _spearman(d, [content[i][k] for i in s])
            fall[key] += rho <= -RHO
            rise[key] += rho >= RHO
    return {"reference_wells": len(ratios),
            "reference_over_vehicle_marker": [round(float(np.min(ratios)), 3), round(float(np.median(ratios)), 3),
                                              round(float(np.max(ratios)), 3)],
            "reference_wells_below_half": int(sum(r < 0.5 for r in ratios)),
            "vehicle_wells": len(veh), "vehicle_wells_off": int(off),
            "vehicle_pair_marker_ratio": [round(float(np.median(pair)), 2), round(float(np.max(pair)), 2)],
            "vehicle_pairs_beyond_factor": int(sum(p > OFF for p in pair)), "vehicle_pairs": len(pair),
            "vehicle_cells_vs_area_median": [round(float(min(rel_n)), 2), round(float(max(rel_n)), 2)],
            "vehicle_wells_cell_number_off": int(off_n),
            "series_falling": fall, "series_rising": rise}


def _balance(r: dict) -> int:
    """Dose series that fall with the concentration minus series that rise (both readouts)."""
    return sum(r["series_falling"].values()) - sum(r["series_rising"].values())


def _checks(ev: dict) -> dict:
    """What a seed must satisfy, on the planted values, to get a level 7."""
    t, n = ev["turned"], ev["by_name"]
    return {
        # the controls behave in the true orientation ...
        "turned_reference_wells_inhibited": t["reference_wells_below_half"] == t["reference_wells"]
        and t["reference_over_vehicle_marker"][2] < 0.2,
        "turned_vehicle_wells_alike": t["vehicle_wells_off"] == 0 and t["vehicle_pairs_beyond_factor"] == 0,
        # ... and do not when the layout is joined by name: the wells called reference are not inhibited as a
        # group (the rule of the reference solution), and several wells called vehicle differ from the rest
        "by_name_reference_wells_not_inhibited": 2 * n["reference_wells_below_half"] < n["reference_wells"]
        and n["reference_over_vehicle_marker"][1] >= 0.5,
        "by_name_vehicle_wells_differ": n["vehicle_wells_off"] >= 4,
        # the dose series tell the same story (the statistic of the second solution): series that fall with the
        # concentration minus series that rise, both readouts
        "series_tell_the_orientation": _balance(t) >= 5 and _balance(t) - _balance(n) >= 5,
    }


def evidence(wells: list) -> dict:
    """Both readings on the planted values of a plate (wells in layout order, A01..H12)."""
    lay = [{"role": w.role, "compound": w.compound, "dose": w.dose, "row": w.row, "col": w.col} for w in wells]
    plate = [{"n": w.n_cells, "m": w.marker_mean if w.marker_mean > 0 else w.extra["marker_level"],
              "focus": w.defocus == 0,
              "resp_m": (w.marker_mean if w.marker_mean > 0 else w.extra["marker_level"]) / (sim.M0 * w.stain_factor),
              "resp_n": w.n_cells / (sim.N0 * w.count_factor)} for w in wells]
    n = len(wells)
    ev = {"turned": _reading(lay, plate), "by_name": _reading(lay, [plate[n - 1 - i] for i in range(n)])}
    ev["checks"] = _checks(ev)
    return ev


def survey(seed: int) -> dict:
    """The checks of a seed without rendering the plate (planned cell numbers and marker levels)."""
    rng = np.random.default_rng(884000 + seed)
    cs = sim.draw_compounds(rng, 3)
    sim.assign_rows(rng, cs, 3)
    shapes = sim4.reshape_compounds(rng, cs)
    wells, _ = sim.layout(rng, cs)
    sim4.plant_wells4(rng, wells, cs, shapes)
    sim.plant_defocus(rng, wells, cs, 3)
    return evidence(wells)


def _rows(tsv: str) -> tuple[list, list]:
    lines = [ln.split("\t") for ln in tsv.splitlines() if ln.strip()]
    return lines[0], lines[1:]


def generate_l7(seed: int, level: int, workdir: Path, notes_fn):
    if level != 7:
        raise ValueError("l7.generate_l7 builds level 7 only")
    pl: dict = {}
    b = level4.generate_l4(seed, level, workdir, notes_fn, plate_out=pl)     # the level-4 plate of the seed
    wells = pl["wells"]
    names = [w.well for w in wells]
    shows = {n: turned(n) for n in names}                  # image name -> plate well it shows
    by_plate = {w.well: w for w in wells}
    if sorted(shows.values()) != sorted(names) or any(turned(shows[n]) != n for n in names):
        raise RuntimeError("level 7: the turn is not a one-to-one map of the plate")

    # ---- does the plate tell the two orientations apart? (planted values) -------------------------------
    ev = evidence(wells)
    failed = [k for k, ok in ev["checks"].items() if not ok]
    if failed:
        raise RuntimeError(f"level 7, seed {seed}: the planted plate does not single out its orientation "
                           f"({failed}); by name {ev['by_name']}, turned {ev['turned']}")

    # ---- images: every pair filed under the name of the opposite well -----------------------------------
    level4_file = {k: v for k, v in b.input_files.items() if k.startswith("images/")}
    for n in names:
        for ch in (1, 2):
            src = level4_file[f"images/{shows[n]}_ch{ch}.tif"]
            if not isinstance(src, Path):
                raise RuntimeError("level 7: unexpected image source of the level-4 build")
            b.input_files[f"images/{n}_ch{ch}.tif"] = src

    # ---- annotations: the fields the lab counted at level 4, under the names of the images that show them --
    in_focus = {n: by_plate[shows[n]].defocus == 0 for n in names}
    annotated4 = sorted(pl["annotated"])                   # plate wells, chosen for coverage (sim4.choose_annotated)
    level4_ann = {k: b.input_files.pop(k) for k in [k for k in b.input_files if k.startswith("annotations/")]}
    if sorted(level4_ann) != [f"annotations/{p}.tsv" for p in annotated4]:
        raise RuntimeError("level 7: unexpected annotation files of the level-4 build")
    if any(by_plate[p].defocus > 0 for p in annotated4):   # level 4 never annotates an out-of-focus field
        raise RuntimeError("level 7: the level-4 build annotated an out-of-focus field")
    annotated = sorted(turned(p) for p in annotated4)      # image names
    for n in annotated:                                    # the level-4 file: same clicks, same manual count
        b.input_files[f"annotations/{n}.tsv"] = level4_ann[f"annotations/{shows[n]}.tsv"]
    if len(annotated) != 8 or len(set(annotated)) != 8 or not all(in_focus[n] for n in annotated):
        raise RuntimeError("level 7: eight in-focus annotated fields expected")
    for n in annotated:
        if len(pl["points"][shows[n]]) != by_plate[shows[n]].n_cells:
            raise RuntimeError("level 7: manual count and planted count of an annotated field differ")

    # ---- planted truth: per image under the image name, per compound unchanged --------------------------
    truth = json.loads(b.tests_files["reference/truth.json"])
    plate_wells = truth["wells"]
    truth["wells"] = {}
    for n in names:
        t = dict(plate_wells[shows[n]])                    # what the image holds: the plate well it shows
        t["annotated"] = int(n in annotated)
        t["shows_plate_well"] = shows[n]
        t["layout_entry_of_this_name"] = {k: plate_wells[n][k] for k in ("role", "compound", "dose_uM")}
        truth["wells"][n] = t
        if t["n_cells"] != by_plate[shows[n]].n_cells or t["qc_fail"] != int(not in_focus[n]):
            raise RuntimeError("level 7: per-image truth and planted plate disagree")
    truth["plate"]["out_of_focus_plate_wells"] = truth["plate"]["out_of_focus"]
    truth["plate"]["out_of_focus"] = sorted(turned(w) for w in truth["plate"]["out_of_focus"])     # image names
    truth["plate"]["annotated_fields"] = annotated
    # what the layout says about the names the annotation files now carry, next to what they hold
    untreated = {area: float(np.median([w.n_cells for w in wells if w.role == "vehicle" and w.defocus == 0
                                        and sim.is_edge(w.row, w.col) == area])) for area in (True, False)}
    named = []
    for n in annotated:
        lay, w = by_plate[n], by_plate[shows[n]]           # layout entry of the name / the field that was counted
        named.append({"file": n, "layout_role": lay.role, "layout_compound": lay.compound, "layout_dose_uM": lay.dose,
                      "manual_count": w.n_cells, "shows_plate_well": w.well, "field_role": w.role,
                      "field_compound": w.compound, "field_dose_uM": w.dose,
                      "untreated_cells_of_plate_area": untreated[sim.is_edge(w.row, w.col)],
                      "count_over_untreated": round(w.n_cells / untreated[sim.is_edge(w.row, w.col)], 3)})
    truth["plants"] = {"plate_orientation": {
        "turned_by_degrees": 180, "image_shows_plate_well": shows, "annotated_plate_wells": annotated4,
        "annotated_fields_as_the_layout_names_them": named, "evidence_on_planted_values": ev}}
    b.tests_files["reference/truth.json"] = json.dumps(truth, indent=1)
    # the level-4 verifier, unchanged, as a module; score.py runs it under each naming of the per-well rows
    b.tests_files["score_l4.py"] = b.tests_files["score.py"]
    b.tests_files["score.py"] = ASSETS / "score_l7.py"

    # ---- solutions ----------------------------------------------------------------------------------------
    b.solution_files = {"solve.sh": ASSETS / "solve_l7.sh", "hcs_l7_voting.py": ASSETS / "hcs_l7_voting.py",
                        "hcs_l7_cellpose.py": ASSETS / "hcs_l7_cellpose.py"}

    # ---- battery: the level-4 battery on the level-7 data, plus the level-7 probes and attacks ----------
    def refiled(files: dict, wm_by: str = "image", val_by: str = "image") -> dict:
        """A level-4 submission written from planted values, for the level-7 data: every per-well row moves
        with its image. `*_by` = "plate" files the rows under the plate well instead (the other naming)."""
        head, rows = _rows(files["well_measurements.tsv"])
        of_plate = {r[0]: r[1:] for r in rows}             # level 4: plate well = image name
        key = {"image": lambda n: n, "plate": lambda n: shows[n]}
        wm = ["\t".join(head)] + ["\t".join([key[wm_by](n)] + of_plate[shows[n]]) for n in names]
        val = ["field\tmanual_count\tpipeline_count"] + [
            f"{key[val_by](n)}\t{by_plate[shows[n]].n_cells}\t{of_plate[shows[n]][0]}" for n in annotated]
        return dict(files, **{"well_measurements.tsv": "\n".join(wm) + "\n", "validation.tsv": "\n".join(val) + "\n"})

    run4, run7 = "python3 /alegen_script/hcs_l4_voting.py ", "python3 /alegen_script/hcs_l7_voting.py "
    l4_py = {"hcs_l4_voting.py": ASSETS / "hcs_l4_voting.py"}
    l7_py = {"hcs_l7_voting.py": ASSETS / "hcs_l7_voting.py"}
    ceil = b.pass_threshold - 0.15

    probes = []
    for p in b.probes:
        if p.name == "independent_solution":
            probes.append(Probe(p.name, "second implementation, no code shared with the reference solution: "
                                "Cellpose nuclei model with settings chosen on the annotated fields, as at level 4; "
                                "the orientation of the plate is found from the direction of the dose series "
                                "(rank correlation of both readouts with the concentration, readouts normalised "
                                "without the control wells), not from the controls", p.min_score, script=p.script))
        elif p.script is None:
            probes.append(Probe(p.name, p.description + " (per-well rows under the image names)", p.min_score,
                                files=refiled(p.files)))
        else:
            raise RuntimeError(f"level 7: unexpected script probe {p.name} in the level-4 battery")
    truth_files = next(p.files for p in probes if p.name == "planted_truth")
    head, rows = _rows(truth_files["well_measurements.tsv"])
    for r in rows:                                          # the refiled probe is the planted truth of this level
        t = truth["wells"][r[0]]
        if int(r[1]) != t["n_cells"] or int(r[3]) != 1 - t["qc_fail"]:
            raise RuntimeError("level 7: the planted-truth probe does not match the planted per-image truth")
    level4_truth = next(p.files for p in b.probes if p.name == "planted_truth")
    probes += [
        Probe("per_well_tables_by_plate_well", "the turn is recognised and both per-well tables are filed by the "
              "plate well each image shows: planted truth with every row under the opposite well name. Accepted "
              "like the tables by image name: the verifier scores each table under both namings", 0.999,
              files=refiled(level4_truth, "plate", "plate")),
        Probe("well_table_by_plate_well_validation_by_image", "the well table is filed by plate well, the "
              "validation table keeps the names of the annotation files: planted truth. The tie between the two "
              "tables (pipeline_count = n_cells of the same image) holds under that naming", 0.999,
              files=refiled(level4_truth, "plate", "image")),
    ]

    attacks = [
        Attack("l6_reference", "the reference solution of levels 4-6, unchanged, on the level-7 inputs: every "
               "decision point of the lower levels is handled, images and layout are joined by well name", ceil,
               script=run4.strip(), script_files=l4_py),
        Attack("layout_by_name_everything_else_right", "T_plate_orientation ignored and nothing else wrong: planted "
               "cell numbers, planted marker levels and planted focus flags of every image, validation table from "
               "the planted counts, potencies as the reference plate model and curve fit give them when the layout "
               "is joined by well name", ceil,
               script=run7 + "--orientation name --measurements /alegen_script/planted_wells.tsv",
               script_files=dict(l7_py, **{"planted_wells.tsv": truth_files["well_measurements.tsv"]})),
        Attack("controls_swapped_only", "the odd controls are noticed but the layout is not turned: the reference "
               "solution with the layout joined by well name and the roles of the vehicle and the reference wells "
               "exchanged (the wells called reference look untreated and are used as vehicle wells)", ceil,
               script=run7 + "--orientation name --controls swapped", script_files=l7_py),
    ]
    for a in b.attacks:
        if a.script is None:                                # truth-derived files follow the images
            attacks.append(Attack(a.name, a.description, a.max_score, files=refiled(a.files)))
        elif a.script.startswith(run4):                     # "ignores exactly one element": the level-7 step is done
            attacks.append(Attack(a.name, a.description + " (plate controls checked, layout applied turned)",
                                  a.max_score, script=run7 + a.script[len(run4):], script_files=l7_py))
        else:                                               # lower-level methods run as they were written
            attacks.append(Attack(a.name, a.description + " (level-7 inputs, layout joined by well name)",
                                  a.max_score, script=a.script, script_files=a.script_files))
    b.probes, b.attacks = probes, attacks

    # ---- task card ----------------------------------------------------------------------------------------
    t, n = ev["turned"], ev["by_name"]

    def controls(r):
        ro, pr, ce = r["reference_over_vehicle_marker"], r["vehicle_pair_marker_ratio"], r["vehicle_cells_vs_area_median"]
        return (f"marker per cell of the {r['reference_wells']} in-focus wells called reference relative to the "
                f"vehicle wells of their column {ro[0]:.2f}-{ro[2]:.2f} (median {ro[1]:.2f}; "
                f"{r['reference_wells_below_half']} of {r['reference_wells']} below half); the two wells called "
                f"vehicle in a column differ in marker per cell by a factor of {pr[0]:.2f} (median) and at most "
                f"{pr[1]:.2f} ({r['vehicle_pairs_beyond_factor']} of {r['vehicle_pairs']} columns beyond {OFF}); cell "
                f"number of the wells called vehicle {ce[0]:.2f}-{ce[1]:.2f} of the median of their plate area "
                f"(inner wells / outer ring; {r['vehicle_wells_cell_number_off']} of {r['vehicle_wells']} beyond a "
                f"factor {OFF}); wells called vehicle whose planted response is beyond a factor {OFF} on either "
                f"readout {r['vehicle_wells_off']} of {r['vehicle_wells']}; dose series whose response falls / rises "
                f"with the concentration (rank correlation beyond {RHO}; eight series per readout): cell number "
                f"{r['series_falling']['count']} / {r['series_rising']['count']}, marker per cell "
                f"{r['series_falling']['marker']} / {r['series_rising']['marker']}")

    def entry(role, compound, dose):
        return role if role in ("vehicle", "reference") else f"{compound} at {dose:.3g} uM"

    kinds = [f"{x['file']}.tsv, {x['manual_count']} cells (the layout calls {x['file']} "
             f"{entry(x['layout_role'], x['layout_compound'], x['layout_dose_uM'])}; the field is plate well "
             f"{x['shows_plate_well']}, {entry(x['field_role'], x['field_compound'], x['field_dose_uM'])})" for x in named]
    low = [x for x in named if x["count_over_untreated"] < 0.5]
    low_ctrl = [x for x in low if x["layout_role"] in ("vehicle", "reference")]
    b.traps = list(b.traps) + [Trap(
        "T_plate_orientation",
        "The plate lay on the microscope stage turned by 180 degrees. The acquisition software names the images by "
        "stage position, so the image pair named for well (row r, column c) shows plate well (row 9 - r, column "
        "13 - c): A01 shows the well the layout calls H12. plate_layout.tsv is what the wet lab dispensed and is "
        "unchanged; no text mentions orientation, the stage holder, checking the controls or a mismatch. An "
        "analysis that joins images and layout by well name uses compound wells as vehicle wells and mid-series "
        "wells as reference wells, reads every dose series backwards and credits it to the compound of the "
        "opposite row. Nothing fails: counting, its validation against the manual counts and the focus QC are per "
        "image and are right, curves are fitted and potencies or `inactive` are reported. The lab counted the same "
        "eight physical fields as at level 4; each annotation file is the level-4 file under the name of the image "
        "that shows its field: " + "; ".join(kinds) + ".",
        "looking at the plate controls before they are used: a plate map of both readouts with the roles of the "
        "layout, or the control statistics (reference against vehicle wells, spread of the vehicle wells). A plate "
        "fits into a stage holder in exactly two orientations, and the same statistics with the layout turned by "
        "180 degrees settle it. Measured on the planted values (cell number and mean marker level of every image; "
        "out-of-focus images left out). Layout joined by well name: " + controls(n) + ". Layout turned by 180 "
        "degrees: " + controls(t) + ". The reference wells are the plainest sign: the reference inhibitor is "
        "described in the assay notes as an inhibitor of the marker pathway at 10 uM. A second sign, also in the "
        "data and in no text: the manual counts. The annotation files carry names that mean something else in the "
        f"layout than the fields they hold; {len(low)} of the 8 files hold fewer than half the cells of an "
        f"untreated field of their plate area ({untreated[False]:.0f} cells inner wells, {untreated[True]:.0f} outer "
        "ring), " + (f"{len(low_ctrl)} of them under a name the layout calls vehicle or reference"
                     if low_ctrl else "none of them under a name the layout calls vehicle or reference")
        + " (" + "; ".join(f"{x['file']}: {x['manual_count']} cells, layout "
                           f"{entry(x['layout_role'], x['layout_compound'], x['layout_dose_uM'])}" for x in low) + ").",
        "Verification / skipped QC step")]
    b.judge_explanation_md = b.judge_explanation_md.rstrip() + L7_JUDGE
    b.capabilities = list(b.capabilities) + ["checking the plate controls against the layout before using them"]
    b.oracle_steps = 2
    return b
