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
  * keeps the level-4 choice of annotated well names (the lab picks fields by the layout it believes in) and
    writes into each annotation file the cells of the image that now carries the name. A name whose image is
    out of focus cannot be counted by hand and is replaced by the next field of the same kind;
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

# The one ceiling that is not pass - 0.15 at level 7. Measured 2026-11-12 (docs/families/hcs_dose_response.md).
NO_FOCUS_QC_NOTE = (
    ". Ceiling 0.85 instead of 0.75 at level 7. What this submission scores depends on which outline threshold the "
    "reference solution selects on the annotated fields, and those are other images than at level 4. With 0.15 "
    "(seeds 1 and 10) it finds no nucleus in an out-of-focus field: the well has no marker value, drops out of "
    "the marker fit by itself and enters the count fit with zero cells; the missing QC then costs the QC flags and "
    "the GI50 of one compound: 0.825, which is 0.075 under the pass line. With 0.12 (seed 2) it finds objects "
    "with a wrong marker value and both potencies of the affected compounds are pulled down: 0.544. The level-4 "
    "task has the same property: on the level-4 inputs of seed 1 both thresholds are equally good on the "
    "annotated fields, the reference solution takes the first (0.12, 0.544), and with 0.15 the same submission "
    "scores 0.825")

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
0.54, 0.81, 0.95, 1.02, 0.90) or jump (seed 1, one series: 0.09, 1.48, 0.76, 0.08, 1.59, 1.53, 1.44, 0.93);
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

**One ceiling of the level-4 battery is raised.** `naive_no_focus_qc` (reference solution with the focus QC
switched off) has a ceiling of 0.85 instead of 0.75; its entry in the attack list gives the measurements. It
stays under the pass line (0.825 at most), but not by the usual margin. The cause is in the level-4 weights
(QC flags 0.05, count potency 0.20, one wrong compound of eight costs 0.125) and not in the level-7 decision
point: when the counting step finds nothing in an out-of-focus field, the marker curve is not affected.

**What is weaker than at level 4.** (1) The annotated fields were chosen at level 4 to cover dense and sparse
untreated fields, fragment-rich fields and fields with enlarged nuclei. The file names are the same, the
images under them are not, so which phenotypes the eight manual counts cover depends on the seed (listed in
the decision point above). On seeds 1 and 2 both solutions met the SOP criterion on every annotated field
and the count tolerance on the other wells; this is measured for two seeds, not guaranteed by construction.
(2) A seed whose control wells sit symmetrically under the turn would show well-behaved controls in both
orientations. The generator refuses such seeds (about one in ten of 400 surveyed): it requires, on the
planted values, that with the layout by name fewer than half of the wells called reference are inhibited
and at least four wells called vehicle are off by more than a factor 1.3, and that the balance of falling
over rising dose series separates the two orientations by at least five series.
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
    blurred = sim.plant_defocus(rng, wells, cs, 3)
    annotated = sim4.choose_annotated(rng, wells, cs)
    ev = evidence(wells)
    ev["annotated_names_with_out_of_focus_image"] = sorted(set(annotated) & {turned(w) for w in blurred})
    return ev


def _annotated_names(annotated4: list, wells: list, in_focus: dict, rng: np.random.Generator):
    """The level-4 choice of names; a name whose image is out of focus is replaced by the next field of the
    same kind (another vehicle well of the same plate area; the neighbouring concentration of the series)."""
    by_name = {w.well: w for w in wells}
    names, subs = list(annotated4), []
    for old in annotated4:
        if in_focus[old]:
            continue
        w = by_name[old]
        if w.role == "vehicle":
            pool = [x.well for x in wells if x.role == "vehicle" and in_focus[x.well] and x.well not in names
                    and sim.is_edge(x.row, x.col) == sim.is_edge(w.row, w.col)]
            new = pool[int(rng.integers(0, len(pool)))]
        else:
            series = sorted((x for x in wells if x.compound == w.compound), key=lambda x: x.dose)
            k = series.index(w)
            near = sorted((x for x in series if x is not w), key=lambda x: (abs(series.index(x) - k), -x.dose))
            new = next(x.well for x in near if in_focus[x.well] and x.well not in names)
        names[names.index(old)] = new
        subs.append({"level4_name": old, "annotated_instead": new,
                     "reason": "the image filed under the level-4 name is out of focus"})
    return sorted(names), subs


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

    # ---- annotations: the level-4 names, the cells of the images that now carry them --------------------
    drng = np.random.default_rng(887000 + seed)            # own stream: the level-4 plate of the seed is unchanged
    in_focus = {n: by_plate[shows[n]].defocus == 0 for n in names}
    annotated, substituted = _annotated_names(sorted(pl["annotated"]), wells, in_focus, drng)
    for k in [k for k in b.input_files if k.startswith("annotations/")]:
        del b.input_files[k]
    adir = workdir / "annotations_l7"
    adir.mkdir(parents=True, exist_ok=True)
    for f in adir.glob("*.tsv"):
        f.unlink()
    for n in annotated:
        p = adir / f"{n}.tsv"
        p.write_text(sim4.annotation_tsv(drng, pl["points"][shows[n]]))
        b.input_files[f"annotations/{n}.tsv"] = p
    if len(annotated) != 8 or len(set(annotated)) != 8 or not all(in_focus[n] for n in annotated):
        raise RuntimeError("level 7: eight in-focus annotated fields expected")

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
    truth["plants"] = {"plate_orientation": {
        "turned_by_degrees": 180, "image_shows_plate_well": shows, "annotated_names_level4": sorted(pl["annotated"]),
        "annotated_names_replaced": substituted, "evidence_on_planted_values": ev}}
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
            note, cap = " (plate controls checked, layout applied turned)", a.max_score
            if a.name == "naive_no_focus_qc":               # measured exception, see NO_FOCUS_QC_NOTE
                note, cap = note + NO_FOCUS_QC_NOTE, max(cap, b.pass_threshold - 0.05)
            attacks.append(Attack(a.name, a.description + note, cap, script=run7 + a.script[len(run4):],
                                  script_files=l7_py))
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

    kinds = []
    for a in annotated:
        w = by_plate[shows[a]]
        what = w.role if w.role in ("vehicle", "reference") else f"{w.compound} at {w.dose:.3g} uM"
        kinds.append(f"{a} shows {what} ({w.n_cells} cells)")
    b.traps = list(b.traps) + [Trap(
        "T_plate_orientation",
        "The plate lay on the microscope stage turned by 180 degrees. The acquisition software names the images by "
        "stage position, so the image pair named for well (row r, column c) shows plate well (row 9 - r, column "
        "13 - c): A01 shows the well the layout calls H12. plate_layout.tsv is what the wet lab dispensed and is "
        "unchanged; no text mentions orientation, the stage holder, checking the controls or a mismatch. An "
        "analysis that joins images and layout by well name uses compound wells as vehicle wells and mid-series "
        "wells as reference wells, reads every dose series backwards and credits it to the compound of the "
        "opposite row. Nothing fails: counting, its validation against the manual counts and the focus QC are per "
        "image and are right (the annotation files carry the level-4 names and hold the cells of the images that "
        "now carry those names: " + "; ".join(kinds) + (
            "; replaced because the image under the level-4 name is out of focus: " + ", ".join(
                f"{s['level4_name']} -> {s['annotated_instead']}" for s in substituted) if substituted else "")
        + "), curves are fitted and potencies or `inactive` are reported.",
        "looking at the plate controls before they are used: a plate map of both readouts with the roles of the "
        "layout, or the control statistics (reference against vehicle wells, spread of the vehicle wells). A plate "
        "fits into a stage holder in exactly two orientations, and the same statistics with the layout turned by "
        "180 degrees settle it. Measured on the planted values (cell number and mean marker level of every image; "
        "out-of-focus images left out). Layout joined by well name: " + controls(n) + ". Layout turned by 180 "
        "degrees: " + controls(t) + ". The reference wells are the plainest sign: the reference inhibitor is "
        "described in the assay notes as an inhibitor of the marker pathway at 10 uM.",
        "Verification / skipped QC step")]
    b.judge_explanation_md = b.judge_explanation_md.rstrip() + L7_JUDGE
    b.capabilities = list(b.capabilities) + ["checking the plate controls against the layout before using them"]
    b.oracle_steps = 2
    return b
