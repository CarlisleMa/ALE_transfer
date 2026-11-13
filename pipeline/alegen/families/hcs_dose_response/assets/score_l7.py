"""Verifier for level 7 of the high-content-screening dose-response family.

The scoring is that of level 4. `score_l4.py` next to this file is the level-4 verifier, byte for byte (it is
`tests/score.py` of levels 4 and 6); every gate, component, weight and tolerance comes from it.

What level 7 adds is one decision about the names in the two per-well tables. The plate lay on the stage
turned by 180 degrees: the image pair named for a well shows the well at the opposite position of the plate
(reference/truth.json: plants.plate_orientation.image_shows_plate_well). The planted per-image truth (cell
number, focus state, annotation) is stored under the name of the image. An analyst who found the turn can
file a measurement in either of two ways:

  by image     the row carries the name of the image file (and of the annotation file) it was measured on
  by plate     the row carries the name of the plate well that image shows, i.e. the opposite well name

"One row for every well of plate_layout.tsv" and "`field` is the well id" can be read either way once the
turn is known, no text says which is meant, and none may be added at this level. Both namings are therefore
accepted, for `well_measurements.tsv` and for `validation.tsv` independently: the submission is scored by the
level-4 logic under each of the four combinations and the best total is the score. Under every combination a
row of `validation.tsv` is compared with the planted count of the image it refers to and with the `n_cells`
the submission reports for that same image. `dose_response.tsv` has no well names and is scored once.
`info.well_naming` in details.json says which combination was used and gives all four totals.

The potencies decide 0.65 of the score and do not depend on the naming, so no naming can lift a submission
that analysed the plate with the layout the wrong way round.
"""
import sys

sys.path.insert(0, "/tests")
import _scoring  # noqa: E402
from _scoring import CandidateError, Ctx, run  # noqa: E402

# score_l4.py is written as a script: its last line is run(main, pass_threshold=...). Importing it with
# _scoring.run replaced records that call instead of executing it, so its main() can be run per naming.
_l4 = {}


def _record(main, pass_threshold=0.9):
    _l4["main"], _l4["pass_threshold"] = main, pass_threshold


_scoring.run = _record
import score_l4  # noqa: E402
_scoring.run = run
if "main" not in _l4:
    raise RuntimeError("score_l4.py did not hand over its main()")

KEY_COLUMN = {"well_measurements.tsv": "well", "validation.tsv": "field"}


class Naming(Ctx):
    """The submission as the level-4 verifier sees it, with the well names of the tables in `by_plate`
    translated from plate wells to the images that show them. Everything else is read unchanged."""

    def __init__(self, output, reference, by_plate, image_of):
        super().__init__(output, reference)
        self.by_plate, self.image_of = by_plate, image_of

    def table(self, rel, required, *args, **kwargs):
        rows = super().table(rel, required, *args, **kwargs)
        col = KEY_COLUMN.get(rel)
        if col and rel in self.by_plate:
            for r in rows:
                w = score_l4.well_id(r.get(col))
                if w in self.image_of:
                    r[col] = self.image_of[w]
        return rows


def main(ctx):
    shows = ctx.ref_json("truth.json")["plants"]["plate_orientation"]["image_shows_plate_well"]
    image_of = {plate: image for image, plate in shows.items()}
    namings = []
    for name, by_plate in (("both tables by image name", ()),
                           ("both tables by plate well", ("well_measurements.tsv", "validation.tsv")),
                           ("well_measurements.tsv by plate well, validation.tsv by image name", ("well_measurements.tsv",)),
                           ("well_measurements.tsv by image name, validation.tsv by plate well", ("validation.tsv",))):
        r = Naming(ctx.output, ctx.reference, by_plate, image_of)
        try:
            _l4["main"](r)
            total, failure = r.total(), None
        except CandidateError as e:           # a contract violation under this naming
            total, failure = 0.0, str(e)
        namings.append({"name": name, "total": total, "failure": failure, "ctx": r})
    best = max(namings, key=lambda x: (x["failure"] is None, x["total"]))    # equal totals: the first listed
    if best["failure"] is not None:           # no naming under which the tables satisfy the contract
        raise CandidateError(namings[0]["failure"])
    used = best["ctx"]
    ctx.gates, ctx.components = used.gates, used.components
    ctx.info.update(used.info)
    ctx.info["well_naming"] = {"used": best["name"],
                               "totals": {n["name"]: round(n["total"], 6) for n in namings},
                               "contract_violations": {n["name"]: n["failure"] for n in namings if n["failure"]}}


run(main, pass_threshold=_l4["pass_threshold"])
