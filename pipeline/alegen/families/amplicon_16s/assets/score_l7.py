"""Verifier for level 7 of the 16S V4 amplicon family.

The scoring is that of levels 1-6. `score_l3.py` next to this file is the level-3 verifier, byte for byte
(it is `tests/score.py` of the lower levels); every gate, component, weight and tolerance comes from it.

What level 7 adds is one decision about the count table. The reads of one extraction blank and of one study
specimen are filed under each other's `sample_id` (reference/truth.json: plants.control_swap). An analyst
who found that can deliver the two columns of `asv_table.tsv` in either of two consistent ways:

  as filed    each column holds the read pairs of the FASTQ pair filed under that sample_id
              (this is how the planted counts are stored)
  by sample   each column holds the read pairs of the library that sample really is, i.e. the same
              two columns the other way round

Once the exchange is known, "the read pairs of that sample" (instruction) and "of that library" (SOP item 2)
can be read either way, and no text says which is meant. Both are therefore accepted: the submission is
scored twice by the level-3 logic, once as submitted and once with the two columns exchanged (both
together, nothing else touched), and the better total is the score. `info.count_table_reading` in
details.json says which reading was used and gives both totals.
"""
import sys

sys.path.insert(0, "/tests")
import _scoring  # noqa: E402
from _scoring import CandidateError, Ctx, run  # noqa: E402

# score_l3.py is written as a script: its last line is run(main, pass_threshold=...). Importing it with
# _scoring.run replaced records that call instead of executing it, so its main() can be run per reading.
_l3 = {}


def _record(main, pass_threshold=0.9):
    _l3["main"], _l3["pass_threshold"] = main, pass_threshold


_scoring.run = _record
import score_l3  # noqa: E402,F401
_scoring.run = run
if "main" not in _l3:
    raise RuntimeError("score_l3.py did not hand over its main()")


class Reading(Ctx):
    """The submission as the level-3 verifier sees it, optionally with the values of the two exchanged
    sample columns of asv_table.tsv swapped. Everything else is read from the output directory unchanged."""

    def __init__(self, output, reference, swap=None):
        super().__init__(output, reference)
        self.swap = swap

    def table(self, rel, required, *args, **kwargs):
        rows = super().table(rel, required, *args, **kwargs)
        if self.swap and rel == "asv_table.tsv":
            a, b = self.swap
            for r in rows:
                if a in r and b in r:
                    r[a], r[b] = r[b], r[a]
        return rows


def main(ctx):
    plant = ctx.ref_json("truth.json")["plants"]["control_swap"]
    pair = (plant["labelled_blank"], plant["labelled_specimen"])
    readings = []
    for name, swap in (("as_submitted", None), ("two_columns_exchanged", pair)):
        r = Reading(ctx.output, ctx.reference, swap)
        try:
            _l3["main"](r)
            total, failure = r.total(), None
        except CandidateError as e:           # a contract violation is one in either reading
            total, failure = 0.0, str(e)
        readings.append({"name": name, "total": total, "failure": failure, "ctx": r})
    best = max(readings, key=lambda x: x["total"])          # equal totals: the table as submitted
    if best["failure"] is not None:
        raise CandidateError(best["failure"])
    used = best["ctx"]
    ctx.gates, ctx.components = used.gates, used.components
    ctx.info.update(used.info)
    ctx.info["count_table_reading"] = {
        "exchanged_sample_ids": list(pair), "used": best["name"],
        "total_as_submitted": round(readings[0]["total"], 6),
        "total_two_columns_exchanged": round(readings[1]["total"], 6)}


run(main, pass_threshold=_l3["pass_threshold"])
