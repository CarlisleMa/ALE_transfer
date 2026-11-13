"""Score the saved output of finished Harbor trials again with the verifier of another task.

Used for the level-4 / level-5 comparison: a level-5 task has the inputs and truth of its level-4 source but
no audit deliverable, so its scorer has one component less. To compare like with like, the outputs that the
models produced on level 4 are scored again with the level-5 scorer (the extra audit.tsv is ignored).

  python -m alegen.rescore <task_id of the verifier to use> <label> <trial_dir> [<trial_dir> ...]

Each trial_dir is results/jobs/<job>/<trial>; its saved output is <trial_dir>/artifacts/task/output.
Appends to results/rescored/<label>.json:  {"<job>/<trial>": {task, verifier, original, rescored, lowest}}.
"""
from __future__ import annotations

import json
import shutil
import sys
import tempfile
from pathlib import Path

from . import validate

ROOT = Path(__file__).resolve().parents[2]


def main(argv: list[str]) -> int:
    if len(argv) < 3:
        print(__doc__)
        return 2
    verifier_task, label, trials = argv[0], argv[1], [Path(a).resolve() for a in argv[2:]]
    tag = validate.tags(verifier_task)[1]
    out_f = ROOT / "results" / "rescored" / f"{label}.json"
    out_f.parent.mkdir(parents=True, exist_ok=True)
    res = json.loads(out_f.read_text()) if out_f.exists() else {}
    for t in trials:
        out = t / "artifacts" / "task" / "output"
        rj = t / "result.json"
        if not rj.exists():
            continue
        j = json.loads(rj.read_text())
        rw = (j.get("verifier_result") or {}).get("rewards")
        original = rw.get("reward") if isinstance(rw, dict) else rw
        work = Path(tempfile.mkdtemp(prefix="alegen-rescore-"))
        try:
            score, details = (None, {"failure": "no saved output"}) if not out.exists() else validate.verify(tag, out, work)
        finally:
            shutil.rmtree(work, ignore_errors=True)
        comps = details.get("components") or []
        low = min(comps, key=lambda c: c["value"]) if comps else None
        key = f"{t.parent.name}/{t.name}"
        res[key] = {"task": str(j.get("task_name", "")).split("/")[-1], "verifier": verifier_task, "original": original,
                    "rescored": score,
                    "lowest": None if not low or low["value"] == 1.0 else f"{low['name']} = {low['value']:.2f}",
                    "failure": details.get("failure") or details.get("verifier_error")}
        print(key, "original", original, "rescored", score, res[key]["lowest"] or "", res[key]["failure"] or "", flush=True)
    out_f.write_text(json.dumps(res, indent=1, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
