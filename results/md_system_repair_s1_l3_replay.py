"""Re-validate md-system-repair-s0001-l3 with the final verifier without repeating the MD runs.

The complete `make` of this task (first block of results/md_system_repair_s1_l3.log) was run before
three verifier constants were tightened. No solution, probe or attack script changed afterwards, so
this script runs the normal `validate.validate` (all gates, Harbor included) but answers every
`_run_script` call with the output directory that the same script produced in that run
(kept under results/work/md-system-repair-s0001-l3.first_verifier/). Images, verifier runs, generic attacks and the Harbor
oracle/nop runs are real and fresh. Usage (from pipeline/):  ../.venv/bin/python ../results/md_system_repair_s1_l3_replay.py
"""
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "pipeline"))
from alegen import harbor, validate  # noqa: E402
from alegen.families import md_system_repair as fam  # noqa: E402

kept = ROOT / "results" / "work" / "md-system-repair-s0001-l3.first_verifier"
b = fam.generate(seed=1, level=3, workdir=ROOT / "results" / "gen" / "md_system_repair-s0001-l3")
task_dir = harbor.write_task(b, ROOT / "tasks")
work = ROOT / "results" / "work" / b.task_id


def replay(agent_tag, script, task_dir, out, timeout, files=None):
    src = kept / Path(out).relative_to(work)
    if not src.is_dir():
        raise RuntimeError(f"no kept output for {src}")
    if Path(out).exists():
        shutil.rmtree(out)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(src, out)
    return subprocess.CompletedProcess([], 0, "", "")


validate._run_script = replay
rep = validate.validate(b, task_dir, ROOT / "results", harbor=True, keep_work=True)
print({"task": b.task_id, "ok": rep["ok"], "minutes": rep["validation_minutes"], "mode": "replayed MD outputs, fresh verifier + Harbor"})
