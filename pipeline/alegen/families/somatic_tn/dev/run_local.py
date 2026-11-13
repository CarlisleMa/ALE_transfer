"""Development harness (not shipped in tasks): run a script of the family in the tools image against a
generated task directory and score the result on the host.

  ../.venv/bin/python -m alegen.families.somatic_tn.dev.run_local <task_dir> <name> "<script>"
"""
import json
import subprocess
import sys
import time
from pathlib import Path

from alegen import dockerlib
from alegen.families import somatic_tn

ASSETS = Path(somatic_tn.__file__).parent / "assets"


def main():
    task, name, script = Path(sys.argv[1]).resolve(), sys.argv[2], sys.argv[3]
    import os
    scratch = Path(os.environ.get("SOMATIC_TN_SCRATCH", "/tmp/somatic_tn_dev")) / task.name
    out = scratch / name / "output"
    tag = dockerlib.tool_image(conda=somatic_tn.CONDA)
    t = time.time()
    r = dockerlib.run_in_image(tag, f"mkdir -p /task/output && cd /task && {script}",
                               mounts={task / "environment" / "input": "/task/input", ASSETS: "/solution",
                                       ASSETS.parent / "dev": "/dev_scripts"},
                               copy_out={"/task/output": out}, cpus=2)
    mins = (time.time() - t) / 60
    if r.returncode != 0:
        print("SCRIPT FAILED", r.stdout[-1500:], r.stderr[-3000:])
    logs = scratch / name / "logs"
    s = subprocess.run([sys.executable, str(ASSETS / "score.py"), "--output", str(out), "--reference",
                        str(task / "tests" / "reference"), "--logs", str(logs)], capture_output=True, text=True,
                       env={"PYTHONPATH": str(task / "tests")})
    d = json.loads((logs / "details.json").read_text()) if (logs / "details.json").exists() else {}
    if s.returncode != 0 or not d:
        print("SCORER:", s.stdout[-500:], s.stderr[-2000:])
    if r.stdout.strip():
        print(r.stdout[-3000:])
    comp = {c["name"]: c["value"] for c in d.get("components", [])}
    info = d.get("info", {})
    print(f"== {task.name} {name}: score={d.get('score')} runtime_min={mins:.2f} failure={d.get('failure')}")
    print("   components", comp)
    print("   min_precision", info.get("min_precision"), "min_recall", info.get("min_recall"), "worst_af_err",
          info.get("worst_af_median_abs_error"), "pairing_wrong", info.get("pairing_wrong"))
    for k, v in info.items():
        if isinstance(v, dict):
            print("   ", k, {q: v[q] for q in ("tp", "fp", "fn", "gray_called", "precision", "recall", "af_median_abs_error")})
            if "-v" in sys.argv:
                print("      FP", v["false_positives"], "\n      FN", v["missed"])


if __name__ == "__main__":
    main()
