"""Fallback when no image may be built (shared disk below the 8-GB rule): run the checks of validate.py
(G2 oracle, G5 determinism, G3 probes, G4 attacks) for one task WITHOUT building its agent/verifier images.
Scripts run in the family's generation tools image (same conda package list as the agent image) with the task's
input directory mounted read-only at /task/input; scoring runs tests/score.py with the host venv.
Not a substitute for G1 (images build), G6 (Harbor) and G7 (login shell)."""
import json, shutil, subprocess, sys, tempfile, time, os
from pathlib import Path
from alegen import dockerlib, harbor, validate
from alegen.families import trio_liftover as T
from alegen.families.trio_liftover import refdata

seed, level, base = int(sys.argv[1]), int(sys.argv[2]), Path(sys.argv[3])
b = T.generate(seed, level, base / "gen")
d = harbor.write_task(b, base / "tasks")
tools = dockerlib.tool_image(conda=refdata.TOOLS)      # exists already; nothing is built
work = base / "work" / b.task_id
shutil.rmtree(work, ignore_errors=True)
pyscore = base / "pyscore"
pyscore.mkdir(parents=True, exist_ok=True)
shutil.copyfile(Path(validate.__file__).parent / "scoring.py", pyscore / "_scoring.py")
ok_all = True

def run_script(script, out, files=None):
    sdir = Path(tempfile.mkdtemp(prefix="emu-"))
    validate._materialize(files or {}, sdir)
    (sdir / "run.sh").write_text(script)
    try:
        return dockerlib.run_in_image(tools, "mkdir -p /task/output && cd /task && bash /alegen_script/run.sh",
                                      mounts={sdir: "/alegen_script", d / "solution": "/solution",
                                              d / "environment" / "input": "/task/input"},
                                      copy_out={"/task/output": out}, timeout=3600)
    finally:
        shutil.rmtree(sdir, ignore_errors=True)

def score(out, logs):
    out.mkdir(parents=True, exist_ok=True)
    shutil.rmtree(logs, ignore_errors=True)
    subprocess.run([sys.executable, str(d / "tests" / "score.py"), "--output", str(out), "--reference",
                    str(d / "tests" / "reference"), "--logs", str(logs)], capture_output=True, text=True,
                   env=dict(os.environ, PYTHONPATH=str(pyscore)))
    if not (logs / "reward.txt").exists():
        return None, {}
    return float((logs / "reward.txt").read_text()), json.loads((logs / "details.json").read_text())

def record(gate, name, ok, **kw):
    global ok_all
    ok_all = ok_all and bool(ok)
    print(f"  [{'ok' if ok else 'FAIL'}] {gate} {name} " + " ".join(f"{k}={v}" for k, v in kw.items()), flush=True)

t0 = time.time()
oo = work / "oracle" / "output"
t1 = time.time()
run_script("bash /solution/solve.sh", oo)
omin = round((time.time() - t1) / 60, 2)
s, det = score(oo, work / "oracle" / "logs")
record("G2*", "oracle", s is not None and s >= 0.99, score=s, runtime_min=omin, failure=det.get("failure") or "")
s2, _ = score(oo, work / "oracle2" / "logs")
record("G5*", "verifier_deterministic", s == s2, first=s, second=s2)
for p in b.probes:
    out = work / "probe" / p.name / "output"
    if p.script is not None:
        run_script(p.script, out, p.script_files)
    else:
        validate._materialize(p.files, out)
    s, det = score(out, work / "probe" / p.name / "logs")
    record("G3*", "probe:" + p.name, s is not None and s >= p.min_score, score=s, floor=p.min_score, failure=det.get("failure") or "")
record("G4*", "has_l3_reference_attack", any(a.name == "l3_reference" for a in b.attacks))
for a in validate.generic_attacks(oo) + list(b.attacks):
    out = work / "attack" / a.name / "output"
    if a.script is not None:
        run_script(a.script.replace("cp -r /task/input/. /task/output/", "cp -r /task/input/. /task/output/"), out, a.script_files)
    else:
        validate._materialize(a.files, out)
    s, det = score(out, work / "attack" / a.name / "logs")
    record("G4*", "attack:" + a.name, s is not None and s <= a.max_score, score=s, ceiling=a.max_score, failure=(det.get("failure") or "")[:80])
print(json.dumps({"task": b.task_id, "ok_without_image_gates": ok_all, "minutes": round((time.time() - t0) / 60, 2)}))
