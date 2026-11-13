"""Development helper: run assets/score.py on the host against an output directory and print the raw metrics."""
import json, subprocess, sys, tempfile
from pathlib import Path
HERE = Path(__file__).resolve().parent
SCORE = str(HERE.parent / "assets" / "score.py")
CORE = HERE.parents[2] / "scoring.py"      # copied into tests/ as _scoring.py by the pipeline
def score(out, ref):
    logs = tempfile.mkdtemp()
    lib = Path(tempfile.mkdtemp())
    (lib / "_scoring.py").write_text(CORE.read_text())
    r = subprocess.run([sys.executable, SCORE, "--output", str(out), "--reference", str(ref), "--logs", logs],
                       capture_output=True, text=True, env={"PYTHONPATH": str(lib), "PATH": "/usr/bin:/bin"})
    f = Path(logs) / "details.json"
    if not f.exists():
        return {"error": (r.stdout + r.stderr)[-800:]}
    return json.loads(f.read_text())
def line(name, d):
    if "error" in d: return f"{name:28s} VERIFIER ERROR {d['error']}"
    if d.get("failure"): return f"{name:28s} score=0 FAIL {d['failure'][:120]}"
    i = d["info"]; c = {x["name"]: x["value"] for x in d["components"]}
    return (f"{name:28s} score={d['score']:.3f} cnt={i['count_agreement_min']:.3f} sex={i['sex_correct']} lab={i['labels_valid']} "
            f"TP={i['de_tp']}/{i['de_planted']} FP={i['de_fp']} sens={i['de_sensitivity']:.3f} fdp={i['de_fdp']:.3f} lfcerr={i['lfc_median_abs_error']:.3f} "
            f"| comp fdp={c['de_false_discoveries']:.2f} sens={c['de_sensitivity']:.2f} lfc={c['de_log2fc']:.2f}")
if __name__ == "__main__":
    root, ref = Path(sys.argv[1]), Path(sys.argv[2])
    for d in sorted(p for p in root.iterdir() if p.is_dir()):
        print(line(d.name, score(d, ref)))
