"""G8, the difficulty gate: keep a generated instance only if a yardstick model fails it while the task is sound.

Gates G1-G7 establish that an instance is correct (oracle and an independent solution reproduce the planted
truth, shortcuts score low). They say nothing about whether it is hard. G8 measures that with the model the
benchmark is aimed at:

    keep      the instance passed G1-G7 and the yardstick failed every one of its attempts (at least MIN)
    drop      the yardstick passed every attempt
    unstable  mixed outcomes: informative for training, but not a reliable "hard" item for evaluation
    pending   fewer than MIN attempts so far

Why "every attempt": single attempts of a flash model vary a lot (the same model scored 1.0 and 0.60-0.66 on
one task). Selecting on one failure would keep instances that are merely unlucky draws, and the pass rate
of the kept set would drift back up on the next run.
Why G1-G7 first: without them a filter of this kind selects broken or ambiguous tasks, because those are
failed most reliably.

The runs themselves are ordinary Harbor trials (`python -m alegen.modelrun`); a repeated attempt uses a job
name that ends in `__r<n>`. This module only reads results.

  python -m alegen.g8 gemini-3.8-flash [--min 2]      # table + results/difficulty/<model>.json
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
JOBS = ROOT / "results" / "jobs"


def attempts(model: str) -> dict[str, list[dict]]:
    """task -> attempts of the model (in job-name order), from the job patterns registered in model_runs.json."""
    cfg = json.loads((ROOT / "results" / "model_runs.json").read_text())
    m = cfg["models"][model]
    skip = m.get("not_counted", {})
    out: dict[str, list[dict]] = {}
    for pat in m["jobs"]:
        for jd in sorted(p for p in JOBS.glob(pat) if p.is_dir() and p.name != "_prepared"):
            for res in sorted(jd.glob("*/result.json")):
                j = json.loads(res.read_text())
                if "task_name" not in j:
                    continue
                task = str(j["task_name"]).split("/")[-1]
                rw = (j.get("verifier_result") or {}).get("rewards")
                rw = rw.get("reward") if isinstance(rw, dict) else rw
                if rw is None or f"{jd.name}/{task}" in skip:
                    continue
                low = None
                d = res.parent / "verifier" / "details.json"
                if d.exists():
                    dj = json.loads(d.read_text())
                    comps = dj.get("components") or []
                    if comps:
                        c = min(comps, key=lambda c: c["value"])
                        low = None if c["value"] == 1.0 else f"{c['name']} = {c['value']:.2f}"
                    elif dj.get("failure"):
                        low = str(dj["failure"])[:60]
                out.setdefault(task, []).append({"job": jd.name, "reward": rw, "lowest": low})
    return out


def validated(task: str) -> bool | None:
    f = ROOT / "results" / "validation" / f"{task}.json"
    return json.loads(f.read_text())["ok"] if f.exists() else None


def verdicts(model: str, min_attempts: int = 2) -> list[dict]:
    rows = []
    for task, att in sorted(attempts(model).items()):
        toml = ROOT / "tasks" / task / "task.toml"
        thr = 0.9
        if toml.exists():
            m = re.search(r"(?m)^pass_threshold\s*=\s*([0-9.]+)", toml.read_text())
            thr = float(m.group(1)) if m else thr
        passes = sum(1 for a in att if a["reward"] >= thr)
        ok = validated(task)
        if ok is not True:
            v = "not validated"
        elif len(att) < min_attempts:
            v = "pending"
        elif passes == 0:
            v = "keep"
        elif passes == len(att):
            v = "drop"
        else:
            v = "unstable"
        rows.append({"task": task, "attempts": len(att), "rewards": [round(a["reward"], 3) for a in att], "passes": passes,
                     "verdict": v, "lowest": [a["lowest"] for a in att]})
    return rows


def table(model: str, min_attempts: int = 2, only: list[str] | None = None) -> list[str]:
    rows = [r for r in verdicts(model, min_attempts) if not only or r["task"] in only]
    lines = [f"| task | attempts of `{model}` | rewards | passes | G8 verdict |", "|---|---|---|---|---|"]
    for r in rows:
        lines.append(f"| `{r['task']}` | {r['attempts']} | {', '.join(f'{x:.3f}' for x in r['rewards'])} | "
                     f"{r['passes']}/{r['attempts']} | {r['verdict']} |")
    return lines


def main(argv: list[str]) -> int:
    if not argv:
        print(__doc__)
        return 2
    model = argv[0]
    mn = int(argv[argv.index("--min") + 1]) if "--min" in argv else 2
    rows = verdicts(model, mn)
    out = ROOT / "results" / "difficulty" / f"{model}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"yardstick": model, "min_attempts": mn, "tasks": rows}, indent=1) + "\n")
    print("\n".join(table(model, mn)))
    counts = {}
    for r in rows:
        counts[r["verdict"]] = counts.get(r["verdict"], 0) + 1
    print("\n" + ", ".join(f"{k}: {v}" for k, v in sorted(counts.items())))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
