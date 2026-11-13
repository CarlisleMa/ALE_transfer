"""Difficulty calibration.

Two numbers per task, reported side by side (ALE itself only has the second, assigned after the fact):

  prior      in [0,1], from the structure of the reference solution (steps, tools, undisclosed
             decision points, runtime, GUI) - available the moment a task is generated, so a
             generator can aim at a difficulty before spending any model budget.
  posterior  1 - mean reward of a reference model over k trials (default: a "flash"-tier model),
             plus the full-pass rate, exactly the two columns of ALE's Table 1.

A task is kept for RL only if it is *learnable*: the weak model fails (posterior high) while a
strong model or the oracle demonstrates that reward is reachable, and partial credit is non-zero
for partial work.

Usage:
  python -m alegen.calibrate prior                       # table from tasks/*/task.toml
  python -m alegen.calibrate harbor <jobs_dir> [label]   # aggregate `harbor run` results
"""
from __future__ import annotations

import json
import sys
import tomllib
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def tier(d: float) -> str:
    return "last-exam" if d >= 0.70 else ("full-spectrum" if d >= 0.45 else "near-term")


def load_tasks(tasks_dir: Path) -> dict[str, dict]:
    out = {}
    for t in sorted(tasks_dir.glob("*/task.toml")):
        cfg = tomllib.loads(t.read_text())
        md = cfg.get("metadata", {})
        out[t.parent.name] = {
            "family": md.get("family"), "subdomain": md.get("subdomain"), "level": md.get("level"),
            "prior": md.get("difficulty_prior"), "tier_prior": md.get("difficulty_tier_prior"),
            "traps": md.get("traps", []), "software": md.get("software", []),
            "pass_threshold": md.get("pass_threshold", 0.9), "needs_gui": md.get("needs_gui", False),
            "oracle_steps": md.get("oracle_steps"), "timeout_sec": cfg.get("agent", {}).get("timeout_sec"),
        }
    return out


def prior_table(tasks_dir: Path) -> str:
    rows = ["| task | subdomain | level | traps | GUI | prior | tier (prior) |", "|---|---|---|---|---|---|---|"]
    for name, t in load_tasks(tasks_dir).items():
        rows.append(f"| `{name}` | {t['subdomain']} | {t['level']} | {len(t['traps'])} | "
                    f"{'yes' if t['needs_gui'] else 'no'} | {t['prior']} | {t['tier_prior']} |")
    return "\n".join(rows)


def harbor_results(jobs_dir: Path) -> dict[str, list[float | None]]:
    """task name -> list of rewards (None = trial errored / no reward)."""
    out = defaultdict(list)
    for res in jobs_dir.rglob("result.json"):
        try:
            j = json.loads(res.read_text())
        except Exception:
            continue
        if not isinstance(j, dict) or "task_name" not in j:
            continue
        name = str(j["task_name"]).split("/")[-1]
        vr = j.get("verifier_result") or {}
        rw = vr.get("rewards") if isinstance(vr, dict) else None
        if isinstance(rw, dict):
            rw = rw.get("reward")
        out[name].append(float(rw) if rw is not None else None)
    return out


def posterior_table(tasks_dir: Path, jobs_dir: Path, label: str) -> str:
    tasks = load_tasks(tasks_dir)
    res = harbor_results(jobs_dir)
    rows = [f"| task | prior | {label} mean score | {label} full-pass | trials (errors) | posterior | tier |",
            "|---|---|---|---|---|---|---|"]
    scores = []
    for name, t in tasks.items():
        r = res.get(name, [])
        ok = [x for x in r if x is not None]
        if not ok:
            rows.append(f"| `{name}` | {t['prior']} | not run | not run | 0 ({len(r)}) | - | {t['tier_prior']} (prior) |")
            continue
        mean = sum(ok) / len(ok)
        passed = sum(1 for x in ok if x >= t["pass_threshold"]) / len(ok)
        scores.append(mean)
        post = round(1 - mean, 3)
        rows.append(f"| `{name}` | {t['prior']} | {mean:.3f} | {passed:.0%} | {len(ok)} ({len(r) - len(ok)}) | "
                    f"{post} | {tier(post)} |")
    if scores:
        rows.append(f"| **mean over {len(scores)} tasks** | | **{sum(scores) / len(scores):.3f}** | | | | |")
    return "\n".join(rows)


def main(argv: list[str]) -> int:
    tasks_dir = ROOT / "tasks"
    if len(argv) >= 1 and argv[0] == "prior":
        print(prior_table(tasks_dir))
        return 0
    if len(argv) >= 2 and argv[0] == "harbor":
        print(posterior_table(tasks_dir, Path(argv[1]), argv[2] if len(argv) > 2 else "model"))
        return 0
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
