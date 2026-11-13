"""Build docs/TASKS.md from what is actually on disk (tasks/*/task.toml + results/validation/*.json),
so the report never states a result that was not measured.

  python -m alegen.report            # writes docs/TASKS.md and prints it
"""
from __future__ import annotations

import json
import re
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def _trials(job: str) -> list[dict]:
    """One dict per trial of a Harbor job directory under results/jobs/."""
    out = []
    for res in sorted((ROOT / "results" / "jobs" / job).glob("*/result.json")):
        try:
            j = json.loads(res.read_text())
        except Exception:
            continue
        if "task_name" not in j:
            continue
        vr = j.get("verifier_result") or {}
        rw = vr.get("rewards")
        rw = rw.get("reward") if isinstance(rw, dict) else rw
        ar = j.get("agent_result") or {}
        ex = j.get("exception_info") or {}
        turns = None
        tj = res.parent / "agent" / "trajectory.json"
        if tj.exists():
            try:
                turns = sum(1 for st in json.loads(tj.read_text()).get("steps") or [] if st.get("source") == "agent")
            except Exception:
                pass
        low, parts = None, None
        for d in sorted((res.parent / "verifier").glob("details.json")):
            dj = json.loads(d.read_text())
            comps = dj.get("components") or []
            if comps:
                c = min(comps, key=lambda c: c["value"])
                low = None if c["value"] == 1.0 else f"{c['name']} = {c['value']:.2f}"
            if dj.get("aggregation"):        # project task: one component per part
                parts = {c["name"]: c["value"] for c in comps}
        secs = None
        if j.get("started_at") and j.get("finished_at"):
            from datetime import datetime
            f = lambda x: datetime.fromisoformat(str(x).replace("Z", "+00:00"))
            secs = int((f(j["finished_at"]) - f(j["started_at"])).total_seconds())
        out.append({"task": str(j["task_name"]).split("/")[-1], "job": job, "reward": rw, "secs": secs, "turns": turns,
                    "in": ar.get("n_input_tokens"), "cache": ar.get("n_cache_tokens"), "out": ar.get("n_output_tokens"),
                    "cost": ar.get("cost_usd"), "exc": ex.get("exception_type") if isinstance(ex, dict) else None,
                    "low": low, "parts": parts})
    return out


def model_runs() -> list[str]:
    """Harbor runs with real models (results/model_runs.json names the job directories of each model).

    A trial listed under `not_counted` there (keyed `<job>/<task>`) is shown with what Harbor recorded and
    the reason, and is left out of the mean, the pass count and the side-by-side table."""
    cfg_f = ROOT / "results" / "model_runs.json"
    if not cfg_f.exists():
        return []
    cfg = json.loads(cfg_f.read_text())
    lines = ["", "## Harbor runs with real models", "", cfg.get("_note", ""), ""]
    per_model, all_valid, attempted = {}, {}, {}
    for model, m in cfg["models"].items():
        jobs = []                      # an entry may be a glob: one Harbor job per trial is `<model>__*`
        for pat in m["jobs"]:
            hits = sorted(p.name for p in (ROOT / "results" / "jobs").glob(pat)      # __r<n> = repeated attempt (G8 table)
                          if p.is_dir() and p.name != "_prepared" and not re.search(r"__r\d+$", p.name))
            jobs += hits or [pat]
        trials = [t for job in jobs for t in _trials(job)]
        skip = m.get("not_counted", {})
        key = lambda t: f"{t['job']}/{t['task']}"
        valid = [t for t in trials if t["reward"] is not None and key(t) not in skip]
        base = set(cfg.get("base_tasks") or [t["task"] for t in valid])
        scored = [t for t in valid if t["task"] in base]
        others = [t for t in trials if t not in valid]
        per_model[model] = {t["task"]: t for t in scored}
        all_valid[model] = {t["task"]: t for t in valid}
        attempted[model] = {t["task"] for t in trials}
        lines += [f"### `{model}` ({m.get('agent', 'terminus-2')})", "", m.get("note", ""), "",
                  "| task | reward | pass | wall time | agent turns | input tokens (cached) | output tokens | cost estimate | lowest component |",
                  "|---|---|---|---|---|---|---|---|---|"]
        for t in sorted(scored, key=lambda t: t["task"]):
            lines.append(
                f"| `{t['task']}` | {t['reward']:.3f} | {'yes' if t['reward'] >= 0.9 else 'no'} | "
                f"{t['secs'] // 60} min {t['secs'] % 60} s | {t['turns']} | {(t['in'] or 0) / 1e6:.2f} M ({(t['cache'] or 0) / 1e6:.2f} M) | "
                f"{(t['out'] or 0) / 1e3:.0f} k | ${t['cost'] or 0:.2f} | {t['low'] or '-'} |")
        if scored:
            rw = [t["reward"] for t in scored]
            lines.append(f"| **{len(rw)} tasks** | **mean {sum(rw) / len(rw):.3f}** | **{sum(1 for x in rw if x >= 0.9)}/{len(rw)}** | "
                         f"{sum(t['secs'] for t in scored) // 60} min | {sum(t['turns'] or 0 for t in scored)} | "
                         f"{sum(t['in'] or 0 for t in scored) / 1e6:.1f} M | {sum(t['out'] or 0 for t in scored) / 1e3:.0f} k | "
                         f"${sum(t['cost'] or 0 for t in scored):.2f} | |")
        if others or m.get("not_run"):
            lines += ["", "Not counted (no valid score):", ""]
            for t in sorted(others, key=lambda t: t["task"]):
                rec = "no reward" if t["reward"] is None else f"reward {t['reward']:.3f}"
                wall = f"{t['secs'] // 60} min" if t["secs"] is not None else "unfinished"
                why = skip.get(key(t)) or "ended before the verifier ran."
                lines.append(f"- `{t['task']}` (job `{t['job']}`; Harbor recorded {rec}, `{t['exc'] or 'no exception'}`, {wall}, "
                             f"{t['turns'] or 0} agent turns, ${t['cost'] or 0:.2f}): {why}")
            for task, why in (m.get("not_run") or {}).items():
                lines.append(f"- `{task}`: {why}")
            if others:
                lines.append(f"\nCost estimate of all trials of this model, counted or not: ${sum(t['cost'] or 0 for t in trials):.2f}.")
        lines.append("")
    # side by side with the pilot model
    pil = ROOT / "results" / "pilot"
    names = sorted({t for d in per_model.values() for t in d})
    if names:
        lines += ["### Side by side", "",
                  "Claude Haiku: pilot harness (`docker exec`), every attempt listed. The other columns: Harbor, one attempt; "
                  "`no valid score` refers to the list under that model's table.", "",
                  "| task | Claude Haiku (pilot) | " + " | ".join(f"`{m}`" for m in per_model) + " |",
                  "|---|---|" + "---|" * len(per_model)]
        for n in names:
            hk = []
            for d in sorted(pil.glob(n + "*")):
                rf = d / "pilot_result.json"
                if rf.exists():
                    hk.append(f"{json.loads(rf.read_text())['reward']:.3f}")
            lines.append(f"| `{n}` | {', '.join(hk) or 'not run'} | " +
                         " | ".join(f"{per_model[m][n]['reward']:.3f}" if n in per_model[m]
                                    else ("no valid score" if n in attempted[m] else "not run") for m in per_model) + " |")
    lines += _level5(cfg, all_valid)
    g8cfg = cfg.get("g8")
    if g8cfg:
        from . import g8
        lines += ["", "### G8: difficulty gate (repeated attempts of the yardstick model)", "", g8cfg.get("note", ""), ""]
        lines += g8.table(g8cfg["yardstick"], g8cfg.get("min_attempts", 2), g8cfg.get("tasks"))
    return lines


def _level5(cfg: dict, valid: dict) -> list[str]:
    """Level 4 against level 5 (scaffolding removed) and single task against project, per model."""
    l5 = cfg.get("level5")
    if not l5:
        return []
    resc_f = ROOT / "results" / "rescored" / "l4_outputs_with_l5_scorer.json"
    resc = {}
    if resc_f.exists():
        for k, v in json.loads(resc_f.read_text()).items():
            resc[(k.split("/")[0], v["task"])] = v["rescored"]
    fmt = lambda x: "not run" if x is None else f"{x:.3f}"
    cell = lambda t: "not run" if t is None else f"{t['reward']:.3f} ({t['secs'] // 60} min, {t['turns']} turns" + (f"; {t['low']}" if t["low"] else "") + ")"
    lines = ["", "### Hardest authored level against the derived levels 5, 6 and 7", "", l5.get("note", ""), "",
             "| family | model | authored level | its output under the level-5/6 scorer | level 5 | level 6 | level 7 |",
             "|---|---|---|---|---|---|---|"]
    for pair in l5["pairs"]:
        for model, tr in valid.items():
            a, b, c, d = (tr.get(pair.get(k, "")) for k in ("l4", "l5", "l6", "l7"))
            if a is None and b is None and c is None and d is None:
                continue
            same = fmt(a["reward"]) if a and not pair.get("scorer_changed") else fmt(resc.get((a["job"], a["task"]))) if a else "not run"
            na = lambda key, x: cell(x) if pair.get(key) else "-"
            lines.append(f"| {pair['family']} | `{model}` | {fmt(a['reward'] if a else None)} | {same} | {na('l5', b)} | {na('l6', c)} | {na('l7', d)} |")
    pj = l5.get("project")
    if pj:
        lines += ["", "### Single task against project (several instances in one task, score = minimum)", "", pj.get("note", ""), "",
                  "| model | single task (`" + pj["part_task"] + "`) | project (`" + pj["task"] + "`) | parts | time | turns |",
                  "|---|---|---|---|---|---|"]
        for model, tr in valid.items():
            a, b = tr.get(pj["part_task"]), tr.get(pj["task"])
            if b is None:
                continue
            parts = ", ".join(f"{k[-1]} {v:.3f}" for k, v in sorted((b.get("parts") or {}).items()))
            lines.append(f"| `{model}` | {fmt(a['reward'] if a else None)} | {b['reward']:.3f} | {parts} | {b['secs'] // 60} min | {b['turns']} |")
    return lines


def main() -> int:
    tasks = sorted((ROOT / "tasks").glob("*/task.toml"))
    val_dir = ROOT / "results" / "validation"
    out = ["# Tasks and validation status", "",
           "Generated by `python -m alegen.report` from `tasks/*/task.toml` and `results/validation/*.json`.",
           "A blank cell means that check has not been run; nothing here is estimated.", "",
           "## Overview", "",
           "| task | ALE subdomain | software | undisclosed decision points | GUI | prior difficulty | prior tier |",
           "|---|---|---|---|---|---|---|"]
    detail = ["", "## Release gates", "",
              "G2 oracle, G3 probes (must score high), G4 attacks (must score low), G5 determinism, G6 Harbor, G7 login shell.",
              "For level-4 tasks the column `level-3 method` is the score of the level-3 reference method on the level-4 inputs.", "",
              "| task | all gates | oracle | oracle min | probes (min score) | attacks: count (highest score / its ceiling) | level-3 method | harbor oracle | harbor nop |",
              "|---|---|---|---|---|---|---|---|---|"]
    attack_rows = ["", "## Attack suite, per task", ""]
    for t in tasks:
        name = t.parent.name
        cfg = tomllib.loads(t.read_text())
        md = cfg.get("metadata", {})
        out.append(f"| `{name}` | {md.get('subdomain', '')} | {', '.join(md.get('software', []))} | "
                   f"{', '.join(md.get('traps', [])) or 'none'} | {'yes' if md.get('needs_gui') else 'no'} | "
                   f"{md.get('difficulty_prior', '')} | {md.get('difficulty_tier_prior', '')} |")
        vf = val_dir / f"{name}.json"
        if not vf.exists():
            detail.append(f"| `{name}` | not validated | | | | | | | |")
            continue
        rep = json.loads(vf.read_text())
        ch = rep["checks"]

        def get(gate, prefix=""):
            return [c for c in ch if c["gate"] == gate and c["name"].startswith(prefix)]

        orc = get("G2")
        probes = get("G3")
        attacks = get("G4")
        ho = get("G6", "harbor_oracle")
        hn = get("G6", "harbor_nop")
        pmin = min((c.get("score") if c.get("score") is not None else -1) for c in probes) if probes else None
        attacks = [c for c in attacks if c["name"].startswith("attack:")]
        l3 = get("G4", "attack:l3_reference")
        worst = max(attacks, key=lambda c: (c.get("score") or 0)) if attacks else None   # strongest shortcut
        detail.append(
            f"| `{name}` | {'PASS' if rep['ok'] else 'FAIL'} | {orc[0].get('score') if orc else ''} | "
            f"{rep.get('oracle_runtime_min', '')} | {len(probes)} ({pmin}) | "
            f"{len(attacks)} ({worst.get('score') if worst else ''} / {worst.get('ceiling') if worst else ''} on `{worst['name'][7:] if worst else ''}`) | "
            f"{l3[0].get('score') if l3 else ''} | {ho[0].get('rewards') if ho else ''} | {hn[0].get('rewards') if hn else ''} |")
        attack_rows += [f"### `{name}`", "", "| check | score | bound | ok | note |", "|---|---|---|---|---|"]
        for c in probes:
            attack_rows.append(f"| probe `{c['name'][6:]}` | {c.get('score')} | >= {c.get('floor')} | "
                               f"{'yes' if c['ok'] else 'NO'} | {str(c.get('failure') or '')[:80]} |")
        for c in attacks:
            attack_rows.append(f"| attack `{c['name'][7:]}` | {c.get('score')} | <= {c.get('ceiling')} | "
                               f"{'yes' if c['ok'] else 'NO'} | {str(c.get('failure') or '')[:80]} |")
        attack_rows.append("")
    # every (family, seed, level) that went through the gates, whether or not its task dir is shipped
    scale = ["", "## Seeds and levels that passed the gates", "",
             "Evidence that a family is a generator, not a single task. `harbor` = G6 was run for that instance.", "",
             "| family | seed | level | all gates | checks | oracle | harbor |", "|---|---|---|---|---|---|---|"]
    for vf in sorted(val_dir.glob("*.json")):
        rep = json.loads(vf.read_text())
        orc = [c for c in rep["checks"] if c["gate"] == "G2"]
        hb = [c for c in rep["checks"] if c["gate"] == "G6"]
        scale.append(f"| {rep.get('family')} | {rep.get('seed')} | {rep.get('level')} | "
                     f"{'PASS' if rep['ok'] else 'FAIL'} | {len(rep['checks'])} | "
                     f"{orc[0].get('score') if orc else ''} | {'yes' if hb else 'no'} |")
    # small-model pilot (results/pilot): measured, one attempt each, not a Harbor run
    pilot = []
    runs_f = ROOT / "results" / "pilot" / "runs.json"
    if runs_f.exists():
        runs = json.loads(runs_f.read_text())
        pilot = ["", "## Small-model pilot (one attempt per task)", "", runs.get("_note", ""), "",
                 "| task | reward | full pass | posterior difficulty (1 - reward) | wall time | tool calls | lowest component |",
                 "|---|---|---|---|---|---|---|"]
        rewards = []
        for name, meta in runs["runs"].items():
            rf = ROOT / "results" / "pilot" / name / "pilot_result.json"
            if not rf.exists():
                pilot.append(f"| `{name}` | not scored | | | | | |")
                continue
            r = json.loads(rf.read_text())
            comps = r.get("components") or []
            low = min(comps, key=lambda c: c[1]) if comps else ("", "")
            rewards.append(r["reward"])
            pilot.append(f"| `{name}` | {r['reward']:.3f} | {'yes' if r.get('passed') else 'no'} | {1 - r['reward']:.3f} | "
                         f"{meta['wall_seconds'] // 60} min {meta['wall_seconds'] % 60} s | {meta['tool_calls']} | "
                         f"{'-' if low[1] == 1.0 else f'{low[0]} = {low[1]:.2f}'} |")
        if rewards:
            pilot.append(f"| **{len(rewards)} tasks** | **mean {sum(rewards) / len(rewards):.3f}** | "
                         f"**{sum(1 for x in rewards if x >= 0.9)}/{len(rewards)}** | | | | |")
        second = runs.get("second_attempts", {})
        if second:
            pilot += ["", second.get("_note", ""), "",
                      "| task | first attempt | second attempt | wall time (2nd) | tool calls (2nd) | lowest component (2nd) |",
                      "|---|---|---|---|---|---|"]
            for name, meta in second.items():
                if name.startswith("_"):
                    continue
                f1 = ROOT / "results" / "pilot" / name / "pilot_result.json"
                f2 = ROOT / "results" / "pilot" / f"{name}__attempt2" / "pilot_result.json"
                if not (f1.exists() and f2.exists()):
                    continue
                r1, r2 = json.loads(f1.read_text()), json.loads(f2.read_text())
                low = min(r2.get("components") or [("", 1.0)], key=lambda c: c[1])
                pilot.append(f"| `{name}` | {r1['reward']:.3f} | {r2['reward']:.3f} | {meta['wall_seconds'] // 60} min "
                             f"{meta['wall_seconds'] % 60} s | {meta['tool_calls']} | {low[0]} = {low[1]:.2f} |")
        # level-4 round: same harness, own time cap; then the level-3 -> level-4 ladder per family
        l4 = runs.get("level4", {})
        if l4.get("runs"):
            pilot += ["", "### Level 4, first attempt", "", l4.get("_note", ""), "",
                      "| task | reward | full pass | posterior difficulty (1 - reward) | wall time | tool calls | lowest component |",
                      "|---|---|---|---|---|---|---|"]
            got = {}
            for name, meta in l4["runs"].items():
                rf = ROOT / "results" / "pilot" / name / "pilot_result.json"
                if not rf.exists():
                    pilot.append(f"| `{name}` | not scored | | | | | |")
                    continue
                r = json.loads(rf.read_text())
                comps = r.get("components") or []
                low = min(comps, key=lambda c: c[1]) if comps else ("", 1.0)
                got[name] = (r, meta)
                pilot.append(f"| `{name}` | {r['reward']:.3f} | {'yes' if r.get('passed') else 'no'} | {1 - r['reward']:.3f} | "
                             f"{meta['wall_seconds'] // 60} min {meta['wall_seconds'] % 60} s | {meta['tool_calls']} | "
                             f"{'-' if low[1] == 1.0 else f'{low[0]} = {low[1]:.2f}'} |")
            if got:
                rw = [r["reward"] for r, _ in got.values()]
                pilot.append(f"| **{len(rw)} tasks** | **mean {sum(rw) / len(rw):.3f}** | "
                             f"**{sum(1 for x in rw if x >= 0.9)}/{len(rw)}** | | | | |")
                pilot += ["", "### Difficulty ladder: the same family and seed at level 3 and level 4", "",
                          "| family | level 3 reward | level 3 time | level 3 tool calls | level 4 reward | level 4 time | level 4 tool calls |",
                          "|---|---|---|---|---|---|---|"]
                for name, (r4, m4) in got.items():
                    n3 = name[:-1] + "3"
                    f3 = ROOT / "results" / "pilot" / n3 / "pilot_result.json"
                    m3 = runs["runs"].get(n3)
                    if not (f3.exists() and m3):
                        continue
                    r3 = json.loads(f3.read_text())
                    pilot.append(f"| {name.rsplit('-s', 1)[0]} | {r3['reward']:.3f} | {m3['wall_seconds'] // 60} min {m3['wall_seconds'] % 60} s | "
                                 f"{m3['tool_calls']} | {r4['reward']:.3f} | {m4['wall_seconds'] // 60} min {m4['wall_seconds'] % 60} s | {m4['tool_calls']} |")
        for key, title in (("full_limit", "### Level 3 re-run with the task's own time limit"),
                           ("level4_second", "### Level 4, second attempt")):
            extra = runs.get(key, {})
            names = [k for k in extra if not k.startswith("_")]
            if not names:
                continue
            pilot += ["", title, "", extra.get("_note", ""), "",
                      "| task | reward | wall time | tool calls | lowest component |", "|---|---|---|---|---|"]
            for name in names:
                meta = extra[name]
                rf = ROOT / "results" / "pilot" / f"{name}__{meta['label']}" / "pilot_result.json"
                if not rf.exists():
                    continue
                r = json.loads(rf.read_text())
                low = min(r.get("components") or [("", 1.0)], key=lambda c: c[1])
                pilot.append(f"| `{name}` | {r['reward']:.3f} | {meta['wall_seconds'] // 60} min {meta['wall_seconds'] % 60} s | "
                             f"{meta['tool_calls']} | {'-' if low[1] == 1.0 else f'{low[0]} = {low[1]:.2f}'} |")
    models = model_runs()
    text = "\n".join(out + detail + models + pilot + scale + attack_rows) + "\n"
    (ROOT / "docs").mkdir(exist_ok=True)
    (ROOT / "docs" / "TASKS.md").write_text(text)
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
