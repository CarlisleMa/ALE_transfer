"""Refresh the result tables of README.md from results/ (blocks between <!-- BEGIN:name --> and <!-- END:name -->).

  cd pipeline && ../.venv/bin/python tools/readme_tables.py

Tables: `base` (ten base tasks x models, first attempt), `levels` (authored level and derived levels 5-7 x models,
first attempt), `g8` (yardstick attempts on the hardest version of every family, by seed), `project`.
Numbers come from the same readers as docs/TASKS.md (alegen.report, alegen.g8); nothing is typed by hand.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "pipeline"))
from alegen import g8, report  # noqa: E402

CFG = json.loads((ROOT / "results" / "model_runs.json").read_text())
MODELS = list(CFG["models"])
SHORT = {"gemini-3.8-flash": "Gemini 3.8 flash", "kimi-k3": "Kimi k3", "deepseek-v4.1-flash": "DeepSeek v4.1 flash",
         "qwen3.8-flash": "Qwen 3.8 flash", "glm-5.3-flash": "GLM 5.3 flash", "minimax-m3.1-flash": "MiniMax M3.1 flash"}
FAM = {"amplicon-16s": "amplicon", "md-system-repair": "md", "somatic-tn": "somatic", "crispr-screen": "crispr",
       "flow-compensation": "flow", "hcs-dose-response": "hcs", "rnaseq-confound": "rnaseq", "trio-liftover": "trio",
       "variant-review": "variant_review", "xlms-model-selection": "xlms"}


def first_attempts() -> tuple[dict, dict]:
    """model -> task -> trial (first attempt, valid score); model -> set of tasks attempted without a valid score."""
    valid, noscore = {}, {}
    for model, m in CFG["models"].items():
        jobs = []
        for pat in m["jobs"]:
            hits = sorted(p.name for p in (ROOT / "results" / "jobs").glob(pat)
                          if p.is_dir() and p.name != "_prepared" and not re.search(r"__r\d+", p.name))
            jobs += hits or [pat]
        skip = m.get("not_counted", {})
        trials = [t for j in jobs for t in report._trials(j)]
        ok = [t for t in trials if t["reward"] is not None and f"{t['job']}/{t['task']}" not in skip]
        valid[model] = {t["task"]: t for t in ok}
        noscore[model] = {t["task"] for t in trials} - set(valid[model])
    return valid, noscore


def name(task: str) -> str:
    m = re.match(r"(.+?)-(project\d+-)?s(\d+)-l(\d+)$", task)
    fam = FAM.get(m.group(1), m.group(1))
    return f"{fam}（L{m.group(4)}）" + ("" if m.group(3) == "0001" else f" 种子 {int(m.group(3))}")


def cell(valid, noscore, model, task) -> str:
    if task in valid[model]:
        return f"{valid[model][task]['reward']:.3f}"
    return "无有效分数" if task in noscore[model] else "未跑"


def table_base(valid, noscore) -> list[str]:
    pil = ROOT / "results" / "pilot"
    out = ["| 题目 | Claude Haiku（替身，每次） | " + " | ".join(SHORT[m] for m in MODELS) + " |", "|---|---|" + "---|" * len(MODELS)]
    for t in CFG["base_tasks"]:
        hk = [f"{json.loads((d / 'pilot_result.json').read_text())['reward']:.3f}" for d in sorted(pil.glob(t + "*"))
              if (d / "pilot_result.json").exists()]
        out.append(f"| {name(t)} | {'、'.join(hk)} | " + " | ".join(cell(valid, noscore, m, t) for m in MODELS) + " |")
    done = {m: [valid[m][t]["reward"] for t in CFG["base_tasks"] if t in valid[m]] for m in MODELS}
    out.append("| **有效 / 通过 / 平均** | | " + " | ".join(
        f"**{len(v)} / {sum(x >= 0.9 for x in v)} / {sum(v) / len(v):.3f}**" if v else "-" for v in done.values()) + " |")
    return out


def table_levels(valid, noscore) -> list[str]:
    out = ["| 家族与档位 | " + " | ".join(SHORT[m] for m in MODELS) + " |", "|---|" + "---|" * len(MODELS)]
    for pair in CFG["level5"]["pairs"]:
        for key in ("l4", "l5", "l6", "l7"):
            t = pair.get(key)
            if t:
                out.append(f"| {name(t)} | " + " | ".join(cell(valid, noscore, m, t) for m in MODELS) + " |")
    return out


def table_project(valid, noscore) -> list[str]:
    pj = CFG["level5"]["project"]
    out = ["| 模型 | 单题（crispr L4） | 三题项目 | 三部分得分 |", "|---|---|---|---|"]
    for m in MODELS:
        b = valid[m].get(pj["task"])
        if b:
            parts = " / ".join(f"{v:.3f}" for _, v in sorted((b.get("parts") or {}).items()))
            out.append(f"| {SHORT[m]} | {cell(valid, noscore, m, pj['part_task'])} | {b['reward']:.3f} | {parts} |")
    return out


def table_g8() -> list[str]:
    rows = g8.verdicts(CFG["g8"]["yardstick"], CFG["g8"].get("min_attempts", 2))
    by_task = {r["task"]: r for r in rows}
    hardest = CFG["g8"]["hardest"]          # family slug -> level of its hardest validated version
    zh = {"keep": "保留", "drop": "淘汰", "unstable": "不稳定", "pending": "不足两次", "not validated": "未过关卡"}
    out = ["| 家族 | 最难版本 | 实例（种子）：两次得分 → 判定 | 家族结论 |", "|---|---|---|---|"]
    for slug, lvl in hardest.items():
        inst = sorted(t for t in by_task if re.match(rf"{slug}-s\d+-l{lvl}$", t))
        cells, verdicts = [], []
        for t in inst:
            r = by_task[t]
            seed = int(re.search(r"-s(\d+)-", t).group(1))
            cells.append(f"种子 {seed}：{'、'.join(f'{x:.3f}' for x in r['rewards'])} → {zh[r['verdict']]}")
            verdicts.append(r["verdict"])
        verdicts = [v for v in verdicts if v != "pending"] or verdicts      # an instance with one attempt decides nothing
        if not verdicts:
            fam = "未跑"
        elif all(v == "pending" for v in verdicts):
            fam = "不足两次"
        elif all(v == "keep" for v in verdicts):
            fam = "拦住"
        elif all(v == "drop" for v in verdicts):
            fam = "没拦住"
        else:
            fam = "部分拦住" if "keep" in verdicts else "不确定"
        out.append(f"| {FAM[slug]} | level {lvl} | {'；'.join(cells) or '未跑'} | {fam} |")
    return out


def line_cost() -> list[str]:
    f = ROOT / "results" / "billing.json"       # written by tools/host_runs/billing.py from the relay logs
    if not f.exists():
        return ["  （没有账单汇总文件）"]
    b = json.loads(f.read_text())
    parts = [f"{SHORT[m]} {b['models'][m]['billed_usd']:.2f}" for m in MODELS if m in b["models"]]
    return [f"  {'、'.join(parts)}，合计 **{b['total_usd']:.2f} 美元**（截至 {b['as_of_utc']} UTC）。"]


def main() -> int:
    valid, noscore = first_attempts()
    blocks = {"base": table_base(valid, noscore), "levels": table_levels(valid, noscore),
              "project": table_project(valid, noscore), "g8": table_g8(), "cost": line_cost()}
    readme = ROOT / "README.md"
    s = readme.read_text()
    for k, lines in blocks.items():
        pat = re.compile(rf"(<!-- BEGIN:{k} -->\n).*?(<!-- END:{k} -->)", re.S)
        if not pat.search(s):
            print(f"README has no block {k}; table follows\n" + "\n".join(lines) + "\n")
            continue
        s = pat.sub(lambda m: m.group(1) + "\n".join(lines) + "\n" + m.group(2), s)
    readme.write_text(s)
    return 0


if __name__ == "__main__":
    sys.exit(main())
