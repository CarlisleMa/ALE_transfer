"""Project-level task: several instances of one family in one Harbor task; accepted only if every part is.

Why: in the pilots the failures sat in the tasks with the most scored conclusions, which suggests that
reliability over a long chain is a source of difficulty of its own. A project multiplies the number of
conclusions that must all be right, and the length of the session, without touching what a single part asks.
What it does NOT have (yet): coupling. The parts are independent instances (different seeds); nothing a part
produces is consumed by another one. A cross-family project needs a merged tool image and a shared
biological story in the generators.

  python -m alegen.project crispr_screen --seed 1 --parts 3 --level 4 [--no-validate] [--no-harbor]

Layout inside the task:  /task/input/part_X/ (inputs of part X + BRIEF.md = the part's instruction with the
paths rewritten), /task/output/part_X/ (deliverables of part X).  Score = min over the parts; the per-part
scores and their components are in details.json (usable as a dense reward).
"""
from __future__ import annotations

import argparse
import importlib
import json
import re
import sys
from pathlib import Path

from . import descaffold, harbor, validate
from .spec import Attack, Probe, TaskBuild

ROOT = Path(__file__).resolve().parents[2]

# How to run one part's reference code on that part's directories (family-specific, checked by the gates).
ADAPTERS = {
    "crispr_screen": {
        "oracle": "python3 /solution/{part}/crispr_l4_pipeline.py --input /task/input/{part} --output /task/output/{part} "
                  "--work /tmp/work_{part}",
        "alt_file": "solve_l4_alt.py",
        "alt_edit": ('Path("/task/input"), Path("/task/output"), Path("/tmp/l4_alt")',
                     'Path("/task/input/{part}"), Path("/task/output/{part}"), Path("/tmp/l4_alt_{part}")'),
        "alt": "python3 /solution/{part}/solve_l4_alt.py",
        # a part's scripted attack/probe is `python3 /alegen_script/<pipeline> <flags>`
        "script": lambda script, part: script.replace("/alegen_script/", f"/alegen_script/{part}/")
                  + f" --input /task/input/{part} --output /task/output/{part} --work /tmp/atk_{part}",
        "who": "computational biologist of a functional-genomics core facility",
        "what": "pooled CRISPR drug-modifier screens that a colleague who has left had started to analyse",
        "unit": "screen",
    },
}

PROJECT_SCORER = r'''#!/usr/bin/env python3
"""Project scorer: runs the scorer of every part on that part's output directory. Project score = minimum."""
import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

PARTS = __PARTS__
PASS = __PASS__

ap = argparse.ArgumentParser()
ap.add_argument("--output", default="/task/output")
ap.add_argument("--reference", default="/tests/reference")   # unused: every part has its own reference
ap.add_argument("--logs", default="/logs/verifier")
a = ap.parse_args()
logs = Path(a.logs)
logs.mkdir(parents=True, exist_ok=True)
here = Path(__file__).resolve().parent
parts = {}
for p in PARTS:
    pl = logs / p
    r = subprocess.run([sys.executable, str(here / p / "score.py"), "--output", str(Path(a.output) / p),
                        "--reference", str(here / p / "reference"), "--logs", str(pl)],
                       env=dict(os.environ, PYTHONPATH=str(here)), capture_output=True, text=True)
    if not (pl / "details.json").exists():      # a scorer crash is a verifier error, not a score of 0
        sys.stderr.write(f"scorer of {p} crashed:\n{r.stdout[-2000:]}\n{r.stderr[-4000:]}\n")
        sys.exit(1)
    parts[p] = json.loads((pl / "details.json").read_text())
score = min(float(d.get("score", 0.0)) for d in parts.values())
failures = [f"{p}: {d['failure']}" for p, d in parts.items() if d.get("failure")]
details = {"pass_threshold": PASS, "score": round(score, 6), "passed": score >= PASS,
           "aggregation": "minimum over the parts", "failure": "; ".join(failures) or None, "gates": [],
           "components": [{"name": p, "value": round(float(d.get("score", 0.0)), 6), "weight": round(1.0 / len(PARTS), 6),
                           "note": "score of this part on its own scorer; the project score is the minimum"}
                          for p, d in parts.items()],
           "info": {"parts": parts}}
(logs / "details.json").write_text(json.dumps(details, indent=1, default=str))
(logs / "reward.txt").write_text(f"{score:.6f}\n")
print(json.dumps({"score": details["score"], "passed": details["passed"], "failure": details["failure"]}))
'''


def _text(src) -> str:
    if isinstance(src, Path):
        return src.read_text()
    return src.decode() if isinstance(src, bytes) else src


def _by_name(items, name):
    for x in items:
        if x.name == name:
            return x
    raise KeyError(name)


def build(family: str, seed: int, n_parts: int, level: int, workdir: Path) -> TaskBuild:
    if family not in ADAPTERS:
        raise SystemExit(f"no project adapter for {family} (have: {sorted(ADAPTERS)})")
    ad = ADAPTERS[family]
    mod = importlib.import_module(f"alegen.families.{family}")
    names = [f"part_{chr(ord('A') + i)}" for i in range(n_parts)]
    parts: dict[str, TaskBuild] = {}
    for i, p in enumerate(names):
        b = mod.generate(seed=seed + i, level=min(level, 4), workdir=workdir / p)
        parts[p] = descaffold.apply(b, level) if level in (5, 6) else b
    first = parts[names[0]]
    if len({b.agent_dockerfile for b in parts.values()}) != 1 or len({b.verifier_dockerfile for b in parts.values()}) != 1:
        raise SystemExit("the parts do not share one tool image")

    def pre(p, d):  # prefix every file of a part
        return {f"{p}/{k}": v for k, v in d.items()}

    inputs, tests, sol = {}, {"score.py": PROJECT_SCORER.replace("__PARTS__", json.dumps(names))
                              .replace("__PASS__", repr(first.pass_threshold))}, {}
    for p, b in parts.items():
        brief = b.instruction_md.replace("/task/input", f"/task/input/{p}").replace("/task/output", f"/task/output/{p}")
        inputs.update(pre(p, b.input_files))
        inputs[f"{p}/BRIEF.md"] = brief.strip() + "\n"
        tests.update(pre(p, b.tests_files))
        files = {k: v for k, v in b.solution_files.items() if k != "solve.sh"}
        old, new = ad["alt_edit"]
        alt = _text(files[ad["alt_file"]])
        if alt.count(old) != 1:
            raise SystemExit(f"{family}: the alternative solution does not contain its path line exactly once")
        files[ad["alt_file"]] = alt.replace(old, new.format(part=p))
        sol.update(pre(p, files))
    loop = " ".join(names)
    sol["solve.sh"] = ("#!/bin/bash\n# Oracle of the project: the reference solution of every part on that part's directories.\n"
                       "# SOLUTION=2 selects the independent second implementation.\nset -euo pipefail\ncd /task\n"
                       f"for p in {loop}; do\n  mkdir -p /task/output/$p\n"
                       '  if [ "${SOLUTION:-1}" = "2" ]; then\n    ' + ad["alt"].format(part="$p") + "\n  else\n    "
                       + ad["oracle"].format(part="$p") + "\n  fi\ndone\n")

    def union(kind: str, name: str, special: dict | None = None):
        """Files of the named file-based probe/attack of every part (special: part -> other item name or None)."""
        out = {}
        for p, b in parts.items():
            n = (special or {}).get(p, name)
            if n is None:
                continue
            item = _by_name(b.probes if kind == "probe" and n == name else b.attacks + b.probes, n)
            if item.files is None:
                raise SystemExit(f"{n} of {p} is not file-based")
            out.update(pre(p, item.files))
        return out

    def scripted(name: str):
        script, files = [], {}
        for p, b in parts.items():
            item = _by_name(b.attacks + b.probes, name)
            script.append(f"mkdir -p /task/output/{p} && " + ad["script"](item.script, p))
            files.update(pre(p, item.script_files or {}))
        return " && ".join(script), files

    last = names[-1]
    beyond = _by_name(first.attacks, "lfc_drug_shifted_beyond_tolerance")
    l3_script, l3_files = scripted("l3_reference")
    probes = [
        Probe("planted_truth", "every part written directly from its planted truth", 0.999, files=union("probe", "planted_truth")),
        Probe("within_tolerance", "every part: the within-tolerance submission of that part", 0.99,
              files=union("probe", "within_tolerance")),
        Probe("independent_solution", "the independent second implementation on every part", 0.99,
              script="SOLUTION=2 bash /solution/solve.sh"),
    ]
    attacks = [
        Attack("l3_reference", "the level-3 reference method on every part", 0.75, script=l3_script, script_files=l3_files),
        Attack("one_part_beyond_tolerance", f"all parts from the planted truth except {last}, which is that part's "
               "`lfc_drug_shifted_beyond_tolerance` submission: the project must not pass on the average", beyond.max_score,
               files=union("probe", "planted_truth", {last: "lfc_drug_shifted_beyond_tolerance"})),
        Attack("one_part_missing", f"all parts from the planted truth, nothing delivered for {last}", 0.0,
               files=union("probe", "planted_truth", {last: None})),
    ]
    sw = re.search(r"Available software.*?(?=\n\n)", first.instruction_md, flags=re.S)
    unit, letters = ad["unit"], ", ".join(p[-1] for p in names)
    instruction = f"""You are the {ad['who']}. {n_parts} {ad['what']} have to be released: {unit}s {letters}.
They are independent experiments with their own samples, files and hand-over notes; they follow the same
facility rules.

Task directory:
- `/task`

For every {unit} X of {letters}:
- inputs (read-only): `/task/input/part_X/`
- the complete brief of that {unit}: `/task/input/part_X/BRIEF.md` (inputs, task, required outputs and
  constraints; it is written for that {unit} alone)
- deliverables: `/task/output/part_X/`

{sw.group(0) if sw else ''}

The project is released as a whole: it is accepted only if every {unit} is accepted on its own.

Constraints:
- Do not modify anything under `/task/input`; write deliverables only under `/task/output/part_X/`.
- Do not ask for confirmation. Execute directly.
"""
    tid = re.sub(r"-s\d+-l\d+$", "", first.task_id) + f"-project{n_parts}-s{seed:04d}-l{level}"
    return TaskBuild(
        task_id=tid, family=family, seed=seed, level=level,
        title=f"Project: {n_parts} x {first.title}",
        summary=f"Project of {n_parts} independent instances (seeds {seed}-{seed + n_parts - 1}) of: {first.summary} "
                "Accepted only if every part is accepted (score = minimum over the parts).",
        subdomain=first.subdomain, software=first.software, instruction_md=instruction,
        agent_dockerfile=first.agent_dockerfile, verifier_dockerfile=first.verifier_dockerfile,
        input_files=inputs, solution_files=sol, tests_files=tests,
        traps=first.traps, attacks=attacks, probes=probes, capabilities=first.capabilities + ["sustained reliability over a long session"],
        judge_explanation_md=("**Project task.** Every part is scored by the unchanged scorer of the single task it comes "
                              "from, on that part's own output directory and hidden reference. The project score is the "
                              "minimum of the part scores, so a project cannot pass on the average; `details.json` lists "
                              "every part with all its components.\n\n" + first.judge_explanation_md),
        pass_threshold=first.pass_threshold, agent_timeout_sec=int(first.agent_timeout_sec * 1.5),
        verifier_timeout_sec=first.verifier_timeout_sec * n_parts,
        oracle_steps=first.oracle_steps * n_parts, oracle_runtime_min=first.oracle_runtime_min * n_parts,
        cpus=first.cpus, memory_mb=first.memory_mb, storage_mb=first.storage_mb,
        extra_metadata={"project_parts": [b.task_id for b in parts.values()], "aggregation": "min"},
    )


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="alegen.project")
    ap.add_argument("family")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--parts", type=int, default=3)
    ap.add_argument("--level", type=int, default=4)
    ap.add_argument("--out", default=str(ROOT / "tasks"))
    ap.add_argument("--results", default=str(ROOT / "results"))
    ap.add_argument("--no-validate", action="store_true")
    ap.add_argument("--no-harbor", action="store_true")
    a = ap.parse_args(argv)
    work = Path(a.results) / "gen" / f"{a.family}-project{a.parts}-s{a.seed:04d}-l{a.level}"
    work.mkdir(parents=True, exist_ok=True)
    b = build(a.family, a.seed, a.parts, a.level, work)
    d = harbor.write_task(b, Path(a.out))
    print(f"wrote {d}")
    if a.no_validate:
        return 0
    rep = validate.validate(b, d, Path(a.results), harbor=not a.no_harbor)
    print(json.dumps({"task": b.task_id, "ok": rep["ok"], "minutes": rep["validation_minutes"]}))
    return 0 if rep["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
