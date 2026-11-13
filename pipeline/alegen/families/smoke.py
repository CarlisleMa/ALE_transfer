"""Smallest possible family: used only to test the pipeline plumbing end to end."""
from __future__ import annotations

import json
import random
from pathlib import Path

from .. import dockerlib
from ..spec import Attack, Probe, TaskBuild

SCORE = '''\
import sys
sys.path.insert(0, "/tests")
from _scoring import run, to_float, close

def main(ctx):
    ref = ctx.ref_json("truth.json")
    ans = ctx.json("answer.json")
    ctx.gate(isinstance(ans, dict) and "total" in ans, "answer.json has key total")
    v = to_float(ans["total"], "total")
    ctx.component("total_matches", 1.0 if close(v, ref["total"], abs_=1e-6) else 0.0, 1.0)

run(main)
'''


def generate(seed: int, level: int, workdir: Path) -> TaskBuild:
    rng = random.Random(seed)
    nums = [rng.randint(1, 999) for _ in range(50)]
    return TaskBuild(
        task_id=f"smoke-s{seed:04d}", family="smoke", seed=seed, level=level,
        title="Smoke test", summary="Sum the numbers.", subdomain="none", software=["python"],
        instruction_md="Sum the integers in `/task/input/numbers.txt` and write "
                       '`/task/output/answer.json` as `{"total": <number>}`.',
        agent_dockerfile=dockerlib.agent_dockerfile(conda=[]),
        verifier_dockerfile=dockerlib.verifier_dockerfile(pip=[]),
        input_files={"numbers.txt": "\n".join(map(str, nums)) + "\n"},
        solution_files={"solve.sh": "#!/bin/bash\nset -e\npython3 -c \"import json;print(json.dumps({'total':sum(int(x) for x in open('/task/input/numbers.txt'))}))\" > /task/output/answer.json\n"},
        tests_files={"score.py": SCORE, "reference/truth.json": json.dumps({"total": sum(nums)})},
        attacks=[Attack("wrong_number", "plausible but wrong total", 0.0, files={"answer.json": '{"total": 1}'})],
        probes=[Probe("float_repr", "same value as float", 0.99, files={"answer.json": json.dumps({"total": float(sum(nums))})})],
        oracle_steps=1, pass_threshold=0.9,
    )
