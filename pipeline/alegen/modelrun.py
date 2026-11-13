"""Run tasks through Harbor with a real model on this build host.

  python -m alegen.modelrun <job_name> <model> <task_dir> [<task_dir> ...]
         [--agent terminus-2] [--n 2] [--k 1] [--ak key=value ...]

  e.g.  MOONSHOT_API_KEY=proxy-injected python -m alegen.modelrun kimi-k3 moonshot/kimi-k3 ../tasks/flow-compensation-s0001-l4

Why not plain `harbor run -p ../tasks`: this host reaches the internet only through a TLS-intercepting
proxy that containers cannot use. Therefore
  - the agent and verifier images are built here (dockerlib adds the proxy CA to a temporary copy of the
    Dockerfile) and injected as `docker_image` into a temporary copy of every task;
  - the agent phase runs with `network_mode = "no-network"` in the container. That is only possible with an
    agent whose model calls are made by the Harbor process on the host (terminus-2). Agents that are
    installed into the container (kimi-cli, gemini-cli, claude-code ...) need a host with normal networking,
    where the task's own `allowlist` setting applies and no image has to be pre-built.
The task content (instruction, inputs, tests, time limits) is not changed.
Results: <repo>/results/jobs/<job_name>/ ; summarise with `python -m alegen.calibrate harbor <dir> <label>`.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

from . import dockerlib, validate

ROOT = Path(__file__).resolve().parents[2]


def _have(tag: str) -> bool:
    return subprocess.run(["docker", "image", "inspect", tag], capture_output=True).returncode == 0


def prepare(task_dir: Path, dest: Path) -> Path:
    tid = task_dir.name
    agent_tag, verifier_tag = validate.tags(tid)
    need_agent, need_verifier = not _have(agent_tag), not _have(verifier_tag)
    if need_agent:
        base = validate._same_recipe_image(task_dir)      # same tool recipe already built here: only the inputs differ
        if not (base and dockerlib.derive_agent_image(base, task_dir / "environment", agent_tag)):
            dockerlib.build_image(task_dir / "environment", agent_tag)
    built = need_agent or need_verifier
    if need_verifier:
        dockerlib.build_image(task_dir / "tests", verifier_tag)
    if built and os.environ.get("ALEGEN_PRUNE_BUILD_CACHE"):
        # this host's image store keeps build cache of about the size of the image for every build
        subprocess.run(["docker", "builder", "prune", "-af"], capture_output=True)
    tmp = dest / tid
    if tmp.exists():
        shutil.rmtree(tmp)
    shutil.copytree(task_dir, tmp)
    toml = (tmp / "task.toml").read_text()
    toml = toml.replace("[verifier.environment]\n", f'[verifier.environment]\ndocker_image = "{verifier_tag}"\n')
    toml = toml.replace("[environment]\n", f'[environment]\ndocker_image = "{agent_tag}"\n')
    toml, n = re.subn(r'(\[agent\]\ntimeout_sec = [0-9.]+\n)network_mode = "allowlist"\nallowed_hosts = .*\n',
                      r'\1network_mode = "no-network"\n', toml)
    if n != 1:
        raise SystemExit(f"{tid}: could not switch the agent phase to no-network")
    (tmp / "task.toml").write_text(toml)
    return tmp


def main(argv: list[str]) -> int:
    if len(argv) < 3:
        print(__doc__)
        return 2
    job, model = argv[0], argv[1]
    agent, n, k, aks, tasks = "terminus-2", "2", "1", [], []
    it = iter(argv[2:])
    for a in it:
        if a == "--agent":
            agent = next(it)
        elif a == "--n":
            n = next(it)
        elif a == "--k":
            k = next(it)
        elif a == "--ak":
            aks.append(next(it))
        else:
            tasks.append(Path(a).resolve())
    jobs = ROOT / "results" / "jobs"
    dest = jobs / "_prepared" / job
    dest.mkdir(parents=True, exist_ok=True)
    for t in tasks:
        print("prepared", prepare(t, dest), flush=True)
    cmd = ["harbor", "run", "-p", str(dest if len(tasks) > 1 else dest / tasks[0].name), "-a", agent, "-m", model,
           "-o", str(jobs), "--job-name", job, "-n", n, "-k", k, "--yes"]
    for a in aks:
        cmd += ["--ak", a]
    env = dict(os.environ, PATH=os.path.expanduser("~/.local/bin") + ":" + os.environ["PATH"])
    print(" ".join(cmd), flush=True)
    return subprocess.run(cmd, env=env).returncode


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
