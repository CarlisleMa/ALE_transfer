"""Pilot harness: let an external agent attempt a task through `docker exec`, then grade it with the
task's own verifier image. Used for difficulty pilots when the Harbor agent/model of choice is not
available on the build host (no model API key). It is NOT a substitute for `harbor run`: the agent
under test has shell access only (no desktop MCP tools) and runs on the host, not in Harbor's sandbox.

Isolation protocol used for the audited pilots (scripts in `pipeline/pilot_isolation/`):
  1. `start` every task of the batch (agent image only; the instruction is copied into the
     container as /task/INSTRUCTION.md), remove all verifier images and the build cache;
  2. `lock.sh`: encrypt the repository, the scratch directory and any archive of them, verify the
     archive by a full decrypt-and-compare, delete the plaintext, and record the access times of
     archive and key (ext4 relatime: a read changes them);
  3. run the agents; `collect.sh` copies /task/output and the container's filesystem diff out;
  4. `unlock.sh`: report whether archive or key were read, restore the plaintext;
  5. `score` each collected output with a freshly built verifier image.
What this cannot exclude: reads of the session logs kept by the agent harness itself.
Where the plaintext may not be deleted (the level-4 round), `seal_overlay.sh` / `unseal_overlay.sh` are
used instead of steps 2 and 4: every file's access time is reset, each tree is hidden under an empty
tmpfs mount, and at the end the mounts are checked and every accessed entry is listed. That is
tamper-evident rather than tamper-proof.

  python -m alegen.pilot start <task_dir>                    # build agent image, start container pilot-<task>
  python -m alegen.pilot score <task_dir> [--output <dir>] [--label attempt2]   # grade (default: copy /task/output from the container)
  python -m alegen.pilot stop  <task_dir>
"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

from . import dockerlib, validate

ROOT = Path(__file__).resolve().parents[2]


def _have(tag: str) -> bool:
    return subprocess.run(["docker", "image", "inspect", tag], capture_output=True).returncode == 0


def main(argv: list[str]) -> int:
    if len(argv) < 2 or argv[0] not in ("start", "score", "stop"):
        print(__doc__)
        return 2
    task_dir = Path(argv[1]).resolve()
    tid = task_dir.name
    agent_tag, verifier_tag = validate.tags(tid)
    cname = f"pilot-{tid}"
    if argv[0] == "start":
        if not _have(agent_tag):
            dockerlib.build_image(task_dir / "environment", agent_tag)
        subprocess.run(["docker", "rm", "-f", cname], capture_output=True)
        subprocess.run(["docker", "run", "-d", "--name", cname, "--network", "none", "--cpus", "4",
                        "--memory", "12g", agent_tag, "sleep", "infinity"], check=True, capture_output=True)
        subprocess.run(["docker", "cp", str(task_dir / "instruction.md"), f"{cname}:/task/INSTRUCTION.md"],
                       check=True, capture_output=True)
        print(f"container {cname} is up")
        return 0
    if argv[0] == "stop":
        subprocess.run(["docker", "rm", "-f", cname], capture_output=True)
        return 0
    label = argv[argv.index("--label") + 1] if "--label" in argv else ""
    work = ROOT / "results" / "pilot" / (tid + (f"__{label}" if label else ""))
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    if "--output" in argv:
        src = Path(argv[argv.index("--output") + 1])
        if src.is_dir():
            shutil.copytree(src, work / "output")
    else:
        subprocess.run(["docker", "cp", f"{cname}:/task/output", str(work / "output")], capture_output=True)
    (work / "output").mkdir(exist_ok=True)
    if not _have(verifier_tag):
        dockerlib.build_image(task_dir / "tests", verifier_tag)
    reward, details = validate.verify(verifier_tag, work / "output", work)
    summary = {"task": tid, "reward": reward, "passed": details.get("passed"), "failure": details.get("failure"),
               "components": [(c["name"], c["value"]) for c in details.get("components", [])],
               "info": details.get("info")}
    (work / "pilot_result.json").write_text(json.dumps(summary, indent=1, default=str))
    print(json.dumps({k: summary[k] for k in ("task", "reward", "passed", "failure", "components")}, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
