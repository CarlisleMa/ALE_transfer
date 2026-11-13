"""Release gate for a generated task. A task ships only if ALL of these hold:

  G1  images build                      (agent image + dedicated verifier image)
  G2  oracle solution scores >= 0.99    (task is solvable as written, inside the agent image, offline)
  G3  every probe scores >= its floor   (independent second solution, in-tolerance perturbations)
  G4  every attack scores <= its ceiling (empty, inputs copied, emptied files, family shortcuts,
                                          and the trap-ignoring "naive pipeline")
  G5  verifier is deterministic         (same output scored twice -> same reward)
  G6  Harbor agrees                     (`harbor run -a oracle` == 1, `-a nop` == 0)
  G7  tools resolve in a login shell    (added after an agent pilot found `bash -lc` hid them)

G4 is the answer to the Epoch review of ALE (public answer / lax grader / grader tampering):
the shortcut has to be tried by us before an RL policy finds it.
"""
from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

from . import dockerlib
from .spec import Attack, Probe, TaskBuild


def tags(task_id: str) -> tuple[str, str]:
    return f"alegen-{task_id}-agent:local", f"alegen-{task_id}-verifier:local"


def _same_recipe_image(task_dir: Path) -> str | None:
    """Agent image of another task in the same directory whose environment/Dockerfile is byte-identical
    (same tool recipe; only the inputs differ), if that image exists on this host."""
    mine = (task_dir / "environment" / "Dockerfile").read_bytes()
    for sib in sorted(task_dir.parent.iterdir()):
        df = sib / "environment" / "Dockerfile"
        if sib == task_dir or not df.exists() or df.read_bytes() != mine:
            continue
        tag = tags(sib.name)[0]
        if subprocess.run(["docker", "image", "inspect", tag], capture_output=True).returncode == 0:
            return tag
    return None


def verify(verifier_tag: str, out_dir: Path, work: Path) -> tuple[float | None, dict]:
    """Run the verifier image on an output directory. Returns (reward or None on verifier error, details)."""
    name = "alegen-ver-" + next(tempfile._get_candidate_names())
    out_dir.mkdir(parents=True, exist_ok=True)
    subprocess.run(["docker", "create", "--name", name, "--network", "none", verifier_tag,
                    "bash", "/tests/test.sh"], check=True, capture_output=True)
    work.mkdir(parents=True, exist_ok=True)
    logs = work / "logs"
    try:
        # same transfer Harbor does for a separate verifier: artifacts land at their original path
        subprocess.run(["docker", "cp", f"{out_dir.resolve()}/.", f"{name}:/task/output"],
                       check=True, capture_output=True)
        r = subprocess.run(["docker", "start", "-a", name], capture_output=True, text=True, timeout=3600)
        if logs.exists():
            shutil.rmtree(logs)
        subprocess.run(["docker", "cp", f"{name}:/logs/verifier", str(logs)], capture_output=True)
    finally:
        subprocess.run(["docker", "rm", "-f", name], capture_output=True)
    reward_f = logs / "reward.txt"
    if not reward_f.exists():
        return None, {"verifier_error": (r.stdout + r.stderr)[-2000:]}
    details = json.loads((logs / "details.json").read_text()) if (logs / "details.json").exists() else {}
    return float(reward_f.read_text().strip()), details


def _materialize(files: dict, dst: Path) -> None:
    if dst.exists():
        shutil.rmtree(dst)
    dst.mkdir(parents=True)
    for rel, src in (files or {}).items():
        p = dst / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        if isinstance(src, Path):
            shutil.copyfile(src, p)
        elif isinstance(src, bytes):
            p.write_bytes(src)
        else:
            p.write_text(src)


def _run_script(agent_tag: str, script: str, task_dir: Path, out: Path, timeout: int,
                files: dict | None = None) -> subprocess.CompletedProcess:
    sdir = Path(tempfile.mkdtemp(prefix="alegen-script-"))
    _materialize(files or {}, sdir)
    (sdir / "run.sh").write_text(script)
    try:
        return dockerlib.run_in_image(
            agent_tag, "cd /task && bash /alegen_script/run.sh",
            mounts={sdir: "/alegen_script", task_dir / "solution": "/solution"},
            copy_out={"/task/output": out}, timeout=timeout)
    finally:
        shutil.rmtree(sdir, ignore_errors=True)


def generic_attacks(oracle_out: Path) -> list[Attack]:
    emptied = {}
    header_only = {}
    for p in sorted(oracle_out.rglob("*")):
        if p.is_file():
            rel = str(p.relative_to(oracle_out))
            emptied[rel] = b""
            try:
                lines = p.read_text().splitlines()
                # single-line files have no "header": leave them empty for this attack
                first = lines[0] + "\n" if len(lines) >= 3 else ""
            except Exception:
                first = ""
            header_only[rel] = first
    return [
        Attack("empty_output", "agent wrote nothing", 0.0, files={}),
        Attack("inputs_copied", "visible inputs copied into output/", 0.1,
               script="cp -r /task/input/. /task/output/ 2>/dev/null; true"),
        Attack("zero_byte_files", "every expected file exists but is empty", 0.0, files=emptied),
        Attack("header_only", "every expected text file holds only its first line", 0.1, files=header_only),
    ]


def harbor_check(task_dir: Path, agent_tag: str, verifier_tag: str, work: Path) -> dict:
    """Run the task through Harbor itself with the oracle and nop agents, using the local images."""
    tmp = work / "harbor_task" / task_dir.name
    if tmp.parent.exists():
        shutil.rmtree(tmp.parent)
    shutil.copytree(task_dir, tmp)
    toml = (tmp / "task.toml").read_text()
    toml = toml.replace("[verifier.environment]\n", f'[verifier.environment]\ndocker_image = "{verifier_tag}"\n')
    toml = toml.replace("[environment]\n", f'[environment]\ndocker_image = "{agent_tag}"\n')
    if os.environ.get("ALEGEN_NO_EGRESS_SIDECAR"):
        # sandbox kernels without nftables cannot run Harbor's allowlist sidecar; oracle/nop need no network
        toml = re.sub(r'(\[agent\]\ntimeout_sec = [0-9.]+\n)network_mode = "allowlist"\nallowed_hosts = .*\n',
                      r'\1network_mode = "no-network"\n', toml)
    (tmp / "task.toml").write_text(toml)
    out = {}
    for agent in ("oracle", "nop"):
        jobs = work / "harbor_jobs"
        name = f"{task_dir.name}-{agent}"
        if (jobs / name).exists():
            shutil.rmtree(jobs / name)
        env = dict(os.environ, PATH=os.path.expanduser("~/.local/bin") + ":" + os.environ["PATH"])
        r = subprocess.run(["harbor", "run", "-p", str(tmp), "-a", agent, "-o", str(jobs),
                            "--job-name", name, "--quiet"], capture_output=True, text=True, env=env,
                           timeout=4 * 3600)
        rewards = []
        for res in (jobs / name).rglob("result.json"):
            try:
                j = json.loads(res.read_text())
            except Exception:
                continue
            vr = j.get("verifier_result") if isinstance(j, dict) else None
            if isinstance(vr, dict) and vr.get("rewards") is not None:
                rw = vr["rewards"]
                rewards.append(rw.get("reward") if isinstance(rw, dict) else rw)
        out[agent] = {"rewards": rewards, "returncode": r.returncode,
                      "tail": (r.stdout + r.stderr)[-600:] if not rewards else ""}
    return out


def validate(b: TaskBuild, task_dir: Path, results_dir: Path, harbor: bool = True,
             keep_work: bool = False) -> dict:
    t0 = time.time()
    agent_tag, verifier_tag = tags(b.task_id)
    work = results_dir / "work" / b.task_id
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    rep: dict = {"task_id": b.task_id, "family": b.family, "seed": b.seed, "level": b.level,
                 "checks": [], "ok": True}

    def record(gate: str, name: str, ok: bool, **kw):
        rep["checks"].append({"gate": gate, "name": name, "ok": bool(ok), **kw})
        rep["ok"] = rep["ok"] and bool(ok)
        print(f"  [{'ok' if ok else 'FAIL'}] {gate} {name} " + " ".join(f"{k}={v}" for k, v in kw.items() if k != "details"))

    # G1
    try:
        base = _same_recipe_image(task_dir)
        if not (base and dockerlib.derive_agent_image(base, task_dir / "environment", agent_tag)):
            base = None
            dockerlib.build_image(task_dir / "environment", agent_tag)
        dockerlib.build_image(task_dir / "tests", verifier_tag)
        record("G1", "images_build", True, **({"agent_image_derived_from": base} if base else {}))
    except Exception as e:
        record("G1", "images_build", False, error=str(e)[-1500:])
        _save(rep, results_dir, b.task_id, t0)
        return rep

    # G7 the software the instruction promises "on PATH" must also resolve in a login shell
    probe = ("A=$(command -v python3); B=$(bash -lc 'command -v python3'); echo \"nonlogin=$A login=$B\"; "
             "bash -lc 'python3 -c \"import numpy\"' && echo NUMPY_OK")
    pr = dockerlib.run_in_image(agent_tag, probe, timeout=300)
    line = (pr.stdout or "").strip().splitlines()
    first = line[0] if line else ""
    same = first.startswith("nonlogin=") and first.split(" login=")[0].split("=", 1)[1] == first.split(" login=")[-1] != ""
    record("G7", "tools_in_login_shell", same and "NUMPY_OK" in (pr.stdout or ""), detail=first)

    # G2 oracle
    oracle_out = work / "oracle" / "output"
    t1 = time.time()
    r = _run_script(agent_tag, "bash /solution/solve.sh", task_dir, oracle_out, b.agent_timeout_sec)
    oracle_min = round((time.time() - t1) / 60, 2)
    score, det = verify(verifier_tag, oracle_out, work / "oracle")
    record("G2", "oracle", score is not None and score >= 0.99, score=score, runtime_min=oracle_min,
           failure=det.get("failure") or det.get("verifier_error", "")[-300:],
           stderr=(r.stderr or "")[-300:] if (score or 0) < 0.99 else "")
    rep["oracle_runtime_min"] = oracle_min
    rep["oracle_components"] = det.get("components")

    # G5 determinism
    score2, _ = verify(verifier_tag, oracle_out, work / "oracle2")
    record("G5", "verifier_deterministic", score == score2, first=score, second=score2)

    # G3 probes
    for p in b.probes:
        out = work / "probe" / p.name / "output"
        if p.script is not None:
            _run_script(agent_tag, p.script, task_dir, out, b.agent_timeout_sec, p.script_files)
        else:
            _materialize(p.files, out)
        s, d = verify(verifier_tag, out, work / "probe" / p.name)
        record("G3", "probe:" + p.name, s is not None and s >= p.min_score, score=s, floor=p.min_score,
               failure=d.get("failure") or d.get("verifier_error", "")[-300:])

    # G4 attacks
    if b.level >= 4 and int(b.extra_metadata.get("source_level", 4)) >= 4:   # the level-4 inputs must defeat the level-3 reference method
        record("G4", "has_l3_reference_attack", any(a.name == "l3_reference" for a in b.attacks))
    if b.level >= 7:   # the level-7 inputs must defeat the reference method of the level below
        record("G4", "has_l6_reference_attack", any(a.name == "l6_reference" for a in b.attacks))
    for a in generic_attacks(oracle_out) + list(b.attacks):
        out = work / "attack" / a.name / "output"
        if a.script is not None:
            _run_script(agent_tag, a.script, task_dir, out, b.agent_timeout_sec, a.script_files)
        else:
            _materialize(a.files, out)
        s, d = verify(verifier_tag, out, work / "attack" / a.name)
        # a verifier crash on hostile input is itself a defect
        record("G4", "attack:" + a.name, s is not None and s <= a.max_score, score=s, ceiling=a.max_score,
               failure=(d.get("failure") or d.get("verifier_error", ""))[-200:])

    # G6 harbor
    if harbor:
        h = harbor_check(task_dir, agent_tag, verifier_tag, work)
        o, n = h["oracle"]["rewards"], h["nop"]["rewards"]
        record("G6", "harbor_oracle", bool(o) and all(x is not None and x >= 0.99 for x in o), rewards=o,
               tail=h["oracle"]["tail"])
        record("G6", "harbor_nop", bool(n) and all(x == 0 for x in n), rewards=n, tail=h["nop"]["tail"])

    _save(rep, results_dir, b.task_id, t0)
    if not keep_work and rep["ok"]:
        shutil.rmtree(work, ignore_errors=True)
    return rep


def _save(rep: dict, results_dir: Path, task_id: str, t0: float) -> None:
    rep["validation_minutes"] = round((time.time() - t0) / 60, 2)
    d = results_dir / "validation"
    d.mkdir(parents=True, exist_ok=True)
    (d / f"{task_id}.json").write_text(json.dumps(rep, indent=1, default=str))
