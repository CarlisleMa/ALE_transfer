"""Write a TaskBuild to disk as a Harbor task directory (schema 1.4, separate verifier).

Layout (mirrors ALE's input/ software/ output/ reference/ contract):
    <task>/instruction.md
    <task>/task.toml
    <task>/README.md                      task card for humans (truth provenance, judge rationale)
    <task>/environment/Dockerfile         agent image: tools + /task/input + /task/software
    <task>/environment/input/...          visible inputs (read-only in the container)
    <task>/environment/software/...       wrappers
    <task>/solution/solve.sh              oracle (SOLUTION=<n> selects an independent alternative)
    <task>/tests/Dockerfile               verifier image, bundles score.py + hidden reference/
    <task>/tests/test.sh, score.py, _scoring.py, reference/...
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path

from .spec import FileSrc, TaskBuild

# Model API endpoints an in-container agent CLI may need. Everything else is blocked while the
# agent runs, so "look the answer up" is not an available strategy (ALE only asked politely).
LLM_API_HOSTS = [
    "generativelanguage.googleapis.com", "aiplatform.googleapis.com", "oauth2.googleapis.com",
    "api.anthropic.com", "api.openai.com", "openrouter.ai", "api.deepseek.com",
    "dashscope.aliyuncs.com", "api.moonshot.cn", "open.bigmodel.cn",
]

TEST_SH = """#!/bin/bash
# Verifier entrypoint (runs in the dedicated verifier image; the agent never sees /tests).
mkdir -p /logs/verifier
python3 /tests/score.py --output /task/output --reference /tests/reference --logs /logs/verifier
"""


def difficulty_prior(b: TaskBuild) -> float:
    """A-priori difficulty in [0,1] from structural features of the reference solution."""
    steps = min(b.oracle_steps / 12.0, 1.0)
    tools = min(len(b.software) / 5.0, 1.0)
    traps = min(len(b.traps) / 3.0, 1.0)
    runtime = min(b.oracle_runtime_min / 30.0, 1.0)
    gui = 1.0 if b.needs_gui else 0.0
    return round(0.15 * steps + 0.15 * tools + 0.40 * traps + 0.15 * runtime + 0.15 * gui, 3)


def prior_tier(p: float) -> str:
    return "last-exam" if p >= 0.70 else ("full-spectrum" if p >= 0.45 else "near-term")


def _write(dst: Path, src: FileSrc, executable: bool = False) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if isinstance(src, Path):
        shutil.copyfile(src, dst)
    elif isinstance(src, bytes):
        dst.write_bytes(src)
    else:
        dst.write_text(src)
    if executable or dst.suffix == ".sh":
        dst.chmod(0o755)


_PATH_ENV = "ENV PATH=/opt/tools/bin:$PATH\n"
_PATH_PROFILE = ("# login shells (tmux, `bash -l`) reset PATH in /etc/profile: keep the task's tools visible there too\n"
                 "RUN printf 'export PATH=/opt/tools/bin:$PATH\\n' > /etc/profile.d/00-task-tools.sh\n")


def _login_shell_path(dockerfile: str) -> str:
    """Found by the first agent pilot: `bash -lc` hid every tool of the task, so 'available on PATH'
    was false for agents whose shell is a login shell. Applied to every family at packaging time."""
    if _PATH_ENV in dockerfile and "/etc/profile.d/" not in dockerfile:
        return dockerfile.replace(_PATH_ENV, _PATH_ENV + _PATH_PROFILE, 1)
    return dockerfile


def _toml_str(s: str) -> str:
    return json.dumps(s, ensure_ascii=False)


def _toml_list(xs) -> str:
    return "[" + ", ".join(_toml_str(str(x)) for x in xs) + "]"


def task_toml(b: TaskBuild) -> str:
    prior = difficulty_prior(b)
    md = {
        "benchmark_alignment": "Agents' Last Exam (ALE) / life_sciences",
        "category": "life_sciences",
        "subdomain": b.subdomain,
        "family": b.family,
        "difficulty_tier_prior": prior_tier(prior),
        "difficulty_explanation": (
            f"{b.oracle_steps} tool invocations in the reference solution, {len(b.software)} tools, "
            f"{len(b.traps)} undisclosed decision point(s), ~{b.oracle_runtime_min:g} min oracle runtime"
            + (", GUI software" if b.needs_gui else "")),
    }
    lines = [
        'schema_version = "1.4"',
        'artifacts = ["/task/output"]',
        "",
        "[task]",
        f"name = {_toml_str('ale-ls-gym/' + b.task_id)}",
        'version = "1.0.0"',
        f"description = {_toml_str(b.summary)}",
        'authors = [{ name = "alegen pipeline" }]',
        f"keywords = {_toml_list(['life-sciences', b.family] + b.software)}",
        "",
        "[metadata]",
    ]
    for k, v in md.items():
        lines.append(f"{k} = {_toml_str(v)}")
    lines += [
        f"seed = {b.seed}",
        f"level = {b.level}",
        f"difficulty_prior = {prior}",
        f"pass_threshold = {b.pass_threshold}",
        f"oracle_steps = {b.oracle_steps}",
        f"oracle_runtime_min = {b.oracle_runtime_min}",
        f"needs_gui = {'true' if b.needs_gui else 'false'}",
        f"software = {_toml_list(b.software)}",
        f"capabilities = {_toml_list(b.capabilities)}",
        f"traps = {_toml_list(t.name for t in b.traps)}",
    ]
    for k, v in b.extra_metadata.items():
        # JSON of scalars and flat lists is valid TOML; JSON of a dict is not, and Harbor then
        # fails much later with an unrelated message ("Either datasets or tasks must be provided")
        if isinstance(v, dict) or (isinstance(v, (list, tuple)) and any(isinstance(x, (dict, list, tuple)) for x in v)):
            raise ValueError(f"extra_metadata[{k!r}] must be a scalar or a flat list (flatten nested values)")
        lines.append(f"{k} = {json.dumps(v, ensure_ascii=False)}")
    lines += [
        "",
        "[verifier]",
        f"timeout_sec = {float(b.verifier_timeout_sec)}",
        'environment_mode = "separate"',
        "",
        "[verifier.environment]",
        'network_mode = "no-network"',
        "",
        "[agent]",
        f"timeout_sec = {float(b.agent_timeout_sec)}",
        'network_mode = "allowlist"',
        f"allowed_hosts = {_toml_list(LLM_API_HOSTS)}",
        "",
        "[environment]",
        'network_mode = "public"',
        "build_timeout_sec = 3600.0",
        'os = "linux"',
        f"cpus = {b.cpus}",
        f"memory_mb = {b.memory_mb}",
        f"storage_mb = {b.storage_mb}",
        "",
    ]
    for m in b.mcp_servers:   # GUI-as-tool: desktop actions exposed as ordinary MCP tools (ALE's design)
        lines += ["[[environment.mcp_servers]]", f"name = {_toml_str(m['name'])}",
                  f"transport = {_toml_str(m.get('transport', 'stdio'))}"]
        if m.get("command"):
            lines.append(f"command = {_toml_str(m['command'])}")
        if m.get("args"):
            lines.append(f"args = {_toml_list(m['args'])}")
        if m.get("url"):
            lines.append(f"url = {_toml_str(m['url'])}")
        lines.append("")
    return "\n".join(lines)


def task_card(b: TaskBuild) -> str:
    prior = difficulty_prior(b)
    out = [f"# {b.title}", "", b.summary, "",
           f"- family: `{b.family}`  seed: `{b.seed}`  level: `{b.level}`",
           f"- ALE subdomain: {b.subdomain}",
           f"- software: {', '.join(b.software)}",
           f"- a-priori difficulty: {prior} ({prior_tier(prior)})  |  pass threshold: {b.pass_threshold}",
           f"- capabilities trained: {', '.join(b.capabilities)}",
           "", "## Undisclosed decision points (not shown to the agent)", ""]
    if not b.traps:
        out.append("none (level 1)")
    for t in b.traps:
        out += [f"- **{t.name}** - {t.what}",
                f"  - discoverable by: {t.discoverable_by}",
                f"  - trains against ALE failure mode: {t.ale_failure_mode}"]
    out += ["", "## Judge", "", b.judge_explanation_md.strip(), "",
            "## Verifier attack suite (must score low)", ""]
    for a in b.attacks:
        out.append(f"- `{a.name}` (max {a.max_score}): {a.description}")
    out += ["", "## Must-pass probes (guard against false negatives)", ""]
    for p in b.probes:
        out.append(f"- `{p.name}` (min {p.min_score}): {p.description}")
    return "\n".join(out) + "\n"


def write_task(b: TaskBuild, out_root: Path, scoring_src: Path | None = None) -> Path:
    d = Path(out_root) / b.task_id
    if d.exists():
        shutil.rmtree(d)
    (d / "environment" / "input").mkdir(parents=True)
    (d / "environment" / "software").mkdir(parents=True)
    _write(d / "instruction.md", b.instruction_md.strip() + "\n")
    _write(d / "task.toml", task_toml(b))
    _write(d / "README.md", task_card(b))
    _write(d / "environment" / "Dockerfile", _login_shell_path(b.agent_dockerfile))
    for rel, src in b.input_files.items():
        _write(d / "environment" / "input" / rel, src)
    for rel, src in b.software_files.items():
        _write(d / "environment" / "software" / rel, src, executable=True)
    if not b.software_files:
        _write(d / "environment" / "software" / "README.txt",
               "All tools are on PATH (conda env at /opt/tools).\n")
    if "solve.sh" not in b.solution_files:
        raise ValueError("solution_files must contain solve.sh")
    for rel, src in b.solution_files.items():
        _write(d / "solution" / rel, src)
    _write(d / "tests" / "Dockerfile", b.verifier_dockerfile)
    _write(d / "tests" / "test.sh", TEST_SH)
    scoring_src = scoring_src or (Path(__file__).parent / "scoring.py")
    _write(d / "tests" / "_scoring.py", scoring_src.read_text())
    if "score.py" not in b.tests_files:
        raise ValueError("tests_files must contain score.py")
    for rel, src in b.tests_files.items():
        _write(d / "tests" / rel, src)
    # battery definition for the validator (lives outside the task dir on purpose)
    return d
