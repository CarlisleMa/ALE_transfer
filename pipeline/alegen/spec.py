"""Data model shared by all task families."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Union

FileSrc = Union[bytes, str, Path]  # bytes/str = literal content, Path = copy from disk


@dataclass
class Trap:
    """An expert decision point hidden in the data (never named in the instruction)."""
    name: str            # short id, e.g. "genome_build_mismatch"
    what: str            # what is wrong / what the naive path gets wrong
    discoverable_by: str  # the standard QC step that reveals it (fairness proof)
    ale_failure_mode: str  # ALE failure taxonomy bucket this trains against


@dataclass
class Attack:
    """A shortcut submission the verifier must reject (score <= max_score)."""
    name: str
    description: str
    max_score: float = 0.2
    files: dict[str, FileSrc] | None = None   # literal files to place under /task/output
    script: str | None = None                 # or: bash script run inside the agent image
    script_files: dict[str, FileSrc] | None = None  # helper files mounted at /alegen_script


@dataclass
class Probe:
    """A submission that MUST score >= min_score (guards against false negatives)."""
    name: str
    description: str
    min_score: float = 0.99
    files: dict[str, FileSrc] | None = None
    script: str | None = None                 # e.g. an independent second solution
    script_files: dict[str, FileSrc] | None = None


@dataclass
class TaskBuild:
    task_id: str                 # directory name, e.g. "crispr-screen-s0001"
    family: str
    seed: int
    level: int                   # 1 = no trap, 2 = some traps, 3 = all traps, 4 = level 3 + depth elements (AUTHORING.md)
    title: str
    summary: str
    subdomain: str               # ALE taxonomy subdomain (3.1 biomolecular / 3.2 genomics / 3.3 cell_imaging)
    software: list[str]
    instruction_md: str
    agent_dockerfile: str        # environment/Dockerfile
    verifier_dockerfile: str     # tests/Dockerfile
    input_files: dict[str, FileSrc] = field(default_factory=dict)      # -> environment/input/
    software_files: dict[str, FileSrc] = field(default_factory=dict)   # -> environment/software/
    solution_files: dict[str, FileSrc] = field(default_factory=dict)   # -> solution/ (must hold solve.sh)
    tests_files: dict[str, FileSrc] = field(default_factory=dict)      # -> tests/ (score.py, reference/...)
    traps: list[Trap] = field(default_factory=list)
    attacks: list[Attack] = field(default_factory=list)
    probes: list[Probe] = field(default_factory=list)
    capabilities: list[str] = field(default_factory=list)  # ALE capability tags trained by this task
    judge_explanation_md: str = ""   # why the truth is right, gates, tolerances (goes to docs, not to the agent)
    pass_threshold: float = 0.9
    agent_timeout_sec: int = 7200
    verifier_timeout_sec: int = 900
    cpus: int = 4
    memory_mb: int = 16000
    storage_mb: int = 20000
    oracle_steps: int = 0        # number of distinct tool invocations in the reference solution
    oracle_runtime_min: float = 0.0
    needs_gui: bool = False
    extra_metadata: dict = field(default_factory=dict)
    mcp_servers: list[dict] = field(default_factory=list)  # e.g. the desktop (GUI-as-tool) server
