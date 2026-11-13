"""Dockerfile building blocks shared by all families, plus local build/run helpers.

All agent images start from the same base lines so that Docker layers are shared between tasks.
The task Dockerfiles are plain and portable. `build_image` only adds a CA/proxy preamble to a
*temporary* copy when the build host sits behind a TLS-inspecting proxy (ALEGEN_CA_BUNDLE set or
/root/.ccr/ca-bundle.crt present); the files written into the task directory are never modified.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

MICROMAMBA_VERSION = "2.0.5"

AGENT_BASE = f"""\
FROM python:3.12-slim
ENV DEBIAN_FRONTEND=noninteractive PIP_NO_CACHE_DIR=1 PIP_DISABLE_PIP_VERSION_CHECK=1 \\
    MAMBA_ROOT_PREFIX=/opt/mamba LC_ALL=C.UTF-8
RUN apt-get update && apt-get install -y --no-install-recommends \\
      bzip2 ca-certificates curl file git less procps tmux unzip xz-utils \\
    && rm -rf /var/lib/apt/lists/*
RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/{MICROMAMBA_VERSION} \\
    | tar -xj -C /usr/local bin/micromamba
"""

AGENT_TAIL = """\
WORKDIR /task
COPY input/ /task/input/
COPY software/ /task/software/
RUN mkdir -p /task/output && chmod -R a-w /task/input && chmod -R a+rx /task/software
"""


def agent_dockerfile(conda: list[str], pip: list[str] | None = None, extra: str = "",
                     python: str = "3.12") -> str:
    """Standard agent image: one conda env at /opt/tools that owns python and every tool."""
    pkgs = " ".join([f"python={python}", "numpy", "scipy", "pandas"] + list(conda))
    out = AGENT_BASE
    out += (
        "RUN micromamba create -y -p /opt/tools -c conda-forge -c bioconda "
        f"{pkgs} \\\n    && micromamba clean -a -y\n"
        "ENV PATH=/opt/tools/bin:$PATH\n"
    )
    if pip:
        out += "RUN pip install " + " ".join(pip) + "\n"
    if extra:
        out += extra.rstrip() + "\n"
    out += AGENT_TAIL
    return out


def verifier_dockerfile(pip: list[str] | None = None, conda: list[str] | None = None,
                        extra: str = "") -> str:
    """Dedicated verifier image. It bundles tests/ (score.py + hidden reference); the agent image
    never contains them. Keep it small: python + numpy unless the grader truly needs a domain tool."""
    pip = pip if pip is not None else ["numpy==2.1.3"]
    if conda:
        out = AGENT_BASE
        pkgs = " ".join(["python=3.12", "numpy"] + list(conda))
        out += (
            "RUN micromamba create -y -p /opt/tools -c conda-forge -c bioconda "
            f"{pkgs} \\\n    && micromamba clean -a -y\n"
            "ENV PATH=/opt/tools/bin:$PATH\n"
        )
    else:
        out = "FROM python:3.12-slim\nENV PIP_NO_CACHE_DIR=1 PIP_DISABLE_PIP_VERSION_CHECK=1\n"
    if pip:
        out += "RUN pip install " + " ".join(pip) + "\n"
    if extra:
        out += extra.rstrip() + "\n"
    out += "WORKDIR /task\nCOPY . /tests/\n"
    return out


# --------------------------------------------------------------------------- local build / run

def _ca_bundle() -> Path | None:
    p = os.environ.get("ALEGEN_CA_BUNDLE") or "/root/.ccr/ca-bundle.crt"
    return Path(p) if Path(p).is_file() else None


_PREAMBLE = """\
COPY --from=alegen_ca {name} /etc/ssl/alegen-ca-bundle.crt
ENV SSL_CERT_FILE=/etc/ssl/alegen-ca-bundle.crt REQUESTS_CA_BUNDLE=/etc/ssl/alegen-ca-bundle.crt \\
    CURL_CA_BUNDLE=/etc/ssl/alegen-ca-bundle.crt PIP_CERT=/etc/ssl/alegen-ca-bundle.crt \\
    GIT_SSL_CAINFO=/etc/ssl/alegen-ca-bundle.crt
RUN mkdir -p /etc/apt/apt.conf.d \\
 && echo 'Acquire::https::CAInfo "/etc/ssl/alegen-ca-bundle.crt";' > /etc/apt/apt.conf.d/99alegen-ca \\
 && for f in /etc/apt/sources.list /etc/apt/sources.list.d/*; do \\
      [ -f "$f" ] && sed -i -E 's#http://(deb\\.debian\\.org|archive\\.ubuntu\\.com|security\\.ubuntu\\.com)#https://\\1#g' "$f"; \\
    done; true
"""


def build_image(context: Path, tag: str, dockerfile: Path | None = None, quiet: bool = True) -> None:
    """docker build. Behind a TLS-inspecting proxy, build from a temp Dockerfile with a CA preamble."""
    dockerfile = dockerfile or (context / "Dockerfile")
    text = dockerfile.read_text()
    cmd = ["docker", "build", "-t", tag]
    ca = _ca_bundle()
    tmp = None
    if ca is not None:
        pre = _PREAMBLE.format(name=ca.name)
        text = re.sub(r"(?m)^(FROM .*)$", lambda m: m.group(1) + "\n" + pre, text)
        tmp = tempfile.NamedTemporaryFile("w", suffix=".Dockerfile", delete=False)
        tmp.write(text)
        tmp.close()
        cmd += ["--network", "host", "--build-context", f"alegen_ca={ca.parent}", "-f", tmp.name]
        for k in ("HTTPS_PROXY", "https_proxy"):
            v = os.environ.get("HTTPS_PROXY") or os.environ.get("https_proxy")
            if v:
                cmd += ["--build-arg", f"{k}={v}"]
    else:
        cmd += ["-f", str(dockerfile)]
    cmd.append(str(context))
    try:
        r = subprocess.run(cmd, capture_output=quiet, text=True, errors="replace")
        if r.returncode != 0:
            tail = (r.stdout or "")[-3000:] + (r.stderr or "")[-6000:]
            raise RuntimeError(f"docker build failed for {tag}:\n{tail}")
    finally:
        if tmp:
            os.unlink(tmp.name)


_INPUT_TAIL = ("COPY input/ /task/input/\nCOPY software/ /task/software/\n"
               "RUN mkdir -p /task/output && chmod -R a-w /task/input && chmod -R a+rx /task/software\n")


def derive_agent_image(base_tag: str, context: Path, tag: str) -> bool:
    """Build-host shortcut: the agent image of a task whose environment/Dockerfile is byte-identical to that of
    a task whose image already exists here differs from it only in /task/input and /task/software. Build it on
    top of that image instead of installing the tools again (no network, a few MB instead of GBs).
    Returns False (nothing built) if the Dockerfile does not end with the standard input lines."""
    if not (context / "Dockerfile").read_text().endswith(_INPUT_TAIL):
        return False
    text = (f"FROM {base_tag}\nUSER root\n"
            "RUN chmod -R u+w /task/input /task/software && rm -rf /task/input /task/software\n" + _INPUT_TAIL)
    with tempfile.NamedTemporaryFile("w", suffix=".Dockerfile", delete=False) as tmp:
        tmp.write(text)
    try:
        r = subprocess.run(["docker", "build", "-t", tag, "-f", tmp.name, str(context)], capture_output=True, text=True,
                           errors="replace")
        if r.returncode != 0:
            raise RuntimeError(f"docker build (derived from {base_tag}) failed for {tag}:\n{(r.stderr or '')[-4000:]}")
    finally:
        os.unlink(tmp.name)
    return True


def run_in_image(tag: str, script: str, mounts: dict[Path, str] | None = None,
                 copy_out: dict[str, Path] | None = None, timeout: int = 7200,
                 env: dict[str, str] | None = None, network: str = "none",
                 cpus: int | None = None) -> subprocess.CompletedProcess:
    """Run a bash script in a fresh container (no network by default); optionally copy paths out."""
    name = "alegen-run-" + next(tempfile._get_candidate_names())
    cmd = ["docker", "create", "--name", name, "--network", network]
    if cpus:
        cmd += ["--cpus", str(cpus)]
    for host, cont in (mounts or {}).items():
        cmd += ["-v", f"{Path(host).resolve()}:{cont}:ro"]
    for k, v in (env or {}).items():
        cmd += ["-e", f"{k}={v}"]
    cmd += [tag, "bash", "-c", script]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    try:
        r = subprocess.run(["docker", "start", "-a", name], capture_output=True, text=True,
                           timeout=timeout)
        for cont, host in (copy_out or {}).items():
            host = Path(host)
            if host.exists():
                shutil.rmtree(host)
            host.parent.mkdir(parents=True, exist_ok=True)
            subprocess.run(["docker", "cp", f"{name}:{cont}", str(host)], check=False,
                           capture_output=True)
        return r
    finally:
        subprocess.run(["docker", "rm", "-f", name], capture_output=True)


def tool_image(conda: list[str], pip: list[str] | None = None, python: str = "3.12") -> str:
    """Build (once) a tools-only image a generator can use for tool-backed steps, e.g. to cross-check
    a planted truth with two independent programs. Same base layers as the agent images."""
    import hashlib
    key = hashlib.sha1((" ".join(conda) + "|" + " ".join(pip or []) + "|" + python).encode()).hexdigest()[:12]
    tag = f"alegen-tools-{key}:local"
    have = subprocess.run(["docker", "image", "inspect", tag], capture_output=True)
    if have.returncode == 0:
        return tag
    ctx = Path(tempfile.mkdtemp(prefix="alegen-tools-"))
    try:
        (ctx / "input").mkdir()
        (ctx / "software").mkdir()
        (ctx / "input" / ".keep").write_text("")
        (ctx / "software" / ".keep").write_text("")
        (ctx / "Dockerfile").write_text(agent_dockerfile(conda=conda, pip=pip, python=python))
        build_image(ctx, tag)
    finally:
        shutil.rmtree(ctx, ignore_errors=True)
    return tag
