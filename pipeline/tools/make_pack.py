"""The task pack: the tasks listed in pack/tasks.txt, with the pipeline, the reports and the evidence.

  cd pipeline
  ../.venv/bin/python tools/make_pack.py tables              # refresh the tables of pack/README.md from results/
  ../.venv/bin/python tools/make_pack.py build <out.zip>     # assemble the pack from the committed tree (HEAD)
  ../.venv/bin/python tools/make_pack.py build-parts <dir>   # the same pack as several smaller zips (unpack all into one place)

`build` takes every tracked file except the task directories that are not in the manifest; in the pack,
pack/README.md becomes README.md and the full report (README.md of the repository) becomes REPORT.md.
Numbers come from the same readers as docs/TASKS.md (alegen.g8); nothing is typed by hand.
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "pipeline"))
from alegen import g8  # noqa: E402

CFG = json.loads((ROOT / "results" / "model_runs.json").read_text())
MODELS = list(CFG["models"])
YARD = CFG["g8"]["yardstick"]
SHORT = {"gemini-3.8-flash": "Gemini 3.8 flash", "kimi-k3": "Kimi k3", "deepseek-v4.1-flash": "DeepSeek v4.1 flash",
         "qwen3.8-flash": "Qwen 3.8 flash", "glm-5.3-flash": "GLM 5.3 flash", "minimax-m3.1-flash": "MiniMax M3.1 flash"}
FAM = {"amplicon-16s": "amplicon", "md-system-repair": "md", "somatic-tn": "somatic", "crispr-screen": "crispr",
       "flow-compensation": "flow", "hcs-dose-response": "hcs", "rnaseq-confound": "rnaseq", "trio-liftover": "trio",
       "variant-review": "variant_review", "xlms-model-selection": "xlms"}


def manifest() -> list[str]:
    return [x.strip() for x in (ROOT / "pack" / "tasks.txt").read_text().splitlines() if x.strip() and not x.startswith("#")]


def label(task: str) -> str:
    m = re.match(r"(.+?)-s(\d+)-l(\d+)$", task)
    return f"{FAM[m.group(1)]} L{m.group(3)}" + ("" if m.group(2) == "0001" else f" 种子 {int(m.group(2))}")


def ordered(att: list[dict]) -> list[dict]:
    return sorted(att, key=lambda a: (bool(re.search(r"__r\d+", a["job"])), a["job"]))


def no_score(model: str, task: str) -> bool:
    for pat in CFG["models"][model]["jobs"]:
        for jd in (ROOT / "results" / "jobs").glob(pat):
            for r in jd.glob(f"{task}__*/result.json"):
                if (json.loads(r.read_text()).get("verifier_result") or {}).get("rewards") is None:
                    return True
    return False


def table_scores() -> list[str]:
    att = {m: g8.attempts(m) for m in MODELS}
    verdict = {r["task"]: r["verdict"] for r in g8.verdicts(YARD, CFG["g8"].get("min_attempts", 2))}
    zh = {"keep": "保留", "drop": "淘汰", "unstable": "不稳定", "pending": "仅一次"}
    others = [m for m in MODELS if m != YARD]
    out = ["| 题目 | 事前难度 | 难度系数 | " + SHORT[YARD] + "：每次尝试 | 难度筛选 | " + " | ".join(SHORT[m] for m in others) + " |",
           "|---|---|---|---|---|" + "---|" * len(others)]
    firsts, post = [], []
    for t in manifest():
        toml = (ROOT / "tasks" / t / "task.toml").read_text()
        prior = float(re.search(r"(?m)^difficulty_prior\s*=\s*([0-9.]+)", toml).group(1))
        a = ordered(att[YARD].get(t, []))
        std = [x for x in a if "__t300" not in x["job"]]
        long_ = [x for x in a if "__t300" in x["job"]]
        base = std or long_                      # the difficulty uses the attempts under the task's regular turn limit
        d = 1 - sum(x["reward"] for x in base) / len(base)
        post.append(d)
        firsts.append(a[0]["reward"])
        cell = "、".join(f"{x['reward']:.3f}" for x in std)
        if long_:
            cell += ("；" if std else "") + "、".join(f"{x['reward']:.3f}" for x in long_) + "（300 轮）"
        rest = []
        for m in others:
            b = ordered(att[m].get(t, []))
            rest.append(f"{b[0]['reward']:.3f}" if b else ("无分" if no_score(m, t) else "–"))
        out.append(f"| {label(t)} | {prior:.3f} | {d:.2f} | {cell} | {zh[verdict[t]]} | " + " | ".join(rest) + " |")
    n = len(firsts)
    out.append(f"| **{n} 道合计** | | **平均 {sum(post) / n:.2f}** | **第一次尝试：通过 {sum(x >= 0.9 for x in firsts)} / {n}，平均分 "
               f"{sum(firsts) / n:.3f}** | | " + " | ".join(
                   (lambda v: f"通过 {sum(x >= 0.9 for x in v)} / {len(v)}" if v else "–")(
                       [ordered(att[m][t])[0]["reward"] for t in manifest() if att[m].get(t)]) for m in others) + " |")
    return out


def table_gates() -> list[str]:
    out = ["| 题目 | 参考解法 | 独立第二解法 | 探针（全部满分） | 攻击数 | 最高攻击得分 | Harbor oracle / nop | 参考解法耗时 |",
           "|---|---|---|---|---|---|---|---|"]
    for t in manifest():
        v = json.loads((ROOT / "results" / "validation" / f"{t}.json").read_text())
        ch = {c["name"]: c for c in v["checks"]}
        atk = [c for c in v["checks"] if c["name"].startswith("attack:") and c.get("score") is not None]
        prb = [c for c in v["checks"] if c["name"].startswith("probe:")]
        top = max(atk, key=lambda c: c["score"])
        hb = "1.0 / 0.0" if ch.get("harbor_oracle", {}).get("ok") and ch.get("harbor_nop", {}).get("ok") else "未跑"
        assert v["ok"] and all(c.get("ok") for c in v["checks"]), t
        out.append(f"| {label(t)} | {ch['oracle']['score']:.3f} | {ch['probe:independent_solution']['score']:.3f} | {len(prb)} | "
                   f"{len(atk)} | {top['score']:.3f}（`{top['name'][7:]}`） | {hb} | {v.get('oracle_runtime_min', 0):.1f} 分钟 |")
    return out


def tables() -> int:
    f = ROOT / "pack" / "README.md"
    s = f.read_text()
    for k, lines in {"scores": table_scores(), "gates": table_gates()}.items():
        pat = re.compile(rf"(<!-- BEGIN:{k} -->\n).*?(<!-- END:{k} -->)", re.S)
        assert pat.search(s), k
        s = pat.sub(lambda m: m.group(1) + "\n".join(lines) + "\n" + m.group(2), s)
    f.write_text(s)
    return 0


def build(out: Path) -> int:
    keep = set(manifest())
    files = subprocess.run(["git", "ls-files", "-z"], cwd=ROOT, capture_output=True, check=True).stdout.decode().split("\0")
    dirty = subprocess.run(["git", "status", "--porcelain"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    if dirty:
        print("working tree is not clean; the pack is built from the files on disk that git tracks:\n" + dirty[:400])
    n = 0
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for rel in sorted(x for x in files if x):
            parts = rel.split("/")
            if parts[0] == "tasks" and parts[1] not in keep:
                continue
            name = {"README.md": "REPORT.md", "pack/README.md": "README.md"}.get(rel, rel)
            z.write(ROOT / rel, "ale-lifesci-pack/" + name)
            n += 1
    tasks = {x.split("/")[2] for x in zipfile.ZipFile(out).namelist() if x.startswith("ale-lifesci-pack/tasks/")}
    assert tasks == keep, (tasks ^ keep)
    print(f"{out}: {n} files, {len(tasks)} tasks, {out.stat().st_size / 2**20:.1f} MiB")
    return 0


def build_parts(outdir: Path) -> int:
    """The same pack as several smaller zips that unpack into one directory (for downloads that fail on one big file):
    core (reports, pipeline, docs, validation records), evidence (results/jobs), and one zip per task family."""
    import hashlib
    keep = set(manifest())
    files = subprocess.run(["git", "ls-files", "-z"], cwd=ROOT, capture_output=True, check=True).stdout.decode().split("\0")
    parts: dict[str, list[tuple[str, str]]] = {}
    for rel in sorted(x for x in files if x):
        seg = rel.split("/")
        if seg[0] == "tasks":
            if seg[1] not in keep:
                continue
            part = "tasks-" + FAM[re.match(r"(.+?)-s\d+-l\d+$", seg[1]).group(1)]
        elif rel.startswith("results/jobs/"):
            part = "evidence"
        else:
            part = "core"
        name = {"README.md": "REPORT.md", "pack/README.md": "README.md"}.get(rel, rel)
        parts.setdefault(part, []).append((rel, name))
    outdir.mkdir(parents=True, exist_ok=True)
    order = ["core", "evidence"] + sorted(k for k in parts if k.startswith("tasks-"))
    sums = []
    for i, part in enumerate(order, 1):
        out = outdir / f"ale-lifesci-pack15-part{i}of{len(order)}-{part}.zip"
        with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as z:
            for rel, name in parts[part]:
                z.write(ROOT / rel, "ale-lifesci-pack/" + name)
        assert zipfile.ZipFile(out).testzip() is None
        h = hashlib.sha256(out.read_bytes()).hexdigest()
        sums.append(f"{h}  {out.name}")
        print(f"{out.name}: {len(parts[part])} files, {out.stat().st_size / 2**20:.1f} MiB")
    (outdir / "SHA256SUMS.txt").write_text("\n".join(sums) + "\n")
    return 0


if __name__ == "__main__":
    if len(sys.argv) >= 2 and sys.argv[1] == "tables":
        sys.exit(tables())
    if len(sys.argv) == 3 and sys.argv[1] == "build":
        sys.exit(build(Path(sys.argv[2])))
    if len(sys.argv) == 3 and sys.argv[1] == "build-parts":
        sys.exit(build_parts(Path(sys.argv[2])))
    print(__doc__)
    sys.exit(2)
