"""Gate-and-score helpers. This file is copied verbatim into every task's tests/ directory as
`_scoring.py`, so it must stay dependency-free (stdlib only).

Contract with Harbor:
  * candidate problems (missing/malformed output)  -> reward 0, reason recorded
  * verifier bugs (unexpected exception)           -> NO reward file, exit 1 (Harbor reports a
    verifier error instead of a fake zero; an RL pipeline should drop that sample)
Rules for score.py authors:
  * never execute, import, unpickle or `eval` anything the agent wrote
  * cap file sizes and row counts before parsing
  * every numeric comparison goes through `close`/`ramp` with a tolerance justified in the docs
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import math
import sys
import traceback
from pathlib import Path


class CandidateError(Exception):
    """The submission violates the output contract. Score 0."""


def ramp(x: float, lo: float, hi: float) -> float:
    """Piecewise-linear credit: 0 at/below lo, 1 at/above hi (lo < hi), or reversed if lo > hi."""
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return 0.0
    if lo == hi:
        return 1.0 if x >= hi else 0.0
    t = (x - lo) / (hi - lo)
    return max(0.0, min(1.0, t))


def close(a: float, b: float, rel: float = 0.0, abs_: float = 0.0) -> bool:
    if a is None or b is None:
        return False
    if not (math.isfinite(a) and math.isfinite(b)):
        return False
    return abs(a - b) <= max(abs_, rel * abs(b))


def prf(tp: int, fp: int, fn: int, beta: float = 1.0) -> tuple[float, float, float]:
    p = tp / (tp + fp) if tp + fp else 0.0
    r = tp / (tp + fn) if tp + fn else 0.0
    b2 = beta * beta
    f = (1 + b2) * p * r / (b2 * p + r) if p + r else 0.0
    return p, r, f


def to_float(s, what: str = "value") -> float:
    try:
        v = float(str(s).strip())
    except Exception:
        raise CandidateError(f"{what} is not numeric: {str(s)[:40]!r}")
    if not math.isfinite(v):
        raise CandidateError(f"{what} is not finite")
    return v


class Ctx:
    def __init__(self, output: Path, reference: Path):
        self.output = output
        self.reference = reference
        self.components: list[dict] = []
        self.gates: list[dict] = []
        self.info: dict = {}

    # ---- candidate files -------------------------------------------------------------
    def path(self, rel: str, max_bytes: int = 200_000_000) -> Path:
        p = (self.output / rel)
        try:
            rp = p.resolve()
            rp.relative_to(self.output.resolve())
        except Exception:
            raise CandidateError(f"{rel}: path escapes the output directory")
        if not rp.is_file():
            raise CandidateError(f"required output file missing: {rel}")
        size = rp.stat().st_size
        if size == 0:
            raise CandidateError(f"{rel} is empty")
        if size > max_bytes:
            raise CandidateError(f"{rel} is larger than {max_bytes} bytes")
        return rp

    def text(self, rel: str, max_bytes: int = 50_000_000) -> str:
        p = self.path(rel, max_bytes)
        raw = p.read_bytes()
        if rel.endswith(".gz") or raw[:2] == b"\x1f\x8b":
            try:
                raw = gzip.decompress(raw)
            except Exception:
                raise CandidateError(f"{rel} is not valid gzip")
            if len(raw) > 20 * max_bytes:
                raise CandidateError(f"{rel} decompresses to too much data")
        try:
            return raw.decode("utf-8-sig")
        except UnicodeDecodeError:
            raise CandidateError(f"{rel} is not UTF-8 text")

    def table(self, rel: str, required: list[str], delimiter: str = "\t",
              max_rows: int = 2_000_000, exact_header: bool = False) -> list[dict]:
        txt = self.text(rel)
        try:
            reader = csv.DictReader(io.StringIO(txt), delimiter=delimiter)
            header = reader.fieldnames or []
            rows = []
            for i, row in enumerate(reader):
                if i >= max_rows:
                    raise CandidateError(f"{rel} has more than {max_rows} rows")
                rows.append(row)
        except csv.Error as e:
            raise CandidateError(f"{rel} is not valid delimited text: {e}")
        header = [h.strip() for h in header]
        missing = [c for c in required if c not in header]
        if missing:
            raise CandidateError(f"{rel} is missing required column(s) {missing}; header={header[:20]}")
        if exact_header and header != required:
            raise CandidateError(f"{rel} header must be exactly {required}")
        out = []
        for row in rows:
            out.append({(k or "").strip(): (v.strip() if isinstance(v, str) else v) for k, v in row.items()})
        return out

    def json(self, rel: str, max_bytes: int = 5_000_000):
        try:
            return json.loads(self.text(rel, max_bytes))
        except json.JSONDecodeError as e:
            raise CandidateError(f"{rel} is not valid JSON: {e}")

    # ---- reference files (evaluator-controlled: failure here is a verifier bug) -------
    def ref_json(self, rel: str):
        return json.loads((self.reference / rel).read_text())

    def ref_table(self, rel: str, delimiter: str = "\t") -> list[dict]:
        raw = (self.reference / rel).read_bytes()
        if rel.endswith(".gz"):
            raw = gzip.decompress(raw)
        return list(csv.DictReader(io.StringIO(raw.decode()), delimiter=delimiter))

    # ---- scoring ---------------------------------------------------------------------
    def gate(self, ok: bool, reason: str) -> None:
        """Hard precondition. Failure forces the task score to 0."""
        self.gates.append({"ok": bool(ok), "reason": reason})
        if not ok:
            raise CandidateError("gate failed: " + reason)

    def component(self, name: str, value: float, weight: float, note: str = "") -> float:
        value = max(0.0, min(1.0, float(value)))
        self.components.append({"name": name, "value": round(value, 6), "weight": weight, "note": note})
        return value

    def total(self) -> float:
        w = sum(c["weight"] for c in self.components)
        if not self.components or abs(w - 1.0) > 1e-6:
            raise RuntimeError(f"component weights must sum to 1 (got {w})")
        return sum(c["value"] * c["weight"] for c in self.components)


def run(main, pass_threshold: float = 0.9) -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", default="/task/output")
    ap.add_argument("--reference", default="/tests/reference")
    ap.add_argument("--logs", default="/logs/verifier")
    a = ap.parse_args()
    logs = Path(a.logs)
    logs.mkdir(parents=True, exist_ok=True)
    ctx = Ctx(Path(a.output), Path(a.reference))
    details: dict = {"pass_threshold": pass_threshold}
    try:
        try:
            main(ctx)
            score = ctx.total()
            details["failure"] = None
        except CandidateError as e:
            score = 0.0
            details["failure"] = str(e)
        score = max(0.0, min(1.0, float(score)))
        details.update({"score": round(score, 6), "passed": score >= pass_threshold,
                        "gates": ctx.gates, "components": ctx.components, "info": ctx.info})
        (logs / "details.json").write_text(json.dumps(details, indent=1, default=str))
        (logs / "reward.txt").write_text(f"{score:.6f}\n")
        print(json.dumps({"score": details["score"], "passed": details["passed"],
                          "failure": details["failure"]}))
    except Exception:
        # verifier bug: leave no reward so the harness reports an error instead of a zero
        (logs / "verifier_error.txt").write_text(traceback.format_exc())
        traceback.print_exc()
        sys.exit(1)
