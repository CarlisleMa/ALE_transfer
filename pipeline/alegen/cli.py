"""CLI:  python -m alegen make <family> --seed 1 --level 3      generate + validate one task
         python -m alegen make <family> --seed 1 --level 5      level 4 without scaffolding (descaffold.py)
         python -m alegen list                                   list families
"""
from __future__ import annotations

import argparse
import importlib
import json
import pkgutil
import sys
from pathlib import Path

from . import descaffold, families, harbor, validate

ROOT = Path(__file__).resolve().parents[2]


def list_families() -> list[str]:
    return sorted(m.name for m in pkgutil.iter_modules(families.__path__) if not m.name.startswith("__"))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="alegen")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("list")
    mk = sub.add_parser("make")
    mk.add_argument("family")
    mk.add_argument("--seed", type=int, default=1)
    mk.add_argument("--level", type=int, default=3)
    mk.add_argument("--out", default=str(ROOT / "tasks"))
    mk.add_argument("--results", default=str(ROOT / "results"))
    mk.add_argument("--no-validate", action="store_true")
    mk.add_argument("--no-harbor", action="store_true")
    mk.add_argument("--keep-work", action="store_true")
    a = ap.parse_args(argv)
    if a.cmd == "list":
        print("\n".join(list_families()))
        return 0
    mod = importlib.import_module(f"alegen.families.{a.family}")
    work = Path(a.results) / "gen" / f"{a.family}-s{a.seed:04d}-l{a.level}"
    work.mkdir(parents=True, exist_ok=True)
    if a.level in (5, 6):   # the level-4 build of the same seed with the scaffolding (6: and the incident reports) removed
        b = descaffold.apply(mod.generate(seed=a.seed, level=descaffold.TOP_LEVEL[a.family], workdir=work), a.level)
    elif a.level == 7:      # a family's level-7 build (new silent decision points, level-4 texts) made pointer-free
        b = descaffold.apply(mod.generate(seed=a.seed, level=7, workdir=work), 7)
    else:
        b = mod.generate(seed=a.seed, level=a.level, workdir=work)
    d = harbor.write_task(b, Path(a.out))
    print(f"wrote {d}  (prior difficulty {harbor.difficulty_prior(b)})")
    if a.no_validate:
        return 0
    rep = validate.validate(b, d, Path(a.results), harbor=not a.no_harbor, keep_work=a.keep_work)
    print(json.dumps({"task": b.task_id, "ok": rep["ok"], "minutes": rep["validation_minutes"]}))
    return 0 if rep["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
