"""Development harness (not shipped): print the raw metric values of every battery entry of a validated task
from results/work/<task_id>/ (needs `make --keep-work`), as a markdown table row per entry.

  ../.venv/bin/python -m alegen.families.somatic_tn.dev.harvest somatic-tn-s0001-l3 [...]
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[5] / "results" / "work"


def main():
    for task in sys.argv[1:]:
        print(f"\n### {task}\n")
        print("| submission | pairing wrong | min precision | min recall | worst AF err | TP/FP/FN per patient | score |")
        print("|---|---|---|---|---|---|---|")
        entries = [("oracle", ROOT / task / "oracle")]
        for kind in ("probe", "attack"):
            d = ROOT / task / kind
            if d.is_dir():
                entries += [(f"{kind}:{p.name}", p) for p in sorted(d.iterdir())]
        for name, p in entries:
            f = p / "logs" / "details.json"
            if not f.exists():
                print(f"| {name} | (no details) | | | | | |")
                continue
            d = json.loads(f.read_text())
            i = d.get("info", {})
            per = " ".join(f"{v['tp']}/{v['fp']}/{v['fn']}" for k, v in sorted(i.items()) if isinstance(v, dict))
            print(f"| {name} | {len(i.get('pairing_wrong', []) or [])} | {i.get('min_precision', '-')} | "
                  f"{i.get('min_recall', '-')} | {i.get('worst_af_median_abs_error', '-')} | {per or '-'} | "
                  f"{d.get('score')}{' (' + d['failure'][:40] + ')' if d.get('failure') else ''} |")


if __name__ == "__main__":
    main()
