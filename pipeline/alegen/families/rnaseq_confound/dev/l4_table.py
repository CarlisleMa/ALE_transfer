"""Markdown table of the raw metrics of a level-4 `make --keep-work` run.
usage: python l4_table.py results/work/<task_id>
"""
import json, sys
from pathlib import Path
root = Path(sys.argv[1])
rows = [("oracle", root / "oracle")] + [(f"probe:{p.name}", p) for p in sorted((root / "probe").iterdir())] \
     + [(f"attack:{p.name}", p) for p in sorted((root / "attack").iterdir())]
print("| submission | score | agree | sex | labels | TP/FP | sens | FDP | err | audit MCC (missed, false alarms) | erratum MCC |")
print("|---|---|---|---|---|---|---|---|---|---|---|")
for name, d in rows:
    f = d / "logs" / "details.json"
    if not f.exists():
        print(f"| {name} | no details | | | | | | | | | |"); continue
    j = json.loads(f.read_text())
    if j.get("failure"):
        print(f"| {name} | {j['score']:.3f} | gate: {j['failure'][:70]} | | | | | | | | |"); continue
    i = j["info"]
    print(f"| {name} | {j['score']:.3f} | {i['count_agreement_min']:.3f} | {i['sex_correct']} | {i['labels_valid']} | {i['de_tp']}/{i['de_fp']} | "
          f"{i['de_sensitivity']:.3f} | {i['de_fdp']:.3f} | {i['lfc_median_abs_error']:.3f} | {i['audit_mcc']:.2f} ({i['audit_missed']}, {i['audit_false_alarms']}) | {i['erratum_mcc']:.2f} |")
