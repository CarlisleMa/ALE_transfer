"""Measured table for the docs from results/work/<task>/**/logs/details.json and results/validation/<task>.json."""
import json, sys
from pathlib import Path
R = Path("/home/user/ale-lifesci-gym/results")
rows = {}
for task in sys.argv[1:]:
    rep = json.loads((R / "validation" / f"{task}.json").read_text())
    for c in rep["checks"]:
        name = c["name"]
        if name == "oracle":
            d = R / "work" / task / "oracle" / "logs" / "details.json"
        elif name.startswith("probe:"):
            d = R / "work" / task / "probe" / name[6:] / "logs" / "details.json"
        elif name.startswith("attack:"):
            d = R / "work" / task / "attack" / name[7:] / "logs" / "details.json"
        else:
            continue
        info = json.loads(d.read_text()).get("info", {}) if d.exists() else {}
        rows.setdefault(name, {})[task] = (c.get("score"), info.get("mcc"), info.get("failure_tag_accuracy"),
                                           info.get("somatic_with_correct_vaf"),
                                           (info.get("session") or {}), (info.get("snapshots") or {}).get("valid_and_unique"), c["ok"])
    print(task, "ok" if rep["ok"] else "FAILED", "minutes", rep["validation_minutes"], "oracle_min", rep.get("oracle_runtime_min"))
tasks = sys.argv[1:]
print("| submission | " + " | ".join(f"{t[-8:]}: MCC / tag acc / VAF ok / images / score" for t in tasks) + " |")
print("|---|" + "---|" * len(tasks))
for name, per in rows.items():
    cells = []
    for t in tasks:
        v = per.get(t)
        if v is None:
            cells.append("n/a"); continue
        sc, m, ta, vf, ses, snap, ok = v
        if m is None:
            cells.append(f"- / - / - / - / {sc}")
        else:
            cells.append(f"{m} / {ta} / {vf} / {snap} / {sc}" + ("" if ok else " FAIL"))
    print(f"| `{name}` | " + " | ".join(cells) + " |")
