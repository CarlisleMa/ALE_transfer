"""Development harness (not shipped): score the file-based probes/attacks of a generated instance on the host.

  ../.venv/bin/python -m alegen.families.somatic_tn.dev.score_files <seed> <level>
"""
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

from alegen import harbor, validate
from alegen.families import somatic_tn


def main():
    seed, level = int(sys.argv[1]), int(sys.argv[2])
    root = Path(os.environ.get("SOMATIC_TN_SCRATCH", "/tmp/somatic_tn_dev")) / f"files-s{seed}-l{level}"
    b = somatic_tn.generate(seed, level, root / "gen")
    d = harbor.write_task(b, root / "tasks")
    for e in list(b.probes) + list(b.attacks):
        if e.files is None:
            continue
        out = Path(tempfile.mkdtemp(dir=root))
        validate._materialize(e.files, out / "output")
        subprocess.run([sys.executable, str(d / "tests" / "score.py"), "--output", str(out / "output"), "--reference",
                        str(d / "tests" / "reference"), "--logs", str(out / "logs")], capture_output=True,
                       env={"PYTHONPATH": str(d / "tests")})
        det = json.loads((out / "logs" / "details.json").read_text())
        lim = getattr(e, "min_score", None)
        lim = f"floor {lim}" if lim is not None else f"ceiling {e.max_score}"
        print(f"{e.name:34s} score={det['score']:.4f} ({lim}) minP={det['info'].get('min_precision')} "
              f"minR={det['info'].get('min_recall')} af={det['info'].get('worst_af_median_abs_error')} fail={det['failure']}")
        print("     ", {c["name"]: c["value"] for c in det["components"]})
    t = json.loads((d / "tests" / "reference" / "truth.json").read_text())
    for p, v in t["patients"].items():
        kinds = {}
        for r in v["scored"]:
            kinds[r["kind"]] = kinds.get(r["kind"], 0) + 1
        print(p, "scored", len(v["scored"]), kinds, "gray", len(v["gray"]), "mirror", len(v["mirror"]), "artefacts",
              len(v["artefacts"]), "germline", len(v["germline"]), "min scored af",
              min(r["af"] for r in v["scored"]), "min alt frags", min(r["alt_fragments"] for r in v["scored"]))
    print("swapped", t["swapped_patients"], "damaged", t["damaged_patients"], t["damage_mode"], "dup", t["duplication"])
    print("input MB", sum(f.stat().st_size for f in (d / "environment" / "input").rglob("*") if f.is_file()) / 1e6)


if __name__ == "__main__":
    main()
