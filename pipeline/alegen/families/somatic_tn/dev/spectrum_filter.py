"""Development harness (not shipped): what does a blanket substitution-spectrum filter score?

Takes the `naive_no_orientation_filter` output of a validated task (results/work, `make --keep-work`), drops
every call of the damage substitution class below an allele-fraction cut-off (in all tumours, or only in the
tumours the notes list as FFPE / shearing batch B), and scores the result on the host.

  ../.venv/bin/python -m alegen.families.somatic_tn.dev.spectrum_filter somatic-tn-s0001-l3 [...]
"""
import gzip
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[5]
CLASSES = {"deamination": {("C", "T"), ("G", "A")}, "oxidation": {("G", "T"), ("C", "A")}}


def main():
    for task in sys.argv[1:]:
        src = ROOT / "results" / "work" / task / "attack" / "naive_no_orientation_filter" / "output"
        tests = ROOT / "tasks" / task / "tests"
        truth = json.loads((tests / "reference" / "truth.json").read_text())
        cls = set().union(*(CLASSES[m] for m in truth["damage_mode"].split("+")))
        for scope in ("all tumours", "damaged tumours only"):
            for cut in (0.15, 0.20, 0.25, 0.30, 0.35):
                out = Path(tempfile.mkdtemp())
                shutil.copytree(src, out / "output")
                summ = json.loads((out / "output" / "summary.json").read_text())
                for p in truth["patients"]:
                    if scope != "all tumours" and p not in truth["damaged_patients"]:
                        continue
                    f = out / "output" / "somatic" / f"{p}.vcf.gz"
                    keep = []
                    for line in gzip.open(f, "rt"):
                        if not line.startswith("#"):
                            c = line.split("\t")
                            af = float(c[7].split("AF=")[1].split(";")[0])
                            if (c[3], c[4]) in cls and af < cut:
                                continue
                        keep.append(line)
                    with gzip.open(f, "wt") as fh:
                        fh.writelines(keep)
                    recs = [x.split("\t") for x in keep if not x.startswith("#")]
                    summ["patients"][p]["n_snv"] = sum(1 for c in recs if len(c[3]) == len(c[4]) == 1)
                    summ["patients"][p]["n_indel"] = sum(1 for c in recs if len(c[3]) != len(c[4]))
                (out / "output" / "summary.json").write_text(json.dumps(summ))
                subprocess.run([sys.executable, str(tests / "score.py"), "--output", str(out / "output"), "--reference",
                                str(tests / "reference"), "--logs", str(out / "logs")], capture_output=True,
                               env={"PYTHONPATH": str(tests)})
                d = json.loads((out / "logs" / "details.json").read_text())
                i = d["info"]
                per = " ".join(f"{v['tp']}/{v['fp']}/{v['fn']}" for k, v in sorted(i.items()) if isinstance(v, dict))
                print(f"{task} drop {truth['damage_mode']} class below AF {cut} in {scope}: score={d['score']} "
                      f"minP={i['min_precision']} minR={i['min_recall']} TP/FP/FN {per}")
                shutil.rmtree(out)


if __name__ == "__main__":
    main()
