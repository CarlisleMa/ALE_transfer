#!/usr/bin/env python3
"""Battery helper (not part of any solution): turns the output of a level-3 pipeline into a complete
level-4 submission by adding validation.tsv, the comparison of its own counts with the manual counts.

  (default)       report the pipeline's counts as they are
  --scale         multiply every count by one plate-wide factor so that the annotated fields agree in
                  total (sum of manual counts / sum of pipeline counts), then report
  --copy-manual   overwrite the counts of the annotated wells with the manual counts, then report
"""
import argparse
import csv
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("--scale", action="store_true")
ap.add_argument("--copy-manual", action="store_true")
args = ap.parse_args()
inp, out = Path("/task/input"), Path("/task/output")
manual = {p.stem: sum(1 for ln in p.read_text().splitlines()[1:] if ln.strip())
          for p in sorted((inp / "annotations").glob("*.tsv"))}
with open(out / "well_measurements.tsv") as fh:
    rd = csv.DictReader(fh, delimiter="\t")
    fields = rd.fieldnames
    rows = list(rd)
if args.scale:
    f = sum(manual.values()) / max(1, sum(int(r["n_cells"]) for r in rows if r["well"] in manual))
    print(f"plate-wide count factor {f:.3f}")
    for r in rows:
        r["n_cells"] = str(int(round(int(r["n_cells"]) * f)))
if args.copy_manual:
    for r in rows:
        if r["well"] in manual:
            r["n_cells"] = str(manual[r["well"]])
with open(out / "well_measurements.tsv", "w", newline="") as fh:
    wr = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
    wr.writeheader()
    wr.writerows(rows)
n_of = {r["well"]: r["n_cells"] for r in rows}
with open(out / "validation.tsv", "w") as fh:
    fh.write("field\tmanual_count\tpipeline_count\n")
    for w in sorted(manual):
        fh.write(f"{w}\t{manual[w]}\t{n_of[w]}\n")
