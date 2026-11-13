"""Level-4 count-level calibration (perfect quantification = planted counts).

usage: PYTHONPATH=pipeline python calibrate_l4.py 1,2,3 [out_dir]
Runs, inside the family's tool image, the level-4 oracle chain and its careless variants plus the
second solution's statistics on the planted counts, and prints the DE metrics of each together with
the composition of the "first pass" gene list (condition-only model, sample sheet as delivered).
"""
import json, shutil, sys
from pathlib import Path
import numpy as np
from alegen import dockerlib
from alegen.families import rnaseq_confound as fam
from alegen.families.rnaseq_confound import l4

FAM = Path(fam.__file__).parent
SCRIPT = r'''
set -e
P=/assets/l4/takeover_pipeline.py
run() { n=$1; shift; python3 $P --input /data/input --dev-counts /data/counts.tsv --output /out/$n --work /tmp/w "$@" > /out/$n.log 2>&1 || { echo FAILED $n; tail -5 /out/$n.log; }; }
mkdir -p /out
run oracle_fix
run oracle_exclude --labels exclude
run q_binary_identity --identity binary
run q_binary_exclude --identity binary --labels exclude
run no_batch --batch no
run labels_given --labels given
run padj_raw --padj raw
python3 /assets/l4/solve_alt_l4.py --input /data/input --counts /data/counts.tsv --output /out/alt > /out/alt.log 2>&1 || { echo FAILED alt; tail -5 /out/alt.log; }
python3 /devs/first_pass.py /data/counts.tsv /data/input/samples.tsv /out/first_pass.tsv > /out/fp.log 2>&1 || { echo FAILED fp; tail -5 /out/fp.log; }
'''

def metrics(d, truth):
    import csv
    lfc = truth["lfc"]; uns = set(truth["unscored_genes"])
    de = {r["gene_id"]: (float(r["log2fc"]), int(r["significant"])) for r in csv.DictReader(open(d / "de_results.tsv"), delimiter="\t")}
    scored = [g for g in truth["genes"] if g not in uns]; planted = [g for g in scored if lfc[g] != 0]
    tp = sum(1 for g in planted if de[g][1] and de[g][0] * lfc[g] > 0)
    fp = sum(1 for g in scored if de[g][1]) - tp
    errs = sorted(abs(de[g][0] - lfc[g]) for g in planted)
    qc = {r["sample_id"]: r["condition_used"] for r in csv.DictReader(open(d / "sample_qc.tsv"), delimiter="\t")}
    return dict(tp=tp, fp=fp, sens=round(tp / len(planted), 3), fdp=round(fp / max(tp + fp, 1), 3), err=round(errs[len(errs) // 2], 3),
                n_used=sum(v != "excluded" for v in qc.values()))

def run(seed, outroot):
    root = Path(outroot) / f"calib_l4_s{seed}"
    if root.exists(): shutil.rmtree(root)
    (root / "input" / "reference").mkdir(parents=True)
    inst = fam.build_instance(seed, 4, None)
    truth = inst["truth"]
    (root / "input" / "samples.tsv").write_text(inst["samples_tsv"])
    (root / "input" / "reference" / "tx2gene.tsv").write_text(inst["tx2gene"])
    (root / "counts.tsv").write_text(fam._tsv_counts(truth["genes"], inst["sids"], inst["file_counts"]))
    r = dockerlib.run_in_image(l4.ensure_gentools(), SCRIPT, mounts={root: "/data", FAM / "assets": "/assets", FAM / "dev": "/devs"},
                               copy_out={"/out": root / "out"}, cpus=2)
    print(f"=== seed {seed} mixed={truth['mixed']} f={truth['mixture']['fraction']} swapped={truth['swapped']}", r.stdout[-800:], r.stderr[-800:])
    for d in sorted(p for p in (root / "out").iterdir() if p.is_dir()):
        print(f"  {d.name:22s}", metrics(d, truth))
    import csv
    fp = list(csv.DictReader(open(root / "out" / "first_pass.tsv"), delimiter="\t"))
    uns = set(truth["unscored_genes"])
    rank = [r["gene_id"] for r in fp if r["gene_id"] not in uns]
    for n in (30, 40, 50, 60):
        t = sum(truth["lfc"][g] != 0 for g in rank[:n])
        print(f"  first pass top {n}: {t} planted responders, {n - t} others; batch genes among others: "
              f"{sum(truth['lfc'][g] == 0 and truth['batch_lfc'][g] != 0 for g in rank[:n])}")

if __name__ == "__main__":
    out = sys.argv[2] if len(sys.argv) > 2 else "/tmp/rnaseq_confound_calib_l4"
    for s in sys.argv[1].split(","):
        run(int(s), out)
