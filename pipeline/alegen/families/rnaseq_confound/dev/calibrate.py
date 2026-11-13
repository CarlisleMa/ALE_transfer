"""Count-level calibration: perfect quantification (true counts), every analysis variant, real R methods."""
import json, os, shutil, sys, time
from pathlib import Path
import numpy as np
from alegen import dockerlib
from alegen.families import rnaseq_confound as fam
sys.path.insert(0, str(Path(__file__).parent))
from score_dir import score, line

S = Path(__file__).resolve().parent
FAM = Path(fam.__file__).parent
TAG = os.environ.get("CALIB_TAG", "alegen-rnaseq-confound-dev:local")   # dev image = image.agent_dockerfile(image.DEV_CONDA); any task agent image works without the R methods
SCRIPT = r'''
set -e
P=/solution/rnaseq_pipeline.py
for v in "fix yes" "exclude yes" "fix no" "given yes" "given no" "cond_only yes"; do set -- $v
  python3 $P --input /data/input --counts /data/counts.tsv --labels $1 --batch $2 --output /out/py_$1_$2 --work /tmp/w > /dev/null
done
python3 /solution/solve_alt.py --input /data/input --counts /data/counts.tsv --output /out/alt_exclude > /dev/null
if [ "$WITH_R" = "1" ]; then
for m in deseq2 voom edger; do for v in fix_yes exclude_yes fix_no given_yes given_no; do
  o=/out/R_${m}_$v; mkdir -p $o; cp /out/py_$v/gene_counts.tsv /out/py_$v/sample_qc.tsv $o/
  python3 -c "
import pandas as pd
q = pd.read_csv('$o/sample_qc.tsv', sep='\t', dtype=str)
q = q[q.condition_used != 'excluded']
q[['sample_id','condition_used','batch_used']].rename(columns={'condition_used':'condition','batch_used':'batch'}).to_csv('/tmp/meta.tsv', sep='\t', index=False)"
  Rscript /devs/de_methods.R /data/counts.tsv /tmp/meta.tsv $m ${v##*_} $o/de_results.tsv 2>/dev/null
done; done
fi
'''
def run(seed, level, with_r, **kw):
    name = f"calib_s{seed}_l{level}_" + "_".join(f"{k}{v}" for k, v in sorted(kw.items()))
    root = Path(OUT) / name
    if root.exists(): shutil.rmtree(root)
    (root / "input" / "reference").mkdir(parents=True); (root / "reference").mkdir()
    inst = fam.build_instance(seed, level, None, **kw)
    (root / "input" / "samples.tsv").write_text(inst["samples_tsv"])
    (root / "input" / "reference" / "tx2gene.tsv").write_text(inst["tx2gene"])
    (root / "counts.tsv").write_text(fam._tsv_counts(inst["truth"]["genes"], inst["sids"], inst["file_counts"]))
    (root / "reference" / "truth.json").write_text(json.dumps(inst["truth"]))
    r = dockerlib.run_in_image(TAG, SCRIPT, mounts={root: "/data", FAM / "assets": "/solution", FAM / "dev": "/devs"},
                               copy_out={"/out": root / "out"}, env={"WITH_R": "1" if with_r else "0"}, cpus=2)
    if r.returncode != 0: print("CONTAINER FAILED", r.stdout[-1500:], r.stderr[-3000:])
    print(f"=== {name}  reads/sample mean {int(inst['file_counts'].sum(axis=0).mean())}")
    res = {}
    for d in sorted(p for p in (root / "out").iterdir() if p.is_dir()):
        det = score(d, root / "reference"); res[d.name] = det
        print(line(d.name, det))
    return res

OUT = "/tmp/rnaseq_confound_calib"

if __name__ == "__main__":
    # usage: PYTHONPATH=pipeline python calibrate.py 1,2 3 R [target_reads=98000 major=5 ...]
    seeds = [int(x) for x in sys.argv[1].split(",")]
    level = int(sys.argv[2]); with_r = sys.argv[3] == "R"
    kw = {}
    for a in sys.argv[4:]:
        k, v = a.split("=")
        kw[k] = float(v) if "." in v else int(v)
    for s in seeds:
        t = time.time(); run(s, level, with_r, **kw); print("  elapsed", round(time.time() - t), "s", flush=True)
