"""Development helper for level 4: run the oracle, the second solution and careless variants on a
generated task directory inside the family's tool image and score them on the host.

usage: PYTHONPATH=pipeline python run_l4_variants.py <task_dir> <out_dir> [name ...]
"""
import json, subprocess, sys, tempfile, time
from pathlib import Path
from alegen import dockerlib
from alegen.families import rnaseq_confound as fam
from alegen.families.rnaseq_confound import l4

FAM = Path(fam.__file__).parent
CORE = FAM.parents[1] / "scoring.py"

def score(out, ref):
    logs = tempfile.mkdtemp(); lib = Path(tempfile.mkdtemp())
    (lib / "_scoring.py").write_text(CORE.read_text())
    r = subprocess.run([sys.executable, str(FAM / "assets" / "l4" / "score.py"), "--output", str(out), "--reference", str(ref), "--logs", logs],
                       capture_output=True, text=True, env={"PYTHONPATH": str(lib), "PATH": "/usr/bin:/bin"})
    f = Path(logs) / "details.json"
    return json.loads(f.read_text()) if f.exists() else {"error": (r.stdout + r.stderr)[-800:]}

def line(name, d):
    if "error" in d: return f"{name:34s} VERIFIER ERROR {d['error']}"
    if d.get("failure"): return f"{name:34s} score=0 FAIL {d['failure'][:140]}"
    i = d["info"]; c = {x["name"]: x["value"] for x in d["components"]}
    return (f"{name:34s} score={d['score']:.3f} cnt={i['count_agreement_min']:.3f} sex={i['sex_correct']} lab={i['labels_valid']} "
            f"TP={i['de_tp']}/{i['de_planted']} FP={i['de_fp']} sens={i['de_sensitivity']:.3f} fdp={i['de_fdp']:.3f} err={i['lfc_median_abs_error']:.3f} "
            f"audit_mcc={i['audit_mcc']:.2f}(miss {i['audit_missed']}, alarm {i['audit_false_alarms']}) erratum_mcc={i['erratum_mcc']:.2f} "
            f"| " + " ".join(f"{k[:9]}={v:.2f}" for k, v in c.items()))

def main():
    task, out = Path(sys.argv[1]), Path(sys.argv[2])
    truth = json.loads((task / "tests" / "reference" / "truth.json").read_text())
    plan = truth["plan"]
    P = "python3 /solution/takeover_pipeline.py"
    variants = {
        "oracle": P,
        "alt": "python3 /solution/solve_alt_l4.py",
        "oracle_exclude": P + " --labels exclude",
        "l3_reference": "python3 /l3/rnaseq_pipeline.py",
        "trust_inherited": P + " --counts-from inherited --audit all_ok",
        "q_identity_not_applied": P + " --identity binary",
        "naive_labels_as_given": P + " --labels given",
        "kept_quant_library_type": P + f" --libtype {plan['libtype_bad']} --audit-force quant_library_type=ok",
        "kept_merge_read_counts": P + " --quantity TPM --audit-force merge_read_counts=ok",
        "kept_matrix_sample_files": P + f" --dup {plan['dup'][0]}:{plan['dup'][1]} --audit-force matrix_sample_files=ok",
        "kept_matrix_all_reads": P + f" --trunc {plan['trunc'][0]}:{plan['trunc'][1]} --audit-force matrix_all_reads=ok",
        "kept_de_session_term": P + " --batch no --audit-force de_session_term=ok",
        "kept_de_contrast_direction": P + " --contrast reversed --audit-force de_contrast_direction=ok",
        "kept_padj_raw": P + " --padj raw",
    }
    names = sys.argv[3:] or list(variants)
    script = "set -e\nmkdir -p /out\n"
    for n in names:
        script += (f"rm -rf /task/output/* /tmp/takeover_work /tmp/rnaseq_work; S=$(date +%s); ({variants[n]}) > /out/{n}.log 2>&1 || echo FAILED {n}; "
                   f"echo $(( $(date +%s) - S )) > /out/{n}.seconds; mkdir -p /out/{n}; cp -r /task/output/. /out/{n}/\n")
    t = time.time()
    r = dockerlib.run_in_image(l4.ensure_gentools(), script, mounts={task / "environment" / "input": "/task/input", task / "solution": "/solution",
                               FAM / "assets": "/l3"}, copy_out={"/out": out}, cpus=2)
    print(r.stdout[-1500:], r.stderr[-1500:], "container", round(time.time() - t), "s")
    print("planted defects:", plan["defects"], "| mixed:", truth["mixed"], "| swapped:", truth["swapped"])
    for n in names:
        d = out / n
        if n == "l3_reference":      # the careless completion the validator's attack adds
            import csv
            sig = {x["gene_id"]: x["significant"] for x in csv.DictReader(open(d / "de_results.tsv"), delimiter="\t")}
            (d / "audit.tsv").write_text("item\tverdict\n" + "".join(f"{i}\tok\n" for i in truth["audit"]))
            (d / "erratum.tsv").write_text("gene_id\tverdict\n" + "".join(f"{g}\t{'confirmed' if sig[g] == '1' else 'retracted'}\n" for g in truth["mailed"]))
        secs = (out / f"{n}.seconds").read_text().strip() if (out / f"{n}.seconds").exists() else "?"
        print(line(n, score(d, task / "tests" / "reference")), f"[{secs}s]")
        if n in ("oracle", "alt"):
            print("    audit:", {x.split("\t")[0]: x.split("\t")[1] for x in (d / "audit.tsv").read_text().splitlines()[1:]} if (d / "audit.tsv").exists() else "missing")

if __name__ == "__main__":
    main()
