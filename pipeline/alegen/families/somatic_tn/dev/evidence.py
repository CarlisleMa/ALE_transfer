"""Development harness (not shipped): discoverability evidence for the design notes. Reads the oracle's BAMs
(copied out of the container) with the host pysam and the hidden truth, and prints what an analyst's QC shows.

  ../.venv/bin/python -m alegen.families.somatic_tn.dev.evidence <task_dir> <dir with *.bam>
"""
import json
import sys
from collections import Counter
from pathlib import Path

import pysam


def main():
    task, bams = Path(sys.argv[1]), Path(sys.argv[2])
    t = json.loads((task / "tests" / "reference" / "truth.json").read_text())
    contig = t["contig"]
    print("== fraction of read pairs with MAPQ < 20 per target (tumour of the first patient)")
    first = sorted(t["patients"])[0]
    with pysam.AlignmentFile(str(bams / f"{t['patients'][first]['tumour_sample']}.bam")) as bam:
        for i, (s, e) in enumerate(t["targets"]):
            n = low = 0
            for r in bam.fetch(contig, s, e):
                if r.is_secondary or r.is_supplementary:
                    continue
                n += 1
                low += r.mapping_quality < 20
            print(f"   T{i + 1:02d} {s}-{e}: reads={n} mapq<20={low / n:.3f}")
    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
    for p in sorted(t["patients"]):
        d = t["patients"][p]
        with pysam.AlignmentFile(str(bams / f"{d['tumour_sample']}.bam")) as bam:
            def table(pos, alt):
                c = Counter()
                for col in bam.pileup(contig, pos - 1, pos, truncate=True, min_base_quality=20, min_mapping_quality=20,
                                      ignore_overlaps=False, max_depth=100000, stepper="nofilter"):
                    for pr in col.pileups:
                        a = pr.alignment
                        if pr.is_del or pr.is_refskip or a.is_duplicate or a.is_secondary or a.is_supplementary:
                            continue
                        if a.query_sequence[pr.query_position] != alt:
                            continue
                        c["F1R2" if a.is_read1 != a.is_reverse else "F2R1"] += 1
                        c["rev" if a.is_reverse else "fwd"] += 1
                return c
            arts = [table(a["pos"], a["alt"]) for a in d["artefacts"]]
            true = [table(v["pos"], v["alt"]) for v in d["scored"] if v["kind"] == "snv"]
        if arts:
            print(f"== {p}: artefact sites, alt reads as F1R2/F2R1 and fwd/rev")
            print("   " + "  ".join(f"{c['F1R2']}/{c['F2R1']}|{c['fwd']}/{c['rev']}" for c in arts))
            sub = Counter()
            for a in d["artefacts"]:
                sub[f"{a['ref']}>{a['alt']}" if a["ref"] in "CT" else f"{comp[a['ref']]}>{comp[a['alt']]}"] += 1
            print("   artefact substitution classes (pyrimidine reference):", dict(sub))
        mino = sorted(min(c["F1R2"], c["F2R1"]) / max(1, c["F1R2"] + c["F2R1"]) for c in true)
        mins = sorted(min(c["fwd"], c["rev"]) / max(1, c["fwd"] + c["rev"]) for c in true)
        sub = Counter()
        for v in d["scored"]:
            if v["kind"] == "snv":
                sub[f"{v['ref']}>{v['alt']}" if v["ref"] in "CT" else f"{comp[v['ref']]}>{comp[v['alt']]}"] += 1
        print(f"== {p}: {len(true)} scored SNVs: minor-orientation fraction of alt reads min={mino[0]:.2f} "
              f"median={mino[len(mino) // 2]:.2f}; minor-strand fraction min={mins[0]:.2f} median={mins[len(mins) // 2]:.2f}")
        print("   true substitution classes:", dict(sub))


if __name__ == "__main__":
    main()
