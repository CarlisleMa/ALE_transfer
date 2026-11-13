#!/usr/bin/env python3
"""Independent second solution. Shares no code and no tool with the reference solution except samtools:

  alignment   minimap2 -ax sr (no bwa), duplicates are not removed
  pairing     bcftools mpileup | bcftools call on all six libraries, genotype discordance between
              every tumour and every normal (no allele-fraction classes)
  evidence    the tumour alignment is split by read-pair orientation with samtools flag filters; one
              `bcftools mpileup -a AD` run over (tumour F1R2, tumour F2R1, normal) gives allele depths per
              orientation and for the normal; the VCF text is parsed with the standard library (no pysam)
  filters     mapping quality (mpileup -q), base quality (-Q) with BAQ, minimum support, absence from the
              normal, both orientations among the supporting reads
"""
import csv
import itertools
import json
import os
import subprocess
from pathlib import Path

IN = Path("/task/input")
OUT = Path("/task/output")
WORK = Path(os.environ.get("WORK", "/tmp/somatic_alt"))
REF = WORK / "ref.fa"


def sh(cmd: str) -> str:
    return subprocess.run("set -o pipefail; " + cmd, shell=True, check=True, executable="/bin/bash",
                          capture_output=True, text=True).stdout


def normalise(seq, pos0, ref, alt):
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos0 = ref[1:], alt[1:], pos0 + 1
    if len(ref) != len(alt):
        while pos0 > 0 and ref[-1] == alt[-1]:
            pos0 -= 1
            ref, alt = seq[pos0] + ref[:-1], seq[pos0] + alt[:-1]
    return pos0, ref, alt


def records(vcf_text: str):
    """(pos0, ref, [alts], is_indel, {sample: [AD...]}) for every data line of an mpileup/call VCF."""
    names = []
    for line in vcf_text.splitlines():
        if line.startswith("##"):
            continue
        f = line.split("\t")
        if line.startswith("#"):
            names = [Path(x).name for x in f[9:]]
            continue
        fmt = f[8].split(":")
        per = {}
        for name, col in zip(names, f[9:]):
            per[name] = dict(zip(fmt, col.split(":")))
        yield int(f[1]) - 1, f[3].upper(), f[4].upper().split(","), "INDEL" in f[7].split(";"), per


def main() -> None:
    WORK.mkdir(parents=True, exist_ok=True)
    (OUT / "somatic").mkdir(parents=True, exist_ok=True)
    sh(f"cp {IN}/reference/ref.fa {REF} && samtools faidx {REF}")
    seq = "".join(x.strip() for x in open(REF) if not x.startswith(">")).upper()
    contig = open(REF).readline()[1:].split()[0]
    sheet = list(csv.DictReader(open(IN / "manifest.tsv"), delimiter="\t"))
    tum = {r["patient"]: r["sample_id"] for r in sheet if r["sample_type"] == "tumour"}
    normals = sorted(r["sample_id"] for r in sheet if r["sample_type"] == "normal")
    patients = sorted(tum)

    # ---- alignment, orientation split of the tumours
    for r in sheet:
        sid = r["sample_id"]
        sh(f"minimap2 -ax sr -t 2 {REF} {IN}/{r['fastq_1']} {IN}/{r['fastq_2']} 2> {WORK}/{sid}.mm2.log "
           f"| samtools sort -o {WORK}/{sid}.bam - && samtools index {WORK}/{sid}.bam")
    for sid in tum.values():
        b = f"{WORK}/{sid}"
        # F1R2: read 1 forward or read 2 reverse.  F2R1: read 1 reverse or read 2 forward.
        sh(f"samtools view -b -f 64 -F 0xF14 {b}.bam > {b}.a.bam && samtools view -b -f 144 -F 0xF04 {b}.bam > {b}.b.bam "
           f"&& samtools merge -f {b}.F1R2.bam {b}.a.bam {b}.b.bam && samtools index {b}.F1R2.bam")
        sh(f"samtools view -b -f 80 -F 0xF04 {b}.bam > {b}.c.bam && samtools view -b -f 128 -F 0xF14 {b}.bam > {b}.d.bam "
           f"&& samtools merge -f {b}.F2R1.bam {b}.c.bam {b}.d.bam && samtools index {b}.F2R1.bam")

    # ---- pairing: discordance of called genotypes at sites where a normal is not homozygous reference
    all_bams = " ".join(f"{WORK}/{s}.bam" for s in list(tum.values()) + normals)
    geno = sh(f"bcftools mpileup -f {REF} -T {IN}/targets.bed -q 20 -Q 20 -d 5000 -a AD -Ou {all_bams} 2> /dev/null "
              f"| bcftools call -mv -Ov 2> /dev/null")
    disc = {(t, n): [0, 0] for t in tum.values() for n in normals}
    for pos0, ref, alts, is_indel, per in records(geno):
        if is_indel:
            continue
        gt = {k.replace(".bam", ""): v.get("GT", "./.").replace("|", "/") for k, v in per.items()}
        gt = {k: "/".join(sorted(v.split("/"))) for k, v in gt.items()}
        if all(gt[n] in ("0/0", "./.") for n in normals):
            continue
        for (t, n), d in disc.items():
            if "." in gt[t] or "." in gt[n]:
                continue
            d[1] += 1
            d[0] += gt[t] != gt[n]
    rate = {k: (d[0] / d[1] if d[1] else 1.0) for k, d in disc.items()}
    best = min(itertools.permutations(normals), key=lambda perm: sum(rate[(tum[p], n)] for p, n in zip(patients, perm)))
    pairing = dict(zip(patients, best))
    (WORK / "pairing_qc.json").write_text(json.dumps({f"{t}|{n}": [round(rate[(t, n)], 4), disc[(t, n)][1]]
                                                      for t, n in disc}, indent=1))

    # ---- somatic calls
    targets = [(int(f[1]), int(f[2])) for f in csv.reader(open(IN / "targets.bed"), delimiter="\t") if len(f) >= 3]
    summary = {}
    for p in patients:
        t, n = tum[p], pairing[p]
        names = [f"{t}.F1R2.bam", f"{t}.F2R1.bam", f"{n}.bam"]
        pile = sh(f"bcftools mpileup -f {REF} -T {IN}/targets.bed -q 20 -Q 20 -d 100000 -L 100000 -m 3 -F 0.01 -a AD -Ov "
                  + " ".join(f"{WORK}/{x}" for x in names) + " 2> /dev/null")
        calls = {}
        for pos0, ref, alts, is_indel, per in records(pile):
            try:
                ad = [[int(x) for x in per[x]["AD"].split(",")] for x in names]
            except (KeyError, ValueError):
                continue
            t_dp, n_dp = sum(ad[0]) + sum(ad[1]), sum(ad[2])
            for i, alt in enumerate(alts, start=1):
                if alt.startswith("<") or alt in ("*", "."):
                    continue
                f1, f2, n_alt = ad[0][i], ad[1][i], ad[2][i]
                t_alt = f1 + f2
                if t_alt < 6 or t_dp < 25 or t_alt / t_dp < 0.045:
                    continue
                if n_dp < 10 or n_alt / n_dp > 0.025:
                    continue
                if min(f1, f2) < 2 or min(f1, f2) / t_alt < 0.12:
                    continue
                key = normalise(seq, pos0, ref, alt)
                if not any(s <= key[0] < e for s, e in targets):
                    continue
                if key not in calls or calls[key][0] < t_alt:
                    calls[key] = (t_alt, t_dp, f1, f2, n_alt, n_dp)
        # an SNV right next to a stronger indel is the indel drawn as a mismatch at read ends
        indels = [(k[0], v[0]) for k, v in calls.items() if len(k[1]) != len(k[2])]
        for k in [k for k in calls if len(k[1]) == len(k[2])]:
            if any(abs(k[0] - x) <= 10 and calls[k][0] < cnt for x, cnt in indels):
                del calls[k]
        vcf = WORK / f"{p}.vcf"
        with open(vcf, "w") as fh:
            fh.write("##fileformat=VCFv4.2\n" f"##contig=<ID={contig},length={len(seq)}>\n"
                     '##INFO=<ID=AF,Number=1,Type=Float,Description="Tumour allele fraction">\n'
                     '##INFO=<ID=AD_T,Number=2,Type=Integer,Description="Tumour ALT depth in F1R2 and F2R1 pairs">\n'
                     '##INFO=<ID=DP_T,Number=1,Type=Integer,Description="Tumour depth">\n'
                     '##INFO=<ID=AD_N,Number=1,Type=Integer,Description="Normal ALT depth">\n'
                     '##INFO=<ID=DP_N,Number=1,Type=Integer,Description="Normal depth">\n'
                     "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for (pos0, ref, alt), (t_alt, t_dp, f1, f2, n_alt, n_dp) in sorted(calls.items()):
                fh.write(f"{contig}\t{pos0 + 1}\t.\t{ref}\t{alt}\t.\tPASS\t"
                         f"AF={t_alt / t_dp:.4f};AD_T={f1},{f2};DP_T={t_dp};AD_N={n_alt};DP_N={n_dp}\n")
        sh(f"bgzip -c {vcf} > {OUT}/somatic/{p}.vcf.gz")
        summary[p] = {"tumour_sample": t, "normal_sample": n,
                      "n_snv": sum(1 for k in calls if len(k[1]) == len(k[2])),
                      "n_indel": sum(1 for k in calls if len(k[1]) != len(k[2]))}
    with open(OUT / "pairing.tsv", "w") as fh:
        fh.write("patient\ttumour_sample\tnormal_sample\n")
        for p in patients:
            fh.write(f"{p}\t{tum[p]}\t{pairing[p]}\n")
    (OUT / "summary.json").write_text(json.dumps({"patients": summary}, indent=1) + "\n")


if __name__ == "__main__":
    main()
