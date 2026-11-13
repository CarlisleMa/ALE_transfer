#!/usr/bin/env python3
"""Reference analysis of the trio call set:
  build check -> contig rename -> split + left-align on the source assembly -> Picard LiftoverVcf
  -> sort + left-align on GRCh38 -> bcftools csq -> frequency lookup -> SOP inheritance rules.
The non-default flags replay plausible-but-wrong analyses for the validator's attack suite; the
oracle (solve.sh) uses the defaults.
"""
import argparse
import gzip
import json
import subprocess
import sys
from pathlib import Path

import pysam

ap = argparse.ArgumentParser()
ap.add_argument("--input", default="/task/input")
ap.add_argument("--output", default="/task/output")
ap.add_argument("--work", default="/tmp/trio_work")
ap.add_argument("--no-lift", action="store_true", help="treat the provider coordinates as GRCh38")
ap.add_argument("--no-rename", action="store_true", help="run the lift-over on the contig names as delivered")
ap.add_argument("--no-norm", action="store_true", help="no splitting / left-alignment anywhere")
ap.add_argument("--ignore-quality", action="store_true", help="skip the DP/GQ and allele-balance rules")
ap.add_argument("--phase-blind", action="store_true", help="any two heterozygous hits in a gene = compound het")
a = ap.parse_args()

inp, out, work = Path(a.input), Path(a.output), Path(a.work)
res = inp / "resources"
out.mkdir(parents=True, exist_ok=True)
work.mkdir(parents=True, exist_ok=True)
QUALIFYING = ["frameshift", "stop_gained", "splice_acceptor", "splice_donor", "missense"]
MIN_DP, MIN_GQ, AB_LO, AB_HI, AF_CUT = 10, 20, 0.25, 0.75, 0.001   # lab_sop.md rules 2-4


def sh(cmd, ok_to_fail=False):
    r = subprocess.run(cmd, shell=True, executable="/bin/bash", capture_output=True, text=True)
    if r.returncode != 0 and not ok_to_fail:
        sys.exit(f"command failed: {cmd}\n{r.stderr[-3000:]}")
    return r


# ----------------------------------------------------------------------------- references
for b in ("GRCh37", "GRCh38"):
    sh(f"bgzip -dc {res}/{b}.chr22.fa.gz > {work}/{b}.fa && samtools faidx {work}/{b}.fa")
sh(f"samtools dict {work}/GRCh38.fa -o {work}/GRCh38.dict")
fa = {b: pysam.FastaFile(str(work / f"{b}.fa")) for b in ("GRCh37", "GRCh38")}
target_contig = fa["GRCh38"].references[0]

# ----------------------------------------------------------------------------- which assembly is the input on?
vin = pysam.VariantFile(str(inp / "trio.vcf.gz"))
samples = list(vin.header.samples)
in_contigs = list(vin.header.contigs)
match = {"GRCh37": 0, "GRCh38": 0}
n = 0
for rec in vin:
    if n >= 2000:
        break
    n += 1
    for b in match:
        try:
            if fa[b].fetch(target_contig, rec.pos - 1, rec.pos - 1 + len(rec.ref)).upper() == rec.ref.upper():
                match[b] += 1
        except Exception:
            pass
vin.close()
build = {k: round(v / max(n, 1), 3) for k, v in match.items()}
source = "GRCh37" if build["GRCh37"] > build["GRCh38"] else "GRCh38"

# ----------------------------------------------------------------------------- contig names
rename = {c: target_contig for c in in_contigs if c != target_contig}
cur = work / "renamed.bcf"
if rename and not a.no_rename:
    (work / "ren.txt").write_text("".join(f"{k} {v}\n" for k, v in rename.items()))
    sh(f"bcftools annotate --rename-chrs {work}/ren.txt {inp}/trio.vcf.gz -Ob -o {cur}")
else:
    sh(f"bcftools view {inp}/trio.vcf.gz -Ob -o {cur}")
(work / "samples.txt").write_text("\n".join(samples) + "\n")
final = out / "trio.GRCh38.vcf.gz"
unmapped_rows = []

if a.no_lift:
    # careless path: believe the coordinates are GRCh38; REF mismatches are only warned about
    if a.no_norm:
        sh(f"bcftools sort {cur} -Oz -o {final}")
    else:
        sh(f"bcftools norm -m-any -c w -f {work}/GRCh38.fa {cur} -Ou 2>/dev/null | bcftools sort -Oz -o {final}")
else:
    src = work / "src.vcf.gz"
    if a.no_norm or a.no_rename:
        sh(f"bcftools view {cur} -Oz -o {src}")
    else:
        sh(f"bcftools norm -m-any -c e -f {work}/{source}.fa {cur} -Oz -o {src}")
    sh(f"tabix -f -p vcf {src}")
    chain = res / ("hg19ToHg38.chr22.over.chain" if source == "GRCh37" else "hg38ToHg19.chr22.over.chain")
    sh(f"picard LiftoverVcf I={src} O={work}/lifted.vcf.gz REJECT={work}/rejected.vcf.gz CHAIN={chain} "
       f"R={work}/GRCh38.fa > {work}/picard.log 2>&1")
    # Picard writes the samples in alphabetical order: restore the provider's order, then sort and
    # bring every record into left-aligned minimal form on GRCh38
    norm = "" if a.no_norm else f"| bcftools norm -m-any -c e -f {work}/GRCh38.fa -Ou "
    sh(f"bcftools view -S {work}/samples.txt {work}/lifted.vcf.gz -Ou | bcftools sort -Ou {norm}"
       f"| bcftools view -Oz -o {final}")
    back = {v: k for k, v in rename.items()} if not a.no_rename else {}
    r = sh(f"bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' {work}/rejected.vcf.gz")
    for ln in r.stdout.splitlines():
        c, p, ref, alts = ln.split("\t")
        for alt in alts.split(","):
            unmapped_rows.append((back.get(c, c), p, ref, alt))
sh(f"tabix -f -p vcf {final}")
with open(out / "unmapped.tsv", "w") as fh:
    fh.write("chrom_in\tpos_in\tref\talt\n")
    for row in unmapped_rows:
        fh.write("\t".join(row) + "\n")

# ----------------------------------------------------------------------------- consequences
gff = res / "GRCh38.chr22.canonical.gff3"
csq = work / "csq.vcf"
r = sh(f"bcftools csq -l -f {work}/GRCh38.fa -g {gff} {final} -Ov -o {csq} 2> {work}/csq.log", ok_to_fail=True)
if r.returncode != 0:
    # Only reachable on the careless paths: csq stops with "the fasta reference does not match the
    # VCF REF allele". The careless reaction is to drop the offending records and carry on.
    kept = work / "ref_consistent.vcf"
    with gzip.open(final, "rt") as fi, open(kept, "w") as fo:
        for ln in fi:
            if not ln.startswith("#"):
                f = ln.split("\t", 5)
                try:
                    if fa["GRCh38"].fetch(f[0], int(f[1]) - 1, int(f[1]) - 1 + len(f[3])).upper() != f[3].upper():
                        continue
                except Exception:
                    continue
            fo.write(ln)
    sh(f"bcftools csq -l -f {work}/GRCh38.fa -g {gff} {kept} -Ov -o {csq} 2> {work}/csq.log")

af = {}
with gzip.open(res / "popfreq.GRCh38.tsv.gz", "rt") as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    ix = {h: i for i, h in enumerate(hdr)}
    for ln in fh:
        f = ln.rstrip("\n").split("\t")
        af[(f[ix["chrom"]], int(f[ix["pos"]]), f[ix["ref"]], f[ix["alt"]])] = float(f[ix["af"]])

panel = {g.strip() for g in (inp / "panel_genes.txt").read_text().split() if g.strip()}
proband = mother = father = None
for ln in (inp / "pedigree.ped").read_text().splitlines():
    f = ln.split()
    if len(f) >= 6 and not ln.startswith("#") and f[5] == "2" and f[2] != "0" and f[3] != "0":
        proband, father, mother = f[1], f[2], f[3]
if proband is None:
    sys.exit("no affected child with both parents in pedigree.ped")

hits = []          # dicts: key, gene, cls, genotype counts
with open(csq) as fh:
    for ln in fh:
        if ln.startswith("##"):
            continue
        f = ln.rstrip("\n").split("\t")
        if ln.startswith("#"):
            col = {s: i for i, s in enumerate(f)}
            continue
        chrom, pos, ref, alts = f[0], int(f[1]), f[3], f[4].split(",")
        info = dict(x.split("=", 1) for x in f[7].split(";") if "=" in x)
        if "BCSQ" not in info:
            continue
        keys = f[8].split(":")
        call = {}
        for s in (proband, mother, father):
            call[s] = dict(zip(keys, f[col[s]].split(":")))
        for i, alt in enumerate(alts, start=1):
            cons = []
            for e in info["BCSQ"].split(","):
                p = e.split("|")
                if len(p) < 2:
                    continue
                if len(alts) > 1 and not (len(p) >= 7 and p[6].endswith(">" + alt)):
                    continue
                toks = [t.lstrip("*") for t in p[0].split("&")]
                cls = next((q for q in QUALIFYING if q in toks), None)
                if cls and p[1] in panel and (p[1], cls) not in cons:
                    cons.append((p[1], cls))
            if not cons:
                continue
            if af.get((chrom, pos, ref, alt), 0.0) >= AF_CUT:
                continue
            g, usable = {}, True
            for s in (proband, mother, father):
                gt = call[s].get("GT", "./.").replace("|", "/").split("/")
                if "." in gt:
                    usable = False
                    break
                g[s] = sum(1 for x in gt if x == str(i))
                other = sum(1 for x in gt if x not in ("0", str(i)))
                try:
                    dp, gq = int(call[s].get("DP", "0")), int(call[s].get("GQ", "0"))
                except ValueError:
                    dp, gq = 0, 0
                if other:
                    usable = False
                if not a.ignore_quality and (dp < MIN_DP or gq < MIN_GQ):
                    usable = False
            if not usable:
                continue
            if g[proband] == 1 and not a.ignore_quality:
                try:
                    ad = [int(x) for x in call[proband].get("AD", "").split(",")]
                    ab = ad[i] / (ad[0] + ad[i])
                except Exception:
                    continue
                if not AB_LO <= ab <= AB_HI:
                    continue
            for gene, cls in cons:
                hits.append(dict(key=(chrom, pos, ref, alt), gene=gene, cls=cls, g=(g[proband], g[mother], g[father])))

rows, by_gene = [], {}
for h in hits:
    if h["g"] == (1, 0, 0):
        rows.append((h, "de_novo"))
    elif h["g"] == (2, 1, 1):
        rows.append((h, "homozygous_recessive"))
    elif h["g"][0] == 1:
        by_gene.setdefault(h["gene"], []).append(h)
for gene, lst in by_gene.items():
    mat = [h for h in lst if h["g"][1:] == (1, 0)]
    pat = [h for h in lst if h["g"][1:] == (0, 1)]
    if mat and pat:
        rows += [(h, "compound_het") for h in mat + pat]
    elif a.phase_blind and len(lst) >= 2:
        rows += [(h, "compound_het") for h in lst]
rows.sort(key=lambda x: x[0]["key"])
with open(out / "candidates.tsv", "w") as fh:
    fh.write("chrom\tpos\tref\talt\tgene\tconsequence\tmodel\n")
    for h, model in rows:
        fh.write("\t".join(map(str, h["key"])) + f"\t{h['gene']}\t{h['cls']}\t{model}\n")
print(json.dumps({"ref_match": build, "source": source, "unmapped_alleles": len(unmapped_rows),
                  "candidates": len(rows)}))
