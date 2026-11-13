#!/usr/bin/env python3
"""Reference analysis of the trio call set, level 4.

External tool invocations (24), in order:
   1-4   bgzip -dc and samtools faidx for both assemblies
         [Python/pysam: REF concordance of the input against both assemblies]
   5-6   bcftools annotate --rename-chrs, tabix
   7     bcftools +liftover on the records as delivered. Multi-allelic records are NOT split before the
         lift-over: splitting first turns "the other ALT" into the provider's reference allele, which is
         wrong as soon as the reference allele changes.
   8-11  bcftools sort, bcftools norm -m-any on GRCh38, bcftools view -i 'GT="alt"' (records that nobody
         carries any more are dropped), tabix
   12    bcftools query on the reject file -> unmapped.tsv
   13    bcftools query: provider positions covered by a record
         [Python: walk the chain blocks, compare the two assemblies base by base, add the positions
          without a record as homozygous ALT alleles]
   14-19 bcftools view -h, bgzip, tabix, bcftools concat -a, bcftools sort, tabix of the deliverable
   20    bcftools csq -l (every variant on its own)
         [Python: chromosome copy of each variant of the proband from the trio genotypes]
   21-24 bcftools view -h, bcftools csq on one VCF per chromosome copy (maternal, paternal, unassigned)
followed by the SOP in Python (frequency lookup, genotype rules, inheritance models).

The non-default flags replay plausible-but-wrong analyses for the validator's attack suite; the
oracle (solve.sh) uses the defaults.
"""
import argparse
import bisect
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
ap.add_argument("--liftover", default="plugin", choices=["plugin", "picard"],
                help="picard: LiftoverVcf with defaults (records whose REF is not the target base are rejected)")
ap.add_argument("--no-hidden", action="store_true", help="do not add positions without a provider record")
ap.add_argument("--split-first", action="store_true", help="split multi-allelic records before the lift-over")
ap.add_argument("--drop-inverted", choices=["unmapped", "silent"], default=None,
                help="use plus-strand chains only; the records of minus-strand blocks are reported as unmapped / lost")
ap.add_argument("--per-variant", action="store_true", help="consequence of every variant on its own")
ap.add_argument("--merge-copies", action="store_true", help="codon partners merged whatever copy carries them")
ap.add_argument("--provider-genotypes", action="store_true",
                help="inheritance patterns read from the genotypes as written against the provider's reference")
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


# ----------------------------------------------------------------------------- references, build check
for b in ("GRCh37", "GRCh38"):
    sh(f"bgzip -dc {res}/{b}.chr22.fa.gz > {work}/{b}.fa && samtools faidx {work}/{b}.fa")
fa = {b: pysam.FastaFile(str(work / f"{b}.fa")) for b in ("GRCh37", "GRCh38")}
target_contig = fa["GRCh38"].references[0]
vin = pysam.VariantFile(str(inp / "trio.vcf.gz"))
samples = list(vin.header.samples)
in_contigs = list(vin.header.contigs)
match, n = {"GRCh37": 0, "GRCh38": 0}, 0
for rec in vin:
    n += 1
    for b in match:
        if fa[b].fetch(target_contig, rec.pos - 1, rec.pos - 1 + len(rec.ref)).upper() == rec.ref.upper():
            match[b] += 1
vin.close()
build = {k: round(v / max(n, 1), 3) for k, v in match.items()}
source = "GRCh37" if build["GRCh37"] > build["GRCh38"] else "GRCh38"
if source != "GRCh37":
    sys.exit(f"the input does not look like a GRCh37 call set: {build}")
chain = res / "hg19ToHg38.chr22.over.chain"
COMP = str.maketrans("ACGTN", "TGCAN")
minus_targets = []
if a.drop_inverted:
    kept, keep = [], True
    for ln in open(chain):
        f = ln.split()
        if f and f[0] == "chain":
            keep = f[9] == "+"
            if not keep:
                minus_targets.append((int(f[5]), int(f[6])))
        if keep:
            kept.append(ln)
    chain = work / "plus_strand_only.chain"
    chain.write_text("".join(kept))

rename = {c: target_contig for c in in_contigs if c != target_contig}
back = {v: k for k, v in rename.items()}
(work / "ren.txt").write_text("".join(f"{k} {v}\n" for k, v in rename.items()) or f"{target_contig} {target_contig}\n")
src = work / "renamed.vcf.gz"
sh(f"bcftools annotate --rename-chrs {work}/ren.txt {inp}/trio.vcf.gz -Oz -o {src} && tabix -f -p vcf {src}")
(work / "samples.txt").write_text("\n".join(samples) + "\n")

# ----------------------------------------------------------------------------- lift-over of the records
lifted = work / "lifted.vcf.gz"
norm = "" if a.no_norm else f"| bcftools norm -m-any -c e -f {work}/GRCh38.fa -Ou 2>{work}/norm.log "
if a.split_first:
    sh(f"bcftools norm -m-any -c e -f {work}/GRCh37.fa {src} -Oz -o {work}/split_first.vcf.gz "
       f"&& tabix -f -p vcf {work}/split_first.vcf.gz")
    src_lift = work / "split_first.vcf.gz"
else:
    src_lift = src
if a.liftover == "plugin":
    sh(f"bcftools +liftover --no-version {src_lift} -Ou -- -s {work}/GRCh37.fa -f {work}/GRCh38.fa -c {chain} "
       f"--reject {work}/rejected.vcf.gz --reject-type z 2>{work}/liftover.log "
       f"| bcftools sort -Ou 2>/dev/null {norm}| bcftools view -i 'GT=\"alt\"' -Oz -o {lifted}")
else:
    sh(f"samtools dict {work}/GRCh38.fa -o {work}/GRCh38.dict")
    sh(f"bcftools norm -m-any -c e -f {work}/GRCh37.fa {src} -Oz -o {work}/split.vcf.gz && tabix -f -p vcf {work}/split.vcf.gz")
    sh(f"picard LiftoverVcf I={work}/split.vcf.gz O={work}/picard.vcf.gz REJECT={work}/rejected.vcf.gz CHAIN={chain} "
       f"R={work}/GRCh38.fa > {work}/picard.log 2>&1")
    sh(f"bcftools view -S {work}/samples.txt {work}/picard.vcf.gz -Ou | bcftools sort -Ou 2>/dev/null {norm}"
       f"| bcftools view -i 'GT=\"alt\"' -Oz -o {lifted}")
sh(f"tabix -f -p vcf {lifted}")
unmapped_rows = []
r = sh(f"bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' {work}/rejected.vcf.gz")
for ln in r.stdout.splitlines():
    c, p, ref, alts = ln.split("\t")
    if a.drop_inverted == "silent" and any(lo < int(p) <= hi for lo, hi in minus_targets):
        continue
    for alt in alts.split(","):
        unmapped_rows.append((back.get(c, c), p, ref, alt))
with open(out / "unmapped.tsv", "w") as fh:
    fh.write("chrom_in\tpos_in\tref\talt\n")
    for row in unmapped_rows:
        fh.write("\t".join(row) + "\n")

# ----------------------------------------------------------------------------- positions without a record
# The provider's file is complete against the provider's reference: where it has no record, the family
# is homozygous for the provider's base. Inside aligned blocks such a base is a homozygous ALT allele
# on GRCh38 wherever the two assemblies differ.
final = out / "trio.GRCh38.vcf.gz"
n_added = 0
if a.no_hidden:
    sh(f"cp {lifted} {final}")
else:
    covered = set()
    r = sh(f"bcftools query -f '%POS\\t%REF\\n' {src}")
    for ln in r.stdout.splitlines():
        p, ref = ln.split("\t")
        covered.update(range(int(p), int(p) + len(ref)))
    s19 = fa["GRCh37"].fetch(target_contig).upper()
    s38 = fa["GRCh38"].fetch(target_contig).upper()
    rows = []
    with open(chain) as fh:
        for ln in fh:
            f = ln.split()
            if not f:
                continue
            if f[0] == "chain":
                plus, qsize = f[9] == "+", int(f[8])
                t, q = int(f[5]), int(f[10])
                continue
            size = int(f[0])
            x = s19[t:t + size]
            # minus-strand chains count query coordinates from the end of the sequence: compare with the
            # reverse complement of GRCh38 and report on the plus strand
            y = s38[q:q + size] if plus else s38[qsize - q - size:qsize - q].translate(COMP)[::-1]
            if x != y:
                for i in range(size):
                    if x[i] != y[i] and x[i] in "ACGT" and y[i] in "ACGT" and (t + i + 1) not in covered:
                        if plus:
                            rows.append((q + i + 1, y[i], x[i]))
                        else:
                            rows.append((qsize - q - i, y[i].translate(COMP), x[i].translate(COMP)))
            if len(f) == 3:
                t, q = t + size + int(f[1]), q + size + int(f[2])
    n_added = len(rows)
    hdr = sh(f"bcftools view -h {lifted}").stdout
    with open(work / "added.vcf", "w") as fh:
        fh.write(hdr)
        for p, ref, alt in sorted(set(rows)):
            fh.write(f"{target_contig}\t{p}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t" + "\t".join(["1/1"] * len(samples)) + "\n")
    sh(f"bgzip -f {work}/added.vcf && tabix -f -p vcf {work}/added.vcf.gz")
    sh(f"bcftools concat -a {lifted} {work}/added.vcf.gz -Ou | bcftools sort -Oz -o {final} 2>/dev/null")
sh(f"tabix -f -p vcf {final}")

# ----------------------------------------------------------------------------- pedigree, call set in memory
panel = {g.strip() for g in (inp / "panel_genes.txt").read_text().split() if g.strip()}
proband = mother = father = None
for ln in (inp / "pedigree.ped").read_text().splitlines():
    f = ln.split()
    if len(f) >= 6 and not ln.startswith("#") and f[5] == "2" and f[2] != "0" and f[3] != "0":
        proband, father, mother = f[1], f[2], f[3]
if proband is None:
    sys.exit("no affected child with both parents in pedigree.ped")
trio = (proband, mother, father)

recs = []           # dicts: key, copies (P, M, F), usable, balance_ok, swap
with gzip.open(final, "rt") as fh:
    for ln in fh:
        if ln.startswith("##"):
            continue
        f = ln.rstrip("\n").split("\t")
        if ln.startswith("#"):
            col = {s: i for i, s in enumerate(f)}
            continue
        alts = f[4].split(",")
        keys = f[8].split(":")
        call = {s: dict(zip(keys, f[col[s]].split(":"))) for s in trio}
        added = keys == ["GT"]                 # position without a provider record: usable by SOP rule 2
        for i, alt in enumerate(alts, start=1):
            copies, usable = [], True
            for s in trio:
                gt = call[s].get("GT", "./.").replace("|", "/").split("/")
                if "." in gt:
                    usable = False
                copies.append(sum(1 for x in gt if x == str(i)))
                if not added and not a.ignore_quality:
                    try:
                        dp, gq = int(call[s].get("DP", "0")), int(call[s].get("GQ", "0"))
                    except ValueError:
                        dp, gq = 0, 0
                    if dp < MIN_DP or gq < MIN_GQ:
                        usable = False
            balance_ok = True
            if copies[0] == 1 and not a.ignore_quality and not added:
                ad = call[proband].get("AD", "").split(",")
                try:
                    alt_n = int(ad[i])
                    ref_n = 0 if ad[0] == "." else int(ad[0])      # reference allele added by the lift-over
                    balance_ok = AB_LO <= alt_n / (ref_n + alt_n) <= AB_HI
                except Exception:
                    balance_ok = False
            recs.append(dict(key=(f[0], int(f[1]), f[3], alt), copies=tuple(copies), usable=usable,
                             balance_ok=balance_ok, swap="SWAP=1" in f[7].split(";"), multi=len(alts) > 1))

# ----------------------------------------------------------------------------- consequences
gff = res / "GRCh38.chr22.canonical.gff3"


def parse_bcsq(path):
    """{(pos, ref, alt): [(class, gene, set of 'posREF>ALT' the entry is made of)]} with '@pos' pointers resolved."""
    raw, order = {}, []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            info = dict(x.split("=", 1) for x in f[7].split(";") if "=" in x)
            for alt in f[4].split(","):
                k = (int(f[1]), f[3], alt)
                raw[k] = (info.get("BCSQ", ""), len(f[4].split(",")) > 1)
                order.append(k)
    out_ = {}
    by_pos = {}
    for k in order:
        by_pos.setdefault(k[0], []).append(k)
    for k in order:
        text, multi = raw[k]
        me = f"{k[0]}{k[1]}>{k[2]}"
        cons = []
        for e in text.split(","):
            if not e or e == ".":
                continue
            if e.startswith("@"):
                for k2 in by_pos.get(int(e[1:]), []):
                    for e2 in raw[k2][0].split(","):
                        p2 = e2.split("|")
                        if len(p2) >= 7 and me in p2[6].split("+"):
                            cons.append((p2[0], p2[1], set(p2[6].split("+"))))
                continue
            p = e.split("|")
            if len(p) < 2:
                continue
            if multi and not (len(p) >= 7 and p[6].endswith(">" + k[2])):
                continue
            cons.append((p[0], p[1], set(p[6].split("+")) if len(p) >= 7 else {me}))
        out_[k] = cons
    return out_


def classes(cons):
    """[(gene, qualifying class or None)] per gene of a consequence list."""
    res_ = {}
    for text, gene, _ in cons:
        toks = [t.lstrip("*") for t in text.split("&")]
        cls = next((q for q in QUALIFYING if q in toks), None)
        if gene not in res_ or (res_[gene] is None and cls is not None):
            res_[gene] = cls
    return res_


csq_l = work / "csq_l.vcf"
sh(f"bcftools csq -l -f {work}/GRCh38.fa -g {gff} {final} -Ov -o {csq_l} 2> {work}/csq_l.log")
alone = parse_bcsq(csq_l)

# which chromosome copy of the proband carries a variant (SOP rule 5)
for r_ in recs:
    c, m, f_ = r_["copies"]
    r_["copy"] = None if c == 0 else ("both" if c == 2 else ("mat" if m > 0 and f_ == 0 else ("pat" if f_ > 0 and m == 0 else "unk")))
if a.merge_copies:
    for r_ in recs:
        if r_["copy"] is not None:
            r_["copy"] = "both"
snv = [r_ for r_ in recs if r_["copy"] is not None and len(r_["key"][2]) == 1 and len(r_["key"][3]) == 1 and not r_["multi"]]
runs = {"mat": [], "pat": []}
unk = [r_ for r_ in snv if r_["copy"] == "unk"]
extra = []                                   # unassigned variants never share a run with a close neighbour
for r_ in sorted(unk, key=lambda x: x["key"][1]):
    for lst in extra:
        if all(abs(o["key"][1] - r_["key"][1]) > 3 for o in lst):
            lst.append(r_)
            break
    else:
        extra.append([r_])
for i, lst in enumerate(extra or [[]]):
    runs[f"unk{i}"] = lst
hap = {}                                     # key -> {run name: consequence list}
if not a.per_variant:
    hdr = [ln for ln in sh(f"bcftools view -h {final}").stdout.splitlines() if not ln.startswith("#CHROM")]
    for name, own in runs.items():
        members = [r_ for r_ in snv if r_["copy"] == "both"] + ([r_ for r_ in snv if r_["copy"] == name] if name in ("mat", "pat") else own)
        seen, body = set(), []
        for r_ in sorted(members, key=lambda x: x["key"]):
            if r_["key"][1] in seen:         # two alleles of one position cannot share a copy
                continue
            seen.add(r_["key"][1])
            k = r_["key"]
            body.append(f"{k[0]}\t{k[1]}\t.\t{k[2]}\t{k[3]}\t.\tPASS\t.\tGT\t1|1")
        vcf = work / f"copy_{name}.vcf"
        vcf.write_text("\n".join(hdr + ["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCOPY"] + body) + "\n")
        sh(f"bcftools csq -f {work}/GRCh38.fa -g {gff} {vcf} -Ov -o {work}/copy_{name}.csq.vcf 2> {work}/copy_{name}.log")
        got = parse_bcsq(work / f"copy_{name}.csq.vcf")
        for r_ in members:
            k = r_["key"]
            if (k[1], k[2], k[3]) in got:
                hap.setdefault(k, {})[name] = got[(k[1], k[2], k[3])]


def qualifying(r_):
    """{gene: class} of a record under rule 5 (None if it does not qualify)."""
    k = r_["key"]
    own = classes(alone.get((k[1], k[2], k[3]), []))
    per_run = hap.get(k)
    if not per_run:
        return {g: c for g, c in own.items() if c}
    # hap[k] holds one entry per run the record is a member of: its own copy for a variant carried once,
    # every copy for a variant carried twice (which must qualify on all of them)
    res_ = {}
    for g in set(own) | {g for cons in per_run.values() for g in classes(cons)}:
        got = []
        for cons in per_run.values():
            c = classes(cons)
            got.append(c[g] if g in c else own.get(g))
        if got and all(c is not None for c in got):
            res_[g] = "&".join(sorted(set(got)))
    return res_


af = {}
with gzip.open(res / "popfreq.GRCh38.tsv.gz", "rt") as fh:
    hdr_ = fh.readline().rstrip("\n").split("\t")
    ix = {h: i for i, h in enumerate(hdr_)}
    for ln in fh:
        f = ln.rstrip("\n").split("\t")
        af[(f[ix["chrom"]], int(f[ix["pos"]]), f[ix["ref"]], f[ix["alt"]])] = float(f[ix["af"]])

hits = []
for r_ in recs:
    if r_["copies"][0] == 0 or not r_["usable"] or not r_["balance_ok"]:
        continue
    if af.get(r_["key"], 0.0) >= AF_CUT:
        continue
    g = r_["copies"]
    if a.provider_genotypes and r_["swap"]:
        g = tuple(2 - x for x in g)          # the pattern as the provider's file shows it
    for gene, cls in qualifying(r_).items():
        if gene in panel:
            hits.append(dict(key=r_["key"], gene=gene, cls=cls, g=g))

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
                  "added_sites": n_added, "copy_runs": sorted(runs), "candidates": len(rows)}))
