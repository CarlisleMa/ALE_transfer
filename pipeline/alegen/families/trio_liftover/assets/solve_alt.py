#!/usr/bin/env python3
"""Independent second solution: no Picard, no bcftools. Everything is done in Python on top of
pysam file access: REF-concordance build check, own chain-file lift-over, own allele
decomposition / left-alignment, own consequence prediction from the GFF3 + GRCh38 sequence, and a
separately written implementation of the SOP. Used to show that the verifier accepts a different
but correct route."""
import bisect
import gzip
from collections import defaultdict
from pathlib import Path

import pysam

inp, out = Path("/task/input"), Path("/task/output")
res = inp / "resources"
out.mkdir(parents=True, exist_ok=True)
fa19 = pysam.FastaFile(str(res / "GRCh37.chr22.fa.gz"))
fa38 = pysam.FastaFile(str(res / "GRCh38.chr22.fa.gz"))
CONTIG = fa38.references[0]
SEQ = {"19": fa19.fetch(CONTIG).upper(), "38": fa38.fetch(CONTIG).upper()}
COMP = str.maketrans("ACGT", "TGCA")
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {x + y + z: AA[16 * i + 4 * j + k] for i, x in enumerate("TCAG") for j, y in enumerate("TCAG")
         for k, z in enumerate("TCAG")}


def bare(c):
    return c[3:] if c.lower().startswith("chr") else c


def leftalign(seq, pos, ref, alt):
    """0-based pos; returns the left-aligned, minimal allele pair."""
    while True:
        if ref and alt and ref[-1] == alt[-1]:
            ref, alt = ref[:-1], alt[:-1]
        elif not ref or not alt:
            pos -= 1
            ref, alt = seq[pos] + ref, seq[pos] + alt
        else:
            break
    while min(len(ref), len(alt)) > 1 and ref[0] == alt[0]:
        pos, ref, alt = pos + 1, ref[1:], alt[1:]
    return pos, ref, alt


# ---------------------------------------------------------------- input, build check
vin = pysam.VariantFile(str(inp / "trio.vcf.gz"))
samples = list(vin.header.samples)
records = []
for rec in vin:
    calls = []
    for s in samples:
        c = rec.samples[s]
        calls.append((tuple(c["GT"]), tuple(c["AD"]), c["DP"], c["GQ"]))
    records.append((rec.chrom, rec.pos, rec.ref.upper(), [x.upper() for x in rec.alts], calls))
vin.close()
agree = {b: sum(SEQ[b][p - 1:p - 1 + len(r)] == r for _, p, r, _, _ in records) / len(records) for b in SEQ}
if agree["19"] <= agree["38"]:
    raise SystemExit(f"input does not look like GRCh37: {agree}")

# ---------------------------------------------------------------- chain: hg19 -> GRCh38 blocks
blocks = []          # (t0, t1, q0, plus_strand)
with open(res / "hg19ToHg38.chr22.over.chain") as fh:
    for ln in fh:
        f = ln.split()
        if not f:
            continue
        if f[0] == "chain":
            ok = bare(f[2]) == bare(CONTIG) and bare(f[7]) == bare(CONTIG)
            t, q, plus = int(f[5]), int(f[10]), f[9] == "+"
        elif ok:
            size = int(f[0])
            blocks.append((t, t + size, q, plus))
            if len(f) == 3:
                t, q = t + size + int(f[1]), q + size + int(f[2])
blocks.sort()
starts = [b[0] for b in blocks]


def lift(pos, ref):
    """0-based start of an allele whose REF spans [pos, pos+len(ref)); None if it has no unique,
    ungapped plus-strand image."""
    i = bisect.bisect_right(starts, pos) - 1
    if i < 0:
        return None
    t0, t1, q0, plus = blocks[i]
    if pos + len(ref) > t1 or not plus:
        return None
    return q0 + pos - t0


lifted, unmapped = {}, []
for chrom, pos, ref, alts, calls in records:
    for ai, alt in enumerate(alts, start=1):
        p19, r, a = leftalign(SEQ["19"], pos - 1, ref, alt)
        p38 = lift(p19, r)
        if p38 is None or SEQ["38"][p38:p38 + len(r)] != r:
            unmapped.append((chrom, pos, ref, alt))
            continue
        p38, r, a = leftalign(SEQ["38"], p38, r, a)
        per = []
        for gt, ad, dp, gq in calls:
            n_alt = sum(1 for x in gt if x == ai)
            per.append((n_alt, len(gt), ad[0], ad[ai], dp, gq))
        lifted[(p38, r, a)] = per

with open(out / "unmapped.tsv", "w") as fh:
    fh.write("chrom_in\tpos_in\tref\talt\n")
    for u in unmapped:
        fh.write("\t".join(map(str, u)) + "\n")

plain = out / "trio.GRCh38.vcf"
with open(plain, "w") as fh:
    fh.write("##fileformat=VCFv4.2\n")
    fh.write(f"##contig=<ID={CONTIG},length={len(SEQ['38'])},assembly=GRCh38>\n")
    fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    fh.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">\n')
    fh.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n')
    fh.write('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n')
    fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
    for (p, r, a) in sorted(lifted):
        cols = []
        for n_alt, ploidy, ad_ref, ad_alt, dp, gq in lifted[(p, r, a)]:
            gt = "/".join(["0"] * (ploidy - n_alt) + ["1"] * n_alt)
            cols.append(f"{gt}:{ad_ref},{ad_alt}:{dp}:{gq}")
        fh.write(f"{CONTIG}\t{p + 1}\t.\t{r}\t{a}\t.\tPASS\t.\tGT:AD:DP:GQ\t" + "\t".join(cols) + "\n")
pysam.tabix_index(str(plain), preset="vcf", force=True)     # bgzip + .tbi, removes the plain file

# ---------------------------------------------------------------- transcripts
gene_name, tx = {}, {}
with open(res / "GRCh38.chr22.canonical.gff3") as fh:
    for ln in fh:
        if ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        at = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
        if f[2] == "gene":
            gene_name[at["ID"]] = at.get("Name", at["ID"])
        elif f[2] == "mRNA":
            tx[at["ID"]] = dict(gene=gene_name[at["Parent"]], strand=f[6], cds=[], exons=[])
        elif f[2] in ("CDS", "exon") and at.get("Parent") in tx:
            tx[at["Parent"]]["cds" if f[2] == "CDS" else "exons"].append((int(f[3]) - 1, int(f[4])))
for t in tx.values():
    t["cds"].sort()
    t["exons"].sort()
    t["lo"], t["hi"] = t["exons"][0][0], t["exons"][-1][1]
    t["coding"] = [p for s, e in t["cds"] for p in range(s, e)]          # ascending genomic positions
    t["index"] = {p: i for i, p in enumerate(t["coding"])}


def consequence(t, pos, ref, alt):
    """Consequence class of one normalised allele against one transcript (None = not qualifying)."""
    seq, minus = SEQ["38"], t["strand"] == "-"
    if len(ref) == 1 and len(alt) == 1:
        if pos in t["index"]:
            i = t["index"][pos]
            n = len(t["coding"])
            j = (n - 1 - i) if minus else i                       # index in transcript orientation
            c0 = j - j % 3
            idx = [(n - 1 - k) if minus else k for k in range(c0, c0 + 3)]
            if c0 + 3 > n:
                return None
            codon = "".join(seq[t["coding"][k]] for k in idx)
            new = "".join(alt if t["coding"][k] == pos else seq[t["coding"][k]] for k in idx)
            if minus:
                codon, new = codon.translate(COMP), new.translate(COMP)
            a0, a1 = CODON.get(codon), CODON.get(new)
            if a0 is None or a1 is None or a0 == a1:
                return None
            if c0 == 0:
                return None                                       # start codon: start_lost, not on the SOP list
            if a0 == "*":
                return None                                       # stop_lost
            return "stop_gained" if a1 == "*" else "missense"
        touched = [pos]
    else:
        if len(ref) > len(alt):
            touched = list(range(pos + 1, pos + len(ref)))
            inside = all(p in t["index"] for p in touched)
        else:
            touched = [pos, pos + 1]
            inside = pos in t["index"] and pos + 1 in t["index"]
        if inside:
            return "frameshift" if (len(ref) - len(alt)) % 3 else None
    for (_, e1), (s2, _) in zip(t["exons"], t["exons"][1:]):     # intron = [e1, s2)
        five, three = {e1, e1 + 1}, {s2 - 2, s2 - 1}
        if len(ref) != len(alt) and len(ref) < len(alt):
            continue                                              # insertions next to an exon are not designed here
        if any(p in five for p in touched):
            return "splice_acceptor" if minus else "splice_donor"
        if any(p in three for p in touched):
            return "splice_donor" if minus else "splice_acceptor"
    return None


# ---------------------------------------------------------------- SOP
freq = {}
with gzip.open(res / "popfreq.GRCh38.tsv.gz", "rt") as fh:
    cols = fh.readline().split()
    for ln in fh:
        d = dict(zip(cols, ln.split()))
        freq[(int(d["pos"]), d["ref"], d["alt"])] = float(d["af"])
panel = set((inp / "panel_genes.txt").read_text().split())
ped = [ln.split() for ln in (inp / "pedigree.ped").read_text().splitlines() if ln.strip() and not ln.startswith("#")]
child = next(p for p in ped if p[5] == "2" and p[2] != "0" and p[3] != "0")
ci, fi, mi = samples.index(child[1]), samples.index(child[2]), samples.index(child[3])

reviewed = defaultdict(list)       # gene -> [(key, class, child, mother, father)]
tx_sorted = sorted(tx.values(), key=lambda t: t["lo"])
for key, per in lifted.items():
    pos, ref, alt = key
    if freq.get((pos + 1, ref, alt), 0.0) >= 0.001:
        continue
    if any(dp < 10 or gq < 20 or ploidy != 2 for _, ploidy, _, _, dp, gq in per):
        continue
    c, m, f = per[ci][0], per[mi][0], per[fi][0]
    if c == 0:
        continue
    if c == 1:
        depth = per[ci][2] + per[ci][3]
        if depth == 0 or not 0.25 <= per[ci][3] / depth <= 0.75:
            continue
    for t in tx_sorted:
        if t["gene"] in panel and t["lo"] - 10 <= pos <= t["hi"] + 10:
            cls = consequence(t, pos, ref, alt)
            if cls:
                reviewed[t["gene"]].append((key, cls, c, m, f))

rows = []
for gene, lst in reviewed.items():
    from_mother = [x for x in lst if x[2:] == (1, 1, 0)]
    from_father = [x for x in lst if x[2:] == (1, 0, 1)]
    for key, cls, c, m, f in lst:
        if (c, m, f) == (1, 0, 0):
            model = "de_novo"
        elif (c, m, f) == (2, 1, 1):
            model = "homozygous_recessive"
        elif from_mother and from_father and (c, m, f) in ((1, 1, 0), (1, 0, 1)):
            model = "compound_het"
        else:
            continue
        rows.append((key[0] + 1, key[1], key[2], gene, cls, model))
with open(out / "candidates.tsv", "w") as fh:
    fh.write("chrom\tpos\tref\talt\tgene\tconsequence\tmodel\n")
    for r in sorted(rows):
        fh.write(CONTIG + "\t" + "\t".join(map(str, r)) + "\n")
print({"ref_match": agree, "lifted": len(lifted), "unmapped": len(unmapped), "candidates": len(rows)})
