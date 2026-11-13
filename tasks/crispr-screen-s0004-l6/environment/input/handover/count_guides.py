#!/usr/bin/env python3
"""count_guides.py - sgRNA read counts for the VAL-2 screens (K.B.)

usage: count_guides.py LIBRARY.tsv SAMPLES.tsv OUT.tsv FASTQ [FASTQ ...]

One FASTQ per row of the sample sheet, given in the same order as the rows.
Only exact matches to a library spacer are counted.
"""
import gzip
import sys

COMP = str.maketrans("ACGTN", "TGCAN")
SCAFFOLD = "GTTTTAGAGCTA"


def revcomp(s):
    return s.translate(COMP)[::-1]

def candidates(read):
    # the spacer is the 20 nt in front of the scaffold; run B delivered Read 2, so look on both strands
    for r in (read, revcomp(read)):
        k = r.find(SCAFFOLD)
        if k >= 20:
            yield r[k - 20:k]


def main():
    lib_path, sheet_path, out_path = sys.argv[1:4]
    fastqs = sys.argv[4:]
    guides, by_seq = [], {}
    with open(lib_path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        i_id, i_seq, i_gene = head.index("sgRNA"), head.index("sequence"), head.index("gene")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            by_seq[f[i_seq]] = len(guides)
            guides.append((f[i_id], f[i_gene]))
    with open(sheet_path) as fh:
        fh.readline()
        samples = [line.split("\t")[0] for line in fh if line.strip()]
    if len(fastqs) != len(samples):
        sys.exit("need one FASTQ per sample sheet row")
    counts = [[0] * len(samples) for _ in guides]
    for j, path in enumerate(fastqs):
        with gzip.open(path, "rt") as fh:
            for n, line in enumerate(fh):
                if n % 4 != 1:
                    continue
                for spacer in candidates(line.rstrip()):
                    k = by_seq.get(spacer)
                    if k is not None:
                        counts[k][j] += 1
                        break
    with open(out_path, "w") as out:
        out.write("sgRNA\tgene\t" + "\t".join(samples) + "\n")
        for (sg, gene), row in zip(guides, counts):
            out.write(sg + "\t" + gene + "\t" + "\t".join(map(str, row)) + "\n")


if __name__ == "__main__":
    main()
