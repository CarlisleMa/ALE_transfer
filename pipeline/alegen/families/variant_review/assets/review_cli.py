#!/usr/bin/env python3
"""Independent second solution of level 4: IGV's command-line route. No mouse, no alignment export.

  counts       `igvtools count -w 1 --bases --strands read` on each CRAM (duplicates excluded by
               default; a second run with `--minMapQuality 20` for the mapping-quality rule):
               depth, alternate molecules per strand, deletions, the normal
  read level   what allele counts cannot give (base quality of the alternate base, position in the
               read, further mismatches, insert size of the pair, whether two adjacent alleles sit on
               the same reads) is read from SVG images of a 221-bp view around each candidate, written
               through IGV's batch port: every read is a polygon, every mismatch a coloured
               rectangle whose opacity is IGV's base-quality shading (opaque from quality 20), a
               deletion is a line between two blocks, pairs with an insert above the display
               threshold (1,000 bp, insert-size computation switched off) are red
  documents    batch port: colorBy READ_STRAND, goto, saveSession, snapshot

`--counts-only` replays a careless review for the validator (plain `igvtools count --bases`, only
GERMLINE and LOW_SUPPORT judged, no session, no images).
"""
import argparse
import re
import socket
import subprocess
import time
from pathlib import Path

INP, OUT, WORK = Path("/task/input"), Path("/task/output"), Path("/tmp/igv_cli")
HALF = 110
BASE_COLOUR = {(0, 255, 0): "A", (0, 0, 255): "C", (209, 113, 5): "G", (255, 0, 0): "T"}
NAMED = {"red": (255, 0, 0), "blue": (0, 0, 255), "lime": (0, 255, 0), "black": (0, 0, 0), "white": (255, 255, 255),
         "gray": (128, 128, 128), "silver": (192, 192, 192), "green": (0, 128, 0), "yellow": (255, 255, 0),
         "aqua": (0, 255, 255), "fuchsia": (255, 0, 255), "maroon": (128, 0, 0), "navy": (0, 0, 128),
         "olive": (128, 128, 0), "purple": (128, 0, 128), "teal": (0, 128, 128)}
LONG_INSERT = (200, 0, 0)
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:6.1f}s] {msg}", flush=True)


# ------------------------------------------------------------------ igvtools count

def count(cram, tag, *opts):
    """{1-based position: ({base: n} forward, {base: n} reverse)} from `igvtools count --bases`."""
    wig = WORK / f"{tag}.wig"
    cmd = ["igvtools", "count", "-w", "1", "--bases", *opts, str(cram), str(wig), str(INP / "reference" / "ref.fa")]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if not wig.exists():
        raise RuntimeError("igvtools count failed: " + (r.stdout + r.stderr)[-500:])
    keys = ("A", "C", "G", "T", "N", "DEL", "INS")
    out = {}
    for line in open(wig):
        f = line.split("\t")
        if len(f) < 8 or not f[0].isdigit():
            continue
        v = [int(float(x)) for x in f[1:]]
        if len(v) >= 14:
            out[int(f[0])] = (dict(zip(keys, v[:7])), dict(zip(keys, v[7:14])))
        else:                       # without --strands: one set of columns
            out[int(f[0])] = (dict(zip(keys, v[:7])), dict.fromkeys(keys, 0))
    return out


def depth(cnt, pos1):
    f, r = cnt.get(pos1, ({}, {}))
    return sum(f.get(k, 0) + r.get(k, 0) for k in ("A", "C", "G", "T", "N", "DEL"))


# ------------------------------------------------------------------ batch port

class Batch:
    def __init__(self):
        self.sock = socket.create_connection(("127.0.0.1", 60151), timeout=600)
        self.io = self.sock.makefile("rw")

    def __call__(self, cmd):
        self.io.write(cmd + "\n")
        self.io.flush()
        ans = self.io.readline().strip()
        if not ans.startswith("OK") and cmd != "echo":
            log(f"batch: {cmd} -> {ans}")
        return ans


def start_igv():
    """Start IGV and return a batch connection once the three start-up tracks are really there (the
    port answers before IGV has loaded them: a session saved to a scratch file tells)."""
    subprocess.run(["bash", "/task/software/open_igv.sh"], capture_output=True)
    probe = WORK / "ready.xml"
    end = time.time() + 300
    while time.time() < end:
        try:
            b = Batch()
            if b("echo") == "echo":
                b(f"saveSession {probe}")
                txt = probe.read_text() if probe.exists() else ""
                if txt.count("<Resource ") >= 3 and txt.count("AlignmentTrack") >= 2:
                    return b
            b.sock.close()
        except OSError:
            pass
        time.sleep(1.0)
    raise RuntimeError("IGV did not come up with its tracks")


# ------------------------------------------------------------------ SVG

def colour(s):
    if not s or s == "none":
        return None
    m = re.match(r"rgb\((\d+),(\d+),(\d+)\)", s)
    return tuple(map(int, m.groups())) if m else NAMED.get(s)


def parse_svg(path):
    """Shapes of the first (tumour) alignment panel and the view's pixels per base.
    Returns (blocks, gaps, marks, scale):
      blocks  (row, x_left, x_right, colour, strand or None)
      gaps    (row, x1, x2)                 deletion drawn between two blocks
      marks   (row, x, base, shaded)        mismatch; shaded = base quality below IGV's full-colour quality (20)
    """
    txt = Path(path).read_text()
    token = re.compile(r"<g\b([^>]*)>|</g\s*>|<(polygon|rect|line)\b([^>]*?)/?>")
    attr = re.compile(r'([\w:-]+)="([^"]*)"')
    stack, polys, rects, lines = [], [], [], []
    for m in token.finditer(txt):
        if m.group(0).startswith("</g"):
            if stack:
                stack.pop()
            continue
        if m.group(2) is None:
            a = dict(attr.findall(m.group(1)))
            up = stack[-1] if stack else {"fill": None, "op": 1.0, "tx": 0.0, "ty": 0.0}
            g = dict(up)
            if "fill" in a:
                g["fill"] = a["fill"]
            if "fill-opacity" in a:
                g["op"] = float(a["fill-opacity"])
            t = re.search(r"translate\(([-\d.]+),([-\d.]+)\)", a.get("transform", ""))
            if t:
                g["tx"], g["ty"] = up["tx"] + float(t.group(1)), up["ty"] + float(t.group(2))
            stack.append(g)
            continue
        if not stack:
            continue
        g, a = stack[-1], dict(attr.findall(m.group(3)))
        fill = a.get("fill", g["fill"])
        op = float(a["fill-opacity"]) if "fill-opacity" in a else g["op"]
        if m.group(2) == "polygon":
            pts = [float(v) for v in a.get("points", "").split()]
            if len(pts) == 12 and fill != "none":
                polys.append((g["tx"], g["ty"], pts, colour(fill)))
        elif m.group(2) == "rect":
            rects.append((g["tx"], g["ty"], float(a.get("x", 0)), float(a.get("y", 0)), float(a.get("width", 0)),
                          float(a.get("height", 0)), colour(fill), op))
        elif a.get("y1") == a.get("y2"):
            lines.append((g["tx"], g["ty"], float(a["x1"]), float(a["x2"]), float(a["y1"])))
    if not polys:
        raise RuntimeError(f"no alignments in {path}")
    ty = min(p[1] for p in polys)                       # the tumour panel is the upper alignment panel
    tx = next(p[0] for p in polys if p[1] == ty)
    # pixels per base from the reference-sequence track (one rectangle per base, lowest panel)
    seq_ty = max(r[1] for r in rects if r[5] == 10 and r[6] in ((0, 150, 0), (0, 0, 255), (209, 113, 5), (255, 0, 0)))
    xs = sorted({r[2] for r in rects if r[1] == seq_ty and r[5] == 10 and r[3] == 2})
    scale = (xs[-1] - xs[0]) / (len(xs) - 1)
    blocks, gaps, marks = [], [], []
    for px, py, pts, col in polys:
        if py != ty:
            continue
        ys = pts[1::2]
        top, mid = min(ys), (min(ys) + max(ys)) / 2
        if max(ys) - top != 12:
            continue
        body = [pts[i] for i in range(0, 12, 2) if pts[i + 1] == top]
        left, right = min(body), max(body)
        tips = [pts[i] for i in range(0, 12, 2) if pts[i + 1] == mid]
        strand = "r" if min(tips) < left else ("f" if max(tips) > right else None)
        blocks.append((top, left, right, col, strand))
    for lx, ly, x1, x2, y in lines:
        if ly == ty and lx == tx and y >= 44:
            gaps.append((y - 6, min(x1, x2), max(x1, x2)))
    for rx, ry, x, y, w, h, col, op in rects:
        if ry == ty and rx == tx and h == 12 and y >= 44 and col in BASE_COLOUR:
            marks.append((y, x, BASE_COLOUR[col], op < 0.999))
    return blocks, gaps, marks, scale


def reads_of(path):
    """Reads of the tumour panel: dicts with start/end (base index in the view), strand, colour,
    mismatches {index: (base, shaded)} and deletions [(index, length)]."""
    blocks, gaps, marks, scale = parse_svg(path)
    idx = lambda x: int(round(x / scale))       # noqa: E731
    rows = {}
    for b in blocks:
        rows.setdefault(b[0], []).append(b)
    reads = []
    for row, bl in rows.items():
        bl.sort(key=lambda b: b[1])
        row_gaps = [g for g in gaps if g[0] == row]
        cur = None
        for top, left, right, col, strand in bl:
            joined = cur is not None and any(abs(g[1] - cur["px_end"]) <= 2 and abs(g[2] - left) <= 2 for g in row_gaps)
            if joined:
                cur["dels"].append((idx(cur["px_end"]), idx(left) - idx(cur["px_end"])))
                cur["blocks"].append((left, right))
                cur["px_end"], cur["end"] = right, idx(right)
                cur["strand"] = cur["strand"] or strand
            else:
                cur = {"row": row, "start": idx(left), "end": idx(right), "px_end": right, "strand": strand, "colour": col,
                       "blocks": [(left, right)], "dels": [], "mm": {}}
                reads.append(cur)
        for y, x, base, shaded in marks:
            if y != row:
                continue
            for rd in reads:
                if rd["row"] != row:
                    continue
                for left, right in rd["blocks"]:
                    if left - 1 <= x < right:
                        rd["mm"][idx(left) + int(round((x - left) / scale))] = (base, shaded)
    return reads


# ------------------------------------------------------------------ review

def left_normalise(ref, gap, d):
    while gap > 0 and ref[gap - 1] == ref[gap + d - 1]:
        gap -= 1
    return gap


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts-only", action="store_true")
    a = ap.parse_args()
    WORK.mkdir(parents=True, exist_ok=True)
    (WORK / "svg").mkdir(exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    lines = (INP / "reference" / "ref.fa").read_text().splitlines()
    contig, ref = lines[0][1:].split()[0], "".join(lines[1:]).upper()
    cands = []
    for line in open(INP / "candidates.vcf"):
        if not line.startswith("#") and line.strip():
            f = line.rstrip("\n").split("\t")
            cands.append({"id": f[2], "pos": int(f[1]), "ref": f[3], "alt": f[4]})

    strands = [] if a.counts_only else ["--strands", "read"]
    t0 = count(INP / "cram" / "tumor.cram", "tumor", *strands)
    n0 = count(INP / "cram" / "normal.cram", "normal", *strands)
    t20 = {} if a.counts_only else count(INP / "cram" / "tumor.cram", "tumor_q20", "--strands", "read", "--minMapQuality", "20")
    log("igvtools count done")

    for c in cands:
        p, d = c["pos"], len(c["ref"]) - len(c["alt"])
        c["snv"] = d == 0

        def alt(cnt, c=c, p=p, d=d):
            if c["snv"]:
                f, r = cnt.get(p, ({}, {}))
                return f.get(c["alt"], 0), r.get(c["alt"], 0)
            # deletion: every placement of the gap that gives the same sequence counts (the reads may draw
            # it anywhere inside a repeat); each supporting read adds d to the DEL column over that stretch
            g = p                           # 0-based start of the gap in the candidate's (left-aligned) form
            while ref[g] == ref[g + d]:
                g += 1
            tot = [sum(cnt.get(q, ({}, {}))[s].get("DEL", 0) for q in range(p + 1, g + d + 1)) for s in (0, 1)]
            return int(round(tot[0] / d)), int(round(tot[1] / d))
        c["af"], c["ar"] = alt(t0)
        c["n"] = c["af"] + c["ar"]
        c["depth"] = depth(t0, p)
        c["n20"] = sum(alt(t20)) if t20 else c["n"]
        c["n_alt"], c["n_depth"] = sum(alt(n0)), depth(n0, p)
        f, r = t0.get(p, ({}, {}))
        c["cover"] = (sum(f.get(k, 0) for k in ("A", "C", "G", "T", "N", "DEL")),
                      sum(r.get(k, 0) for k in ("A", "C", "G", "T", "N", "DEL")))

    igv = None
    if not a.counts_only:
        igv = start_igv()
        for cmd in ("preference SAM.COMPUTE_ISIZES false", "preference SAM.INSERT_SIZE_THRESHOLD 1000",
                    "colorBy UNEXPECTED_PAIR", "expand", "maxPanelHeight 20000", f"snapshotDirectory {WORK / 'svg'}"):
            igv(cmd)
        for c in cands:
            lo = c["pos"] - HALF
            igv(f"goto {contig}:{lo}-{c['pos'] + HALF}")
            igv(f"snapshot {c['id']}.svg")
            c["reads"], c["origin"] = reads_of(WORK / "svg" / f"{c['id']}.svg"), lo
        log("SVG views written and parsed")

    def supporters(c, pos=None, alt=None, ref_allele=None, view=None):
        """Reads of the SVG view `view` (default: the candidate's own) that support an allele."""
        view = view or c
        pos, alt, ref_allele = pos or c["pos"], alt or c["alt"], ref_allele or c["ref"]
        i, d = pos - view["origin"], len(ref_allele) - len(alt)
        out = []
        for rd in view["reads"]:
            if not rd["start"] <= i < rd["end"]:
                continue
            if d == 0:
                if rd["mm"].get(i, (None,))[0] != alt:
                    continue
                out.append({"read": rd, "low_bq": rd["mm"][i][1], "end": min(i - rd["start"], rd["end"] - 1 - i),
                            "other": len(rd["mm"]) - 1})
            else:
                for gi, gl in rd["dels"]:
                    if gl == d and left_normalise(ref, view["origin"] - 1 + gi, d) == pos:
                        q = gi - rd["start"]
                        out.append({"read": rd, "low_bq": False, "end": min(q, 100 - q), "other": len(rd["mm"])})
                        break
        return out

    for k, c in enumerate(cands):
        n = c["n"]
        if c["n_alt"] >= 2 and c["n_depth"] and c["n_alt"] / c["n_depth"] >= 0.05:
            c["tag"] = "GERMLINE"
        elif n < 4:
            c["tag"] = "LOW_SUPPORT"
        elif a.counts_only:
            c["tag"] = "PASS"
        else:
            sup = supporters(c)
            if len(sup) != n:
                log(f"{c['id']}: {len(sup)} supporting reads in the image, {n} counted by igvtools")
            m = max(1, len(sup))
            if sum(s["low_bq"] for s in sup) >= 0.75 * m:
                c["tag"] = "LOW_BASE_QUALITY"
            elif n - c["n20"] >= 0.75 * n:
                c["tag"] = "LOW_MAPPING_QUALITY"
            elif sum(s["other"] for s in sup) >= 4 * m:
                c["tag"] = "MULTI_MISMATCH"
            elif max(c["af"], c["ar"]) >= 0.9 * n and c["cover"][0] > 0 and c["cover"][1] > 0:
                c["tag"] = "STRAND_BIAS"
            elif sum(s["end"] < 5 for s in sup) >= 0.75 * m:
                c["tag"] = "READ_END"
            elif sum(s["read"]["colour"] == LONG_INSERT for s in sup) >= 0.75 * m:
                c["tag"] = "DISCORDANT_PAIRS"
            else:
                c["tag"] = "PASS"
    if not a.counts_only:           # multi-nucleotide events: the two alleles on the same reads of one view
        for left, right in zip(cands, cands[1:]):
            if not (left["snv"] and right["snv"] and right["pos"] == left["pos"] + 1):
                continue
            if left["tag"] != "PASS" or right["tag"] != "PASS":
                continue
            sa = {id(s["read"]) for s in supporters(left)}
            sb = {id(s["read"]) for s in supporters(left, pos=right["pos"], alt=right["alt"], ref_allele=right["ref"])}
            if sa | sb and len(sa & sb) >= 0.9 * len(sa | sb):
                right["tag"] = "MNV_PARTNER"

    rows = ["variant_id\tcall\ttag\ttumor_alt\ttumor_depth\ttumor_vaf\tnormal_alt\tnormal_depth"]
    for c in cands:
        vaf = c["n"] / c["depth"] if c["depth"] else 0.0
        rows.append(f"{c['id']}\t{'somatic' if c['tag'] == 'PASS' else 'fail'}\t{c['tag']}\t{c['n']}\t{c['depth']}\t"
                    f"{vaf:.4f}\t{c['n_alt']}\t{c['n_depth']}")
    (OUT / "review.tsv").write_text("\n".join(rows) + "\n")
    log(f"reviewed {len(cands)} candidates: {sum(1 for c in cands if c['tag'] == 'PASS')} somatic")
    if a.counts_only:
        return

    # documentation through the batch port
    (OUT / "snapshots").mkdir(exist_ok=True)
    for cmd in ("colorBy READ_STRAND", "maxPanelHeight 1000", f"snapshotDirectory {OUT / 'snapshots'}"):
        igv(cmd)
    for c in cands:
        igv(f"goto {contig}:{c['pos'] - 40}-{c['pos'] + 40}")
        igv(f"snapshot {c['id']}.png")
    igv(f"goto {contig}:{cands[0]['pos'] - 40}-{cands[0]['pos'] + 40}")
    igv(f"saveSession {OUT / 'review_session.xml'}")
    log("session and images written")


if __name__ == "__main__":
    main()
