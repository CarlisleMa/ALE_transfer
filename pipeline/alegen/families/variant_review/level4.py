"""Level 4 of `variant_review`: review on a locked-down workstation where IGV is the only reader.

Elements (AUTHORING.md "Level 4"):
  T  the alignments are CRAM 3.1 and the image holds no reader but IGV 2.19.8 (GUI, batch port,
     igvtools); the Java runtime is cut down to what IGV needs and the bundled libraries lose
     every program entry point except IGV's own.
  P  three situations in which counting alleles gives the wrong call:
       mnv            two adjacent substitutions; on the same molecules (one event, the right-hand
                      candidate is merged) or on disjoint molecules with identical per-strand counts
       long_insert    support confined to read pairs with an insert size of 1.2-1.9 kb
       normal_trace   a true variant with one stray alternate molecule in the normal
     plus mixed artefacts: every artefact class has clean supporters and carries other blemishes
     below their thresholds, true variants carry a few low-quality supporters.

Everything of levels 1-3 is untouched; `generate()` in `__init__.py` calls into this module only
for `level >= 4`.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import struct
import zlib
from pathlib import Path

import numpy as np

from ... import desktop, dockerlib
from ...spec import Attack, Probe, TaskBuild, Trap

ASSETS = Path(__file__).parent / "assets"
PASS_MARK = 0.9


# --------------------------------------------------------------------------- planting

def plan_extras(cands: list[dict], rng) -> None:
    """Level-4 decorations of the candidate plan (consumes the generator's RNG at level 4 only)."""
    som = [c for c in cands if c["cls"] == "somatic_snv"]
    for i in rng.permutation(len(som))[:4]:
        som[int(i)]["normal_trace"] = True


def mix(frags: list, cands: list[dict], rng, Frag, RL: int) -> None:
    """Make the tumour evidence mixed instead of textbook-clean, staying far from every SOP threshold.
    Called after the mechanisms are planted and before unrelated duplicates are added."""
    low_q = lambda: int(rng.integers(5, 13))        # noqa: E731
    low_mq = lambda: int(rng.integers(0, 10))       # noqa: E731
    for c in cands:
        if "del_len" in c:
            continue
        x, cls = c["pos0"], c["cls"]
        sup = [f for f in frags if not f.dup and x in f.subs and f.covering(x)]
        sup = [sup[i] for i in rng.permutation(len(sup))]
        k = len(sup)
        if cls in ("somatic_snv", "mnv"):
            if cls == "somatic_snv" or (c["phase"] == "cis" and c["role"] == "L"):
                # keep true variants clearly two-stranded (the SOP's strand rule must not fire by chance)
                x2 = x + 1 if cls == "mnv" else None
                while True:
                    nf = sum(1 for f in sup if f.covering(x) == "f")
                    nr = len(sup) - nf
                    if nf + nr >= 8 and max(nf, nr) <= 0.72 * (nf + nr):     # and never borderline in number
                        break
                    want = "r" if nf > nr else "f"
                    pool = [f for f in frags if not f.dup and x not in f.subs and f.covering(x) == want
                            and (x2 is None or (f.covering(x2) == want and x2 not in f.subs))]
                    f = pool[int(rng.integers(0, len(pool)))]
                    f.subs[x] = (c["alt"], 37)
                    if x2 is not None:
                        f.subs[x2] = (c["mate_alt"], 37)
                    sup.append(f)
                k = len(sup)
            # true variants are not spotless: a few low-quality or poorly mapped supporters
            nb = int(rng.integers(1, k // 6 + 2))
            for f in sup[:nb]:
                f.subs[x] = (f.subs[x][0], low_q())
            if c.get("role") != "R":
                for f in sup[nb:nb + int(rng.integers(0, k // 8 + 1))]:
                    f.mapq = low_mq()
        elif cls == "low_baseq":
            sup[0].subs[x] = (c["alt"], 37)                     # one good-quality supporter
        elif cls == "low_mapq":
            sup[0].mapq = 60                                    # one well-mapped supporter
            sup[1].subs[x] = (c["alt"], low_q())
        elif cls == "multi_mismatch":
            for f in sup:
                f.xmm += 1
            sup[0].xmm = 0                                      # one clean read among the paralog-like ones
        elif cls == "strand_bias":
            for f in sup[:2]:
                f.subs[x] = (c["alt"], low_q())
            sup[2].mapq = low_mq()
        elif cls == "read_end":
            for f in sup[:2]:
                f.mapq = low_mq()
            flen = int(np.clip(rng.normal(340, 30), 232, 480))  # one supporter with the base mid-read
            f = Frag(x - int(rng.integers(30, 70)), flen, int(rng.integers(0, 2)), bool(rng.integers(0, 2)),
                     f"T:m{c['id'][-3:]}")
            f.subs[x] = (c["alt"], 37)
            frags.append(f)
        elif cls == "long_insert":
            for f in sup[:2]:
                f.subs[x] = (c["alt"], low_q())


# --------------------------------------------------------------------------- helpers

def _asset_module(name: str):
    spec = importlib.util.spec_from_file_location(f"_vr_{name}", ASSETS / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _write_cram(bam: Path, cram: Path, ref_fa: Path, ref_seq: str) -> None:
    """BAM -> CRAM 3.1 whose bytes do not depend on the build host.

    htslib records two host paths in a CRAM: the output path as the 20-byte file id (patched after
    writing) and the absolute path of the reference as @SQ UR (always rewritten, whatever the
    header says). The reference is therefore opened through /proc/self/cwd, which is the same
    string on every host; readers that take the reference by name or MD5 (IGV does) ignore UR."""
    import os

    import pysam
    cram.parent.mkdir(parents=True, exist_ok=True)
    here = os.getcwd()
    try:
        os.chdir(ref_fa.parent)
        with pysam.AlignmentFile(str(bam.resolve())) as b:
            hd = b.header.to_dict()
            hd["SQ"][0]["M5"] = hashlib.md5(ref_seq.encode()).hexdigest()
            with pysam.AlignmentFile(str(cram.resolve()), "wc", header=hd,
                                     reference_filename=f"/proc/self/cwd/{ref_fa.name}",
                                     format_options=[b"version=3.1"]) as out:
                for r in b:
                    out.write(r)
    finally:
        os.chdir(here)
    raw = bytearray(cram.read_bytes())
    if raw[:4] != b"CRAM" or raw[4:6] != b"\x03\x01":
        raise RuntimeError("unexpected CRAM version")
    raw[6:26] = cram.name.encode()[:20].ljust(20, b"\0")
    cram.write_bytes(bytes(raw))
    pysam.index(str(cram))


def _png(width: int, height: int, rows: list[bytes]) -> bytes:
    def chunk(t, d):
        return struct.pack(">I", len(d)) + t + d + struct.pack(">I", zlib.crc32(t + d))
    raw = b"".join(b"\x00" + r for r in rows)
    return (b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0))
            + chunk(b"IDAT", zlib.compress(raw, 6)) + chunk(b"IEND", b""))


def _mock_snapshot(key: str) -> bytes:
    """A structurally valid stand-in for an IGV snapshot (battery only): rows of read-like bars."""
    rng = np.random.default_rng(int(hashlib.sha256(key.encode()).hexdigest()[:12], 16))
    w, h = 1148, 520
    img = np.full((h, w, 3), 250, dtype=np.uint8)
    for row in range(60):
        y = 20 + row * 8
        s = int(rng.integers(0, 300))
        while s < w - 40:
            e = min(w, s + int(rng.integers(300, 480)))
            img[y:y + 6, s:e] = (185, 185, 185)
            if rng.random() < 0.5:
                img[y:y + 6, s + 20:s + 24] = (255, 0, 0)
            s = e + int(rng.integers(8, 60))
    return _png(w, h, [img[i].tobytes() for i in range(h)])


def _blank_snapshot() -> bytes:
    return _png(1148, 800, [bytes([250]) * (1148 * 3)] * 800)


def _mock_session(contig: str, first_pos1: int) -> str:
    """What IGV 2.19.8 writes for this task (structure copied from a session saved by IGV)."""
    def panel(sm):
        p = f"/task/input/cram/{sm}.cram"
        return (f'    <Panel height="1736" name="Panel_{sm}" width="1129">\n'
                f'        <Track attributeKey="{sm}.cram Coverage" autoScale="true" clazz="org.broad.igv.sam.CoverageTrack" '
                f'fontSize="10" id="{p}_coverage" name="{sm}.cram Coverage" snpThreshold="0.2" visible="true"/>\n'
                f'        <Track attributeKey="{sm}.cram" clazz="org.broad.igv.sam.AlignmentTrack" color="185,185,185" '
                f'displayMode="EXPANDED" experimentType="OTHER" fontSize="10" id="{p}" name="{sm}.cram" visible="true">\n'
                '            <RenderOptions colorOption="READ_STRAND" duplicatesOption="FILTER"/>\n'
                '        </Track>\n    </Panel>\n')
    return ('<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n'
            f'<Session genome="/task/input/reference/ref.fa" locus="{contig}:{first_pos1 - 41}-{first_pos1 + 40}" version="8">\n'
            '    <Resources>\n        <Resource path="/task/input/cram/tumor.cram" type="cram"/>\n'
            '        <Resource path="/task/input/cram/normal.cram" type="cram"/>\n'
            '        <Resource path="/task/input/candidates.bed" type="bed"/>\n    </Resources>\n'
            + panel("tumor") + panel("normal") +
            '    <Panel height="99" name="FeaturePanel" width="1129">\n'
            '        <Track attributeKey="Reference sequence" clazz="org.broad.igv.track.SequenceTrack" fontSize="10" '
            'id="Reference sequence" name="Reference sequence" visible="true"/>\n'
            '        <Track attributeKey="candidates.bed" clazz="org.broad.igv.track.FeatureTrack" fontSize="10" '
            'id="/task/input/candidates.bed" name="candidates.bed" visible="true"/>\n    </Panel>\n'
            '    <PanelLayout dividerFractions="0.0106,0.4078,0.8156"/>\n</Session>\n')


def _handwritten_session(contig: str) -> str:
    """The level-3 style of session: a resource list typed by hand."""
    return ('<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n'
            f'<Session genome="/task/input/reference/ref.fa" locus="{contig}:1-2000" version="8">\n <Resources>\n'
            '  <Resource path="/task/input/cram/tumor.cram"/>\n  <Resource path="/task/input/cram/normal.cram"/>\n'
            '  <Resource path="/task/input/candidates.bed"/>\n </Resources>\n</Session>\n')


# --------------------------------------------------------------------------- texts

INSTRUCTION = """\
You are the variant-review analyst of a cancer genomics laboratory. An automated caller produced
somatic variant candidates for one tumour/normal pair; before reporting, every candidate must be
reviewed against the aligned reads following the laboratory's review SOP. You work on the
laboratory's locked-down review workstation.

Task directory:
- `/task`

Inputs (read-only) under `/task/input`:
- `cram/tumor.cram`, `cram/normal.cram` (+ `.crai`) - coordinate-sorted alignments (CRAM 3.1,
  reference-compressed), duplicates marked
- `reference/ref.fa` (+ `.fai`) - the reference the reads were aligned to
- `candidates.vcf` - the {n} candidates to review (`ID` column = variant id); `candidates.bed` holds
  the same loci for display
- `review_sop.md` - the laboratory's review SOP (evidence definition, failure tags and their order,
  documentation of the review)

Available software. **IGV 2.19.8 is the only software on this workstation that can read the
alignment files.** There is no samtools, no pysam or htslib, and no other CRAM/BAM reader; Python
3.12 is available with the standard library and numpy. IGV can be used in three ways (reference:
`/task/software/igv/README.md`):
- the IGV desktop application - `/task/software/open_igv.sh` starts it on the virtual display `:99`
  (1280x800 pixels) with the reference, both CRAM files and the candidate track loaded;
- IGV's batch port (127.0.0.1:60151, text commands, available while the application is running);
- `igvtools` (on `PATH`).

The desktop is operated with the command `desktop` (on `PATH`; `desktop --help`):
`desktop screenshot <file.png> [--region X,Y,W,H]`,
`desktop click <x> <y> [--button left|middle|right] [--clicks N]`, `desktop move <x> <y>`,
`desktop drag <x1> <y1> <x2> <y2>`, `desktop type <text>`, `desktop key <combo> [<combo> ...]` (xdotool
names, e.g. `Return`, `Escape`, `ctrl+a`), `desktop scroll <x> <y> <clicks>`, `desktop launch <command>`,
`desktop windows`, `desktop wait <seconds>`.
If your harness connects MCP servers, the server `desktop` offers the same actions as tools
(`screenshot`, `click`, `move`, `drag`, `type_text`, `key`, `scroll`, `launch`, `list_windows`, `wait`).
There is no internet access.

Your task: review every candidate according to `review_sop.md` and document the review as the SOP
requires.

Required outputs under `/task/output`:
- `review.tsv` - tab-separated, one row per candidate, columns
  `variant_id`, `call`, `tag`, `tumor_alt`, `tumor_depth`, `tumor_vaf`, `normal_alt`, `normal_depth`.
  `call` is `somatic` or `fail`; `tag` is `PASS` for a somatic call, otherwise exactly one failure
  tag from the SOP; the count columns are in molecules as defined in the SOP and
  `tumor_vaf = tumor_alt / tumor_depth`.
- `review_session.xml` - the IGV session of the review, saved by IGV as specified in the SOP.
- `snapshots/<variant_id>.png` - one IGV image per candidate as specified in the SOP.

Constraints:
- Do not modify anything under `/task/input`; write deliverables only under `/task/output`.
- Do not ask for confirmation. Execute directly.
"""

SOP = """\
# Somatic variant review SOP (laboratory document R-07, revision 4)

Review each candidate in the tumour and in the matched normal. The reviewer records one call per
candidate: `somatic`, or `fail` with exactly one failure tag.

## What counts as evidence
- Evidence is counted in sequenced molecules (read pairs). Reads marked as duplicates add no
  evidence and no depth. Secondary and supplementary alignments are ignored.
- Depth at a candidate = molecules whose read covers the candidate position.
- A molecule supports a single-nucleotide candidate if its read shows the alternate base at that
  position, whatever the base quality or mapping quality (those are judged by the tags below).
- For insertions and deletions, judge what the read sequence shows, not where the aligner drew the
  gap: a read supports the candidate if it carries a gap that produces the same sequence as the
  candidate allele.
- "Supporting reads" are the reads of the supporting molecules that cover the candidate.

## Failure tags, in order of precedence (record the first that applies)
1. `GERMLINE` - the alternate allele is seen in the normal in at least 2 molecules and at least 5%
   of the normal depth.
2. `LOW_SUPPORT` - fewer than 4 supporting molecules in the tumour.
3. `LOW_BASE_QUALITY` - at least three quarters of the supporting molecules show the alternate
   base with base quality below 20.
4. `LOW_MAPPING_QUALITY` - at least three quarters of the supporting molecules have mapping
   quality below 20.
5. `MULTI_MISMATCH` - the supporting reads carry, on average, four or more further mismatches to
   the reference within the same read (aligned bases, not counting the candidate itself).
6. `STRAND_BIAS` - at least 90% of the supporting reads are on one strand although reads on both
   strands cover the position.
7. `READ_END` - at least three quarters of the supporting reads show the alternate base within 5
   bases of either end of the read.
8. `DISCORDANT_PAIRS` - at least three quarters of the supporting molecules are read pairs with an
   insert size (absolute template length) above 1,000 bp. (The library was prepared with fragments
   of 200-500 bp.)
9. `MNV_PARTNER` - concerns two single-nucleotide candidates at directly adjacent positions to
   neither of which a tag above applies. If they are carried by the same molecules (at least 90% of
   the molecules that support either candidate support both), they are one multi-nucleotide
   variant: it is reported on the left-hand candidate, which stays `somatic`, and the right-hand
   candidate is recorded as `fail` with this tag. Adjacent candidates carried by different
   molecules are separate variants and are both kept.

A candidate to which no tag applies is called `somatic` (tag `PASS`).

## Documentation of the review
The review is documented from IGV; both items below must be written by IGV itself.
- Session (`review_session.xml`): saved from IGV with the reference `/task/input/reference/ref.fa`,
  both alignment files and the candidate track loaded, the alignments of both alignment tracks
  coloured by read strand, and the view on the first candidate of the list (the displayed range
  contains its position and is at most 200 bp wide).
- Images (`snapshots/<variant_id>.png`): for every candidate one PNG image saved from IGV with the
  view on that candidate (displayed range containing the candidate, at most 200 bp wide) and the
  alignments visible. Each candidate has its own image, at least 800 pixels wide and 300 pixels
  high; an image that shows no data, or the same image filed for two candidates, does not count.
"""

IGV_README = """\
# IGV 2.19.8 on the review workstation

## Desktop application
`/task/software/open_igv.sh` starts IGV on display `:99` (1280x800) with the reference, both CRAM
files and the candidate track (log: `/tmp/igv.log`) and returns once these are loaded, usually
after 10-40 seconds (IGV opens its window and its batch port a few seconds earlier, while it is
still loading). On its first start the script writes the workstation defaults to
`~/igv/prefs.properties` (offline start-up, read down-sampling switched off) and sets IGV's Java
heap in `~/.igv/java_arguments` (2 GB).

## Batch port
While the application is running it listens on `127.0.0.1:60151`. Send one command per line (plain
text, newline-terminated); IGV answers every command with one line (`OK`, or an error text) once
the command has finished. The commands are those of IGV batch scripts, among them:

    echo                               answers "echo"
    new                                unload all tracks
    genome <id|path>                   load a reference genome
    load <path>[,<path>...]            load data files
    goto <locus>                       e.g. goto chr1:1000-1200
    snapshotDirectory <dir>            directory for images
    snapshot [<file>]                  save the IGV view as an image (.png or .svg, by extension)
    saveSession <file.xml>             save the session
    maxPanelHeight <pixels>            maximum height of a track panel in an image
    collapse | expand | squish         display mode of the tracks
    viewaspairs                        draw alignments as pairs
    sort <option> [<locus>]            sort alignments, e.g. base, strand, quality, insertSize
    group <option> [<position>]        group alignments, e.g. strand, base, readGroup
    colorBy <option>                   colour alignments, e.g. READ_STRAND, FIRST_OF_PAIR_STRAND,
                                       INSERT_SIZE, UNEXPECTED_PAIR, NONE
    preference <key> <value>           set an IGV preference (keys as in prefs.properties)
    setSleepInterval <milliseconds>    pause between commands
    exit                               quit IGV

## igvtools
`igvtools <command> ...` runs IGV's command-line utilities (Java heap 1 GB).

    igvtools count [options] <input> <output> <genome>
        coverage of an alignment file. <output>: a .wig file (text), a .tdf file, or both separated
        by a comma. <genome>: a FASTA file, a .chrom.sizes file or a genome id. Options:
          -w, --windowSize <n>       window size in bp (default 25)
          -z, --maxZoom <n>          zoom levels precomputed for .tdf output (default 7)
          -e, --extFactor <n>        extend reads by n bp
          --strands read|first       count the two strands separately (strand of the read / of the
                                     first read of the pair)
          --bases                    count A, C, G, T, N, deletions and insertions separately
          --query <chr:start-end>    restrict to a region (indexed input)
          --minMapQuality <q>        skip alignments below this mapping quality (default 0)
          --includeDuplicates        include alignments flagged as duplicates
    igvtools index <file>            index a feature file
    igvtools toTDF [options] <input> <output.tdf> <genome>
    igvtools version
    (`igvtools sort` is not available for alignment files on this workstation.)
"""

JUDGE = """\
**Where the truth comes from.** The alignments are synthesised read by read, so the generator knows
what every molecule carries. Each candidate gets one planted mechanism with parameters far from
every SOP threshold (the generator re-derives the SOP outcome from the reads it wrote and refuses
the seed if any quantity is within its safety margin of a threshold): the eleven level-3 mechanisms,
plus pairs of adjacent substitutions carried by the same molecules (one event) or by disjoint
molecules with identical per-strand allele counts (two events), variants supported only by read
pairs with a 1.2-1.9 kb insert, and true variants with one stray alternate molecule in the normal.
Artefacts are mixed (each has clean supporters and carries other blemishes below their thresholds)
and true variants carry a few low-quality supporters. `reference/truth.json` stores the planted
call and tag and the realised molecule-level allele fraction.

**Hard gates (score 0).** `review.tsv` missing/unparseable; a candidate missing or listed twice;
`call` not in {somatic, fail}; `tag` outside the SOP vocabulary; `PASS` used with `fail` or a
failure tag with `somatic`; non-numeric count columns.

**Components.**

| weight | component | metric | full credit | zero |
|---|---|---|---|---|
| 0.40 | calls | Matthews correlation of `call` with the planted truth | >= 0.965 | <= 0.86 |
| 0.30 | failure_tags | fraction of planted failures given the planted tag | >= 0.97 | <= 0.85 |
| 0.15 | allele_fraction | fraction of planted somatic variants called `somatic` with `tumor_vaf` within 0.05 of the realised molecule-level fraction | >= 0.96 | <= 0.80 |
| 0.07 | igv_session | session written by IGV (track panels present) with the reference, the three tracks, both alignment tracks coloured by read strand, view on the first candidate (<= 200 bp) | all five | not an IGV session |
| 0.08 | snapshots | fraction of candidates with their own PNG (>= 800x300, shows data, not filed twice) | all | none |

The review components carry 0.85. Session and images cannot lift a failing review over the pass
mark: if the review components reach less than 0.9 of their weight, the total is capped at 0.89.

**Why these bands.** The truth is exact, so a correct review reproduces every call and tag; the
bands give graded reward, they do not absorb noise. With 68 candidates one wrong call gives MCC
0.970 and full credit; ignoring a class of four candidates gives MCC 0.89 (calls 0.27 of their
weight) and tag accuracy 0.895 (0.37). The allele-fraction tolerance covers legitimate counting
conventions (depth at the anchor base versus across a repeat) but not position-exact indel counts.

**What cannot be verified.** That IGV wrote the session and the images: the verifier checks
structure (IGV's track panels and render options; PNG size, content, distinctness), which a careful
imitation satisfies. Both together are worth 0.15 and cannot decide a failing review. The review
table itself cannot be produced without reading the alignments.
"""


def dockerfile() -> str:
    """Agent image of level 4. Own Dockerfile string (the core helper always installs scipy/pandas
    and cannot post-process the conda environment in the same layer)."""
    harden = (ASSETS / "harden_igv.py").read_text()
    if "\nFROM " in harden or "ALEGEN_" in harden:
        raise RuntimeError("harden_igv.py would break the Dockerfile")
    out = dockerlib.AGENT_BASE
    out += ("# IGV 2.19.8 is the only alignment reader of this image. The Java runtime is re-linked with the\n"
            "# modules IGV needs (no compiler, no jshell, no scripting API) and IGV's bundled libraries lose\n"
            "# their JavaScript engine, the alignment sorter of igvtools and every other program entry point.\n"
            "RUN <<'ALEGEN_EOS'\n"
            "set -e\n"
            "micromamba create -y -p /opt/tools -c conda-forge -c bioconda python=3.12 numpy igv=2.19.8\n"
            "cat > /tmp/harden_igv.py <<'ALEGEN_PY'\n"
            + harden.rstrip("\n") + "\n"
            "ALEGEN_PY\n"
            "python3 /tmp/harden_igv.py\n"
            "rm -f /tmp/harden_igv.py\n"
            "micromamba clean -a -y\n"
            "ALEGEN_EOS\n"
            "ENV PATH=/opt/tools/bin:$PATH\n")
    out += "RUN pip install " + " ".join(desktop.PIP) + "\n"
    out += desktop.APT
    out += desktop.CLI_LINK
    out += "RUN ln -s /task/software/igv/igvtools /usr/local/bin/igvtools\n"
    out += dockerlib.AGENT_TAIL
    return out


OPEN_IGV = """#!/bin/bash
# Start IGV on the virtual display with the reference, both CRAM files and the candidate track.
export DISPLAY="${DISPLAY:-:99}"
bash /task/software/desktop/start_desktop.sh
if pgrep -f 'org.broad.igv.ui.Main' >/dev/null 2>&1; then
  echo "IGV is already running (log: /tmp/igv.log); batch port 60151"
  exit 0
fi
mkdir -p "$HOME/igv" "$HOME/.igv"
if [ ! -s "$HOME/igv/prefs.properties" ]; then
  printf 'PROVISIONING_URL_DEFAULT=/task/software/igv/oauth_offline.json\\nIGV.genome.sequence.dir=/task/software/igv/genomes.tsv\\nCHECK_VERSION=false\\nSAM.DOWNSAMPLE_READS=false\\n' > "$HOME/igv/prefs.properties"
fi
[ -s "$HOME/.igv/java_arguments" ] || echo "-Xmx2g" > "$HOME/.igv/java_arguments"
nohup igv -g /task/input/reference/ref.fa /task/input/cram/tumor.cram,/task/input/cram/normal.cram,/task/input/candidates.bed "$@" > /tmp/igv.log 2>&1 &
# IGV opens its window and its batch port before the reference and the tracks are loaded: wait for them.
for _ in $(seq 1 480); do
  grep -q 'Loading resource: .*candidates.bed' /tmp/igv.log 2>/dev/null && break
  sleep 0.5
done
sleep 3
if grep -q 'Loading resource: .*candidates.bed' /tmp/igv.log 2>/dev/null; then
  echo "IGV is running with the tracks loaded (log: /tmp/igv.log); batch port 60151"
else
  echo "IGV has not finished starting after 4 minutes, see /tmp/igv.log" >&2
  exit 1
fi
"""

IGVTOOLS = """#!/bin/sh
# igvtools of the IGV 2.19.8 installation (the IGV package of this workstation has no launcher for it).
exec java --module-path=/opt/tools/lib/igv -Xmx1g -Djava.awt.headless=true \\
    --module=org.igv/org.broad.igv.tools.IgvTools "$@"
"""


# --------------------------------------------------------------------------- task assembly

def build(seed: int, level: int, workdir: Path, ref: str, cands: list[dict], truth: dict, bam_dir: Path,
          vcf: str, bed: str, ref_url: str, ref_sha256: str, contig: str) -> TaskBuild:
    import pysam
    n = len(cands)
    tv = truth["variants"]
    for c in cands:     # finer class labels for the task card / analysis (not used by the verifier)
        if c["cls"] == "mnv":
            tv[c["id"]]["class"] = "mnv_" + c["phase"] + ("_left" if c["role"] == "L" else "_right")
        elif c.get("normal_trace"):
            tv[c["id"]]["class"] = "somatic_snv_normal_trace"
    truth["first_candidate"] = {"contig": contig, "pos": cands[0]["pos0"] + 1}

    fasta = f">{contig}\n" + "\n".join(ref[i:i + 60] for i in range(0, len(ref), 60)) + "\n"
    (workdir / "ref.fa").write_text(fasta)
    pysam.faidx(str(workdir / "ref.fa"))

    # ---- the planted truth must be what the SOP gives on the reads, with a safety margin ----------
    sr = _asset_module("sop_review")
    reads = {}
    for sm in ("tumor", "normal"):
        with pysam.AlignmentFile(str(bam_dir / f"{sm}.bam")) as b:
            reads[sm] = sorted((sr.Read(r.to_string().split("\t")) for r in b), key=lambda r: r.start)
    c4 = [(c["id"], c["pos0"], c["ref"], c["alt"]) for c in cands]

    def run(keep_duplicates=False, exact_indel=False, **opt):
        t = sr.Sample(reads["tumor"], ref, keep_duplicates, exact_indel)
        nn = sr.Sample(reads["normal"], ref, keep_duplicates, exact_indel)
        return sr.review(c4, ref, t, nn, opt)

    rows, res = run()
    bad = [(x["id"], tv[x["id"]]["class"], tv[x["id"]]["tag"], x["tag"]) for x in res if x["tag"] != tv[x["id"]]["tag"]]
    if bad:
        raise RuntimeError(f"planted truth and SOP outcome differ, change the seed: {bad[:5]}")
    for x in res:
        t = tv[x["id"]]
        if (len(x["sup"]), x["depth"], x["n_alt"], x["n_depth"]) != (t["tumor_alt"], t["tumor_depth"], t["normal_alt"], t["normal_depth"]):
            raise RuntimeError(f"realised counts differ for {x['id']}")
        q, why = x["frac"], []
        if 2 <= x["n_alt"] and x["n_alt"] / x["n_depth"] < 0.12:
            why.append("normal fraction")
        if x["n_alt"] >= 2 and t["tag"] != "GERMLINE":
            why.append("normal molecules")
        if t["tag"] not in ("GERMLINE", "LOW_SUPPORT"):
            if q["n"] < 6:
                why.append("support")
            for key in ("low_bq", "low_mapq", "read_end", "long_insert"):
                if 0.40 < q[key] < 0.85:
                    why.append(key)
            if 2.0 < q["other_mm"] < 5.0:
                why.append("other_mm")
            if 0.80 < q["one_strand"] < 0.97:
                why.append("one_strand")
            if x["shared"] is not None and 0.2 < x["shared"] < 0.98:
                why.append("shared")
        if why:
            raise RuntimeError(f"{x['id']} ({t['class']}) is too close to an SOP threshold ({why}); change the seed")
    for sm in ("tumor", "normal"):      # IGV's default down-sampling (100 reads per 50 bp) must never bite
        per_pos = np.bincount(np.array([r.start for r in reads[sm] if not r.flag & 0x400]), minlength=len(ref) + 50)
        if np.convolve(per_pos, np.ones(50, dtype=int)).max() > 95:
            raise RuntimeError("local read density would trigger IGV's default down-sampling; change the seed")

    # ---- CRAM ---------------------------------------------------------------------------------------
    cram_dir = workdir / "cram"
    for sm in ("tumor", "normal"):
        _write_cram(bam_dir / f"{sm}.bam", cram_dir / f"{sm}.cram", workdir / "ref.fa", ref)

    # ---- battery ------------------------------------------------------------------------------------
    ids = [c["id"] for c in cands]
    session = _mock_session(contig, cands[0]["pos0"] + 1)
    docs = {"review_session.xml": session}
    docs.update({f"snapshots/{i}.png": _mock_snapshot(f"{seed}:{i}") for i in ids})
    poor_docs = {"review_session.xml": _handwritten_session(contig)}
    poor_docs.update({f"snapshots/{i}.png": _blank_snapshot() for i in ids})

    def table(fn):
        out = [sr.HEADER]
        for i in ids:
            call, tag, vaf = fn(i)
            t = tv[i]
            out.append(f"{i}\t{call}\t{tag}\t{t['tumor_alt']}\t{t['tumor_depth']}\t{vaf:.4f}\t{t['normal_alt']}\t{t['normal_depth']}")
        return "\n".join(out) + "\n"

    def careless(**opt):
        return {"review.tsv": "\n".join(run(**opt)[0]) + "\n", **docs}

    prng = np.random.default_rng(seed + 991)
    one_wrong = ids[int(prng.integers(0, n))]
    two_wrong = [ids[int(i)] for i in prng.permutation(n)[:2]]

    def flipped(wrong, noise):
        def fn(i):
            t = tv[i]
            if i in wrong:
                return ("fail", "LOW_SUPPORT", t["tumor_vaf"]) if t["call"] == "somatic" else ("somatic", "PASS", t["tumor_vaf"])
            return t["call"], t["tag"], min(1.0, max(0.0, t["tumor_vaf"] + noise()))
        return fn

    exact = table(lambda i: (tv[i]["call"], tv[i]["tag"], tv[i]["tumor_vaf"]))
    probes = [
        Probe("planted_truth", "review written from the planted truth, session and images in IGV's format", 0.999,
              files={"review.tsv": exact, **docs}),
        Probe("within_tolerance", "one candidate reviewed wrongly, allele fractions off by ~0.01", 0.99,
              files={"review.tsv": table(flipped({one_wrong}, lambda: float(prng.normal(0, 0.008)))), **docs}),
        Probe("independent_solution", "GUI route: the IGV desktop application operated with the `desktop` command "
              "(alignment export from the track menu, session and images from the File menu); the oracle is the "
              "command-line route (igvtools count + batch port)", 0.99, script="SOLUTION=2 bash /solution/solve.sh"),
        Probe("review_credit_without_documentation", "correct table, hand-typed session, blank images: keeps the "
              "whole review weight (0.85)", 0.85, files={"review.tsv": exact, **poor_docs}),
    ]
    cap = round(PASS_MARK - 0.15, 2)
    attacks = [
        Attack("l3_reference", "the level-3 reference review (pysam) run on the level-4 workstation: no reader, no output",
               cap, script="python3 /alegen_script/review.py; true", script_files={"review.py": ASSETS / "review.py"}),
        Attack("l3_logic_on_igv_export", "level-3 SOP logic on reads exported from IGV: handles every level-3 decision "
               "point, ignores multi-nucleotide events and discordant pairs", cap,
               files=careless(no_mnv_rule=True, no_pair_rule=True)),
        Attack("igvtools_counts_only", "allele counts from `igvtools count --bases` without its duplicate, strand and "
               "mapping-quality switches: only GERMLINE and LOW_SUPPORT are judged", cap,
               script="python3 /solution/review_cli.py --counts-only"),
        Attack("ignore_mnv_phase", "adjacent candidates reviewed one by one (per-variant logic)", cap,
               files=careless(no_mnv_rule=True)),
        Attack("merge_every_adjacent_pair", "adjacent candidates merged without looking at the molecules", cap,
               files=careless(merge_adjacent=True)),
        Attack("ignore_discordant_pairs", "insert size of the supporting pairs not examined", cap,
               files=careless(no_pair_rule=True)),
        Attack("any_alt_read_in_normal_is_germline", "a single alternate molecule in the normal taken as germline", cap,
               files=careless(normal_any_alt=True)),
        Attack("counts_only_review", "only GERMLINE and LOW_SUPPORT are checked (no read-level inspection)", cap,
               files=careless(counts_only=True, no_mnv_rule=True)),
        Attack("duplicates_counted", "duplicate-flagged reads shown in IGV and counted as evidence", cap,
               files=careless(keep_duplicates=True)),
        Attack("position_exact_indels", "a deletion only counts when the gap is drawn at the candidate position", cap,
               files=careless(exact_indel=True)),
        Attack("trust_the_caller", "every candidate called somatic", 0.35,
               files={"review.tsv": table(lambda i: ("somatic", "PASS", tv[i]["tumor_vaf"])), **docs}),
        Attack("fail_everything", "every candidate failed as LOW_SUPPORT", 0.30,
               files={"review.tsv": table(lambda i: ("fail", "LOW_SUPPORT", 0.0)), **docs}),
        Attack("handwritten_session_blank_snapshots", "correct table, hand-typed session, blank images: loses the "
               "documentation weight (and exactly that, see the probe of the same content)", 0.85,
               files={"review.tsv": exact, **poor_docs}),
        Attack("two_wrong_calls", "two candidates reviewed wrongly (one more than the band absorbs)", 0.895,
               files={"review.tsv": table(flipped(set(two_wrong), lambda: 0.0)), **docs}),
        Attack("vaf_beyond_tolerance", "every allele fraction off by 0.07", 0.89,
               files={"review.tsv": table(lambda i: (tv[i]["call"], tv[i]["tag"], min(1.0, tv[i]["tumor_vaf"] + 0.07))), **docs}),
    ]

    traps = [
        Trap("T_igv_only_reader", "The alignments are CRAM 3.1 and the workstation has no reader except IGV 2.19.8 "
             "(desktop application, batch port, igvtools); the level-3 terminal route (pysam/samtools) does not exist.",
             "stated in the instruction; `python3 -c 'import pysam'` and `samtools` fail", "Approach / tool use"),
        Trap("P_mnv_same_molecule", "Four pairs of adjacent substitutions sit on the same molecules (one event, the "
             "right-hand candidate is merged), four pairs sit on disjoint molecules with identical per-strand allele "
             "counts: per-variant counting cannot tell them apart.",
             "SOP tag 9; in IGV the two coloured bases are on the same reads or on different reads; read names in "
             "an alignment export", "Understanding / per-variant instead of per-haplotype logic"),
        Trap("P_discordant_pair_support", "Four candidates are supported almost only by read pairs with a 1.2-1.9 kb "
             "insert; base counts, strands, qualities and positions in the read look like a real variant.",
             "SOP tag 8; IGV colours the reads by insert size by default (red), template length in read details "
             "and in an alignment export", "Approach / incomplete workflow"),
        Trap("P_normal_trace", "Four true variants have one stray alternate molecule in the normal (about 2% of the "
             "normal depth, below both conditions of the germline rule).",
             "SOP tag 1 states both conditions; the counts are visible in the normal track", "Execution / rule misapplied"),
        Trap("artefact_classes", "Five artefact mechanisms look like real variants in a pile-up count (strand, read "
             "end, base quality, mapping quality, multi-mismatch reads); at this level they are mixed with clean "
             "supporters and with blemishes of other kinds below their thresholds.",
             "review_sop.md defines the tags; per-read attributes in IGV (shading, strand colouring, read details) "
             "or in an alignment export", "Approach / incomplete workflow"),
        Trap("duplicate_only_support", "8-15 supporting reads on both strands that are only 2-3 molecules plus their flagged duplicates.",
             "IGV filters duplicates by default (Duplicates menu of the alignment track); `igvtools count` excludes "
             "them unless asked; SOP counts molecules", "Execution / implementation bug"),
        Trap("indel_representation", "Four true deletions in homopolymers are left-aligned in the candidate list "
             "but right-aligned in the reads; position-exact counting finds zero support.",
             "the gap is visible a few bases downstream in IGV; SOP says to judge the sequence, not the gap position",
             "Understanding / domain knowledge gap"),
    ]

    return TaskBuild(
        task_id=f"variant-review-s{seed:04d}-l{level}", family="variant_review", seed=seed, level=level,
        title="Somatic variant review on a locked-down workstation where IGV is the only alignment reader",
        summary="Review somatic variant candidates read by read in tumour/normal CRAM files following a laboratory "
                "SOP, using IGV (desktop application, batch port or igvtools) as the only software that can read the "
                "alignments, and deliver calls, failure tags, an IGV session and one IGV image per candidate.",
        subdomain="3.2 Genomics & Sequence Analysis",
        software=["IGV (GUI)", "IGV batch port", "igvtools", "desktop CLI / MCP server", "Python"],
        instruction_md=INSTRUCTION.format(n=n),
        agent_dockerfile=dockerfile(),
        verifier_dockerfile=dockerlib.verifier_dockerfile(pip=[]),
        input_files={"cram/tumor.cram": cram_dir / "tumor.cram", "cram/tumor.cram.crai": cram_dir / "tumor.cram.crai",
                     "cram/normal.cram": cram_dir / "normal.cram", "cram/normal.cram.crai": cram_dir / "normal.cram.crai",
                     "reference/ref.fa": workdir / "ref.fa", "reference/ref.fa.fai": workdir / "ref.fa.fai",
                     "candidates.vcf": vcf, "candidates.bed": bed, "review_sop.md": SOP},
        software_files={**desktop.CLI_FILES, "open_igv.sh": OPEN_IGV, "igv/igvtools": IGVTOOLS,
                        "igv/README.md": IGV_README,
                        "igv/oauth_offline.json": json.dumps({"client_id": "offline", "client_secret": "offline",
                                                              "auth_uri": "http://127.0.0.1/auth",
                                                              "token_uri": "http://127.0.0.1/token",
                                                              "project_id": "offline"}),
                        "igv/genomes.tsv": "id\tname\turl\n"},
        solution_files={"solve.sh": ASSETS / "solve_l4.sh", "sop_review.py": ASSETS / "sop_review.py",
                        "igv_gui.py": ASSETS / "igv_gui.py", "review_cli.py": ASSETS / "review_cli.py"},
        tests_files={"score.py": ASSETS / "score_l4.py", "reference/truth.json": json.dumps(truth)},
        traps=traps, attacks=attacks, probes=probes,
        capabilities=["operating professional GUI software (IGV) through desktop actions or its own command line",
                      "read-level evidence inspection", "following a multi-rule SOP with precedence",
                      "per-molecule (haplotype) reasoning", "molecule vs read counting", "indel representation",
                      "GUI-native deliverables (IGV session, snapshots)"],
        judge_explanation_md=JUDGE, pass_threshold=PASS_MARK,
        oracle_steps=14 + 4 * n, oracle_runtime_min=3.0, needs_gui=True, cpus=4, memory_mb=8000, storage_mb=10000,
        mcp_servers=[desktop.MCP_SERVER],
        extra_metadata={"reference_source": ref_url, "reference_sha256": ref_sha256,
                        "level4_elements": ["T", "P"]},
    )
