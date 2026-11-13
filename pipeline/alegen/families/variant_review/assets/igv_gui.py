#!/usr/bin/env python3
"""Reference solution of level 4, GUI route: IGV's desktop application operated with the
`desktop` command (mouse, keyboard, screenshots), nothing else.

  1. start IGV (`/task/software/open_igv.sh`)
  2. three overlapping 20-kb views; in each, the alignment track's pop-up menu
     "Export Alignments..." writes the reads on display to a SAM file (tumour and normal)
  3. apply the SOP to the exported reads (sop_review.py)
  4. colour both alignment tracks by read strand (pop-up menu), go to the first candidate,
     File > Save Session...
  5. for every candidate: go to it, File > Save PNG Image ...

Every step is checked on its result (window titles, a strip of screen pixels, the written file)
and repeated if the result is missing. Pixel positions are those of IGV 2.19.8 as opened by
open_igv.sh on the 1280x800 display.
"""
import hashlib
import re
import struct
import subprocess
import sys
import time
import zlib
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sop_review as sr  # noqa: E402

INP, OUT, WORK = Path("/task/input"), Path("/task/output"), Path("/tmp/igv_gui")
LOCUS_BOX = (345, 59)
TRACK = {"tumor": (700, 300), "normal": (700, 560)}          # a point inside each alignment track
PROBE = (560, 291, 40, 246)         # a column of pixels through both alignment tracks
PROBE_ROWS = ((0, 16), (230, 246))  # its rows inside the tumour and inside the normal alignment track
NAME_PANEL = (5, 234, 150, 14)                               # where the track name "tumor.cram Coverage" is written
MENU_EXPORT, MENU_COLOR, SUB_READ_STRAND = (790, 765), (790, 140), (1077, 216)
MENU_PROBE, SUBMENU_PROBE = (706, 757, 300, 16), (1020, 208, 240, 16)
FILE_MENU, FILE_SAVE_SESSION, FILE_SAVE_PNG = (18, 30), (62, 193), (73, 252)
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:6.1f}s] {msg}", flush=True)


def _desktop_module():
    """The `desktop` command as a module: its actions are then called in this process, which saves one
    interpreter start per mouse/keyboard action (a few hundred of them follow)."""
    import importlib.machinery
    import importlib.util
    import shutil
    path = shutil.which("desktop")
    loader = importlib.machinery.SourceFileLoader("desktop_cli", path)
    mod = importlib.util.module_from_spec(importlib.util.spec_from_loader("desktop_cli", loader))
    loader.exec_module(mod)
    return mod


_CLI = None


def desktop(*args):
    """Run one `desktop` action; `launch` and `windows` go through the command itself."""
    global _CLI
    args = [str(a) for a in args]
    if args[0] in ("launch", "windows"):
        return subprocess.run(["desktop", *args], capture_output=True, text=True, timeout=120).stdout
    if _CLI is None:
        _CLI = _desktop_module()
    import contextlib
    import io
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        _CLI.main(args)
    return buf.getvalue()


def click(xy, button="left", clicks=1):
    desktop("click", xy[0], xy[1], "--button", button, "--clicks", clicks)


def enter_text(xy, text):
    """Replace the content of a text field: click into it, select everything with the keyboard, type."""
    click(xy)
    desktop("key", "End", "shift+Home")
    desktop("type", text)
    desktop("key", "Return")


def dismiss_messages():
    """Close IGV message boxes (e.g. after a mistyped locus) so that they cannot block the next step."""
    for _ in range(5):
        w = find_window("Message") or find_window("Error") or find_window("Warning")
        if w is None:
            return
        log(f"closing a '{w[4]}' window")
        click((w[0] + w[2] // 2, w[1] + w[3] - 46))
        time.sleep(0.5)


def windows():
    """[(x, y, w, h, title)] of the top-level windows."""
    out = []
    for line in desktop("windows").splitlines():
        f = line.split(None, 7)
        if len(f) >= 7 and f[0].startswith("0x"):
            out.append((int(f[2]), int(f[3]), int(f[4]), int(f[5]), f[7] if len(f) > 7 else ""))
    return out


def find_window(title):
    for w in windows():
        if w[4].startswith(title):
            return w
    return None


def wait_window(title, present=True, timeout=10.0):
    end = time.time() + timeout
    while time.time() < end:
        w = find_window(title)
        if (w is not None) == present:
            return w if present else True
        time.sleep(0.25)
    return None


def pixels(region):
    """RGB rows of a part of the screen (small regions only: pure-Python PNG decoding)."""
    path = WORK / "strip.png"
    desktop("screenshot", path, "--region", ",".join(map(str, region)))
    raw = path.read_bytes()
    off, idat, w, h = 8, b"", 0, 0
    while off < len(raw):
        ln, typ = struct.unpack(">I4s", raw[off:off + 8])
        if typ == b"IHDR":
            w, h, depth, ctype = struct.unpack(">IIBB", raw[off + 8:off + 18])
            if depth != 8 or ctype not in (2, 6):
                raise RuntimeError("unexpected screenshot format")
            bpp = 3 if ctype == 2 else 4
        elif typ == b"IDAT":
            idat += raw[off + 8:off + 8 + ln]
        off += 12 + ln
    data = zlib.decompress(idat)
    stride = w * bpp
    rows, prev = [], bytearray(stride)
    for y in range(h):
        ft = data[y * (stride + 1)]
        cur = bytearray(data[y * (stride + 1) + 1:(y + 1) * (stride + 1)])
        for i in range(stride):
            a = cur[i - bpp] if i >= bpp else 0
            b = prev[i]
            c = prev[i - bpp] if i >= bpp else 0
            if ft == 1:
                cur[i] = (cur[i] + a) & 255
            elif ft == 2:
                cur[i] = (cur[i] + b) & 255
            elif ft == 3:
                cur[i] = (cur[i] + (a + b) // 2) & 255
            elif ft == 4:
                p = a + b - c
                pa, pb, pc = abs(p - a), abs(p - b), abs(p - c)
                cur[i] = (cur[i] + (a if pa <= pb and pa <= pc else b if pb <= pc else c)) & 255
        rows.append(bytes(cur))
        prev = cur
    return rows, bpp


def fraction(region, test):
    rows, bpp = pixels(region)
    n = hit = 0
    for r in rows:
        for i in range(0, len(r), bpp):
            n += 1
            hit += test(r[i], r[i + 1], r[i + 2])
    return hit / n


def is_ink(r, g, b):            # anything drawn on the light panel background
    return min(r, g, b) < 235


def is_menu(r, g, b):           # background of a Swing pop-up menu
    return 228 <= r <= 245 and 228 <= g <= 245 and 228 <= b <= 245


def strips():
    rows = pixels(PROBE)[0]
    return [rows[a:b] for a, b in PROBE_ROWS]


def reads_drawn(st):
    def frac(rows):
        n = sum(len(r) // 3 for r in rows)
        return sum(1 for r in rows for i in range(0, len(r), 3) if min(r[i], r[i + 1], r[i + 2]) < 235) / n
    return all(frac(rows) > 0.25 for rows in st)


def goto(locus, timeout=120.0):
    """Type a locus into IGV's locus box and wait until both alignment tracks show reads and stand still."""
    dismiss_messages()
    before = strips()
    enter_text(LOCUS_BOX, locus)
    end = time.time() + timeout
    changed, last = False, None
    while time.time() < end:
        time.sleep(0.15)
        cur = strips()
        if not changed:
            changed = cur != before or time.time() > end - timeout + 6
            last = cur
            continue
        if reads_drawn(cur) and cur == last:
            return True
        last = cur
    return False


def popup(sample):
    """Open the pop-up menu of an alignment track; True once the menu is on screen."""
    for _ in range(4):
        dismiss_messages()
        click(TRACK[sample], button="right")
        time.sleep(0.7)
        if fraction(MENU_PROBE, is_menu) > 0.6:
            return True
        desktop("key", "Escape")
        time.sleep(1.0)
    return False


def fill_dialog(title, path, timeout=60.0):
    """IGV's file dialog: put the full path into the file-name field and confirm."""
    w = wait_window(title, True, 8.0)
    if w is None:
        return False
    path = Path(path)
    if path.exists():
        path.unlink()
    enter_text((w[0] + 198, w[1] + 292), str(path))
    if not wait_window(title, False, 10.0):
        desktop("key", "Escape")
        return False
    end, size = time.time() + timeout, -1
    while time.time() < end:            # written and no longer growing
        time.sleep(0.25)
        if path.exists() and path.stat().st_size > 0:
            if path.stat().st_size == size:
                return True
            size = path.stat().st_size
    return False


def export_alignments(sample, path, lo, hi):
    """Alignment track > pop-up menu > Export Alignments...; accept only an export that spans [lo, hi]."""
    for attempt in range(4):
        if popup(sample):
            click(MENU_EXPORT)
            if fill_dialog("Save Visible Data", path):
                pos = [int(line.split("\t", 4)[3]) for line in open(path) if not line.startswith("@")]
                if len(pos) > 200 and min(pos) <= lo + 150 and max(pos) >= hi - 150:
                    log(f"exported {len(pos)} alignments of {sample} ({lo}-{hi})")
                    return
                log(f"export of {sample} incomplete ({len(pos)} records), repeating")
        desktop("key", "Escape")
        time.sleep(3.0)
    raise RuntimeError(f"could not export the alignments of {sample} for {lo}-{hi}")


def colour_by_read_strand(sample):
    for _ in range(4):
        if popup(sample):
            click(MENU_COLOR)
            time.sleep(0.7)
            if fraction(SUBMENU_PROBE, is_menu) > 0.6:
                click(SUB_READ_STRAND)
                time.sleep(0.8)
                return
        desktop("key", "Escape")
        time.sleep(1.0)
    raise RuntimeError("could not open the colour menu")


def file_menu(item, title, path):
    for attempt in range(4):
        if attempt:
            dismiss_messages()
        click(FILE_MENU)
        time.sleep(0.4)
        click(item)
        if fill_dialog(title, path):
            return True
        desktop("key", "Escape")
        time.sleep(1.0)
    return False


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    (OUT / "snapshots").mkdir(parents=True, exist_ok=True)
    contig, ref = sr.load_reference(INP / "reference" / "ref.fa")
    cands = sr.load_candidates(INP / "candidates.vcf")

    # 1. IGV
    desktop("launch", "/task/software/open_igv.sh")
    if wait_window("IGV", True, 240.0) is None:
        raise RuntimeError("IGV did not start")
    end = time.time() + 240
    while time.time() < end and fraction(NAME_PANEL, lambda r, g, b: max(r, g, b) < 110) < 0.01:
        time.sleep(1.0)                 # until the track names are written: tracks are loaded
    time.sleep(2.0)
    log("IGV is up")

    # 2. read-level evidence: export what the alignment tracks display, 20 kb at a time
    sams = {"tumor": [], "normal": []}
    step = 20000
    for lo in range(1, len(ref) + 1, step):
        a, b = max(1, lo - 200), min(len(ref), lo + step - 1)
        if not goto(f"{contig}:{a}-{b}"):
            raise RuntimeError(f"alignments did not appear for {a}-{b}")
        for sample in ("tumor", "normal"):
            path = WORK / f"{sample}_{a}_{b}.sam"
            export_alignments(sample, path, a, b)
            sams[sample].append(path)

    # 3. the review
    tumor = sr.Sample(sr.load_sam(sams["tumor"]), ref)
    normal = sr.Sample(sr.load_sam(sams["normal"]), ref)
    rows, res = sr.review(cands, ref, tumor, normal)
    (OUT / "review.tsv").write_text("\n".join(rows) + "\n")
    log(f"reviewed {len(cands)} candidates: {sum(1 for x in res if x['tag'] == 'PASS')} somatic")

    # 4. session: both tracks coloured by read strand, view on the first candidate
    def window(pos0):
        return f"{contig}:{pos0 + 1 - 40}-{pos0 + 1 + 40}"
    session = OUT / "review_session.xml"
    for attempt in range(3):
        goto(window(cands[0][1]))
        colour_by_read_strand("tumor")
        colour_by_read_strand("normal")
        goto(window(cands[0][1]))
        ok = file_menu(FILE_SAVE_SESSION, "Save Session", session)
        txt = session.read_text() if ok and session.exists() else ""
        if txt.count('colorOption="READ_STRAND"') == 2:
            break
        log("session not as required, repeating")
    else:
        raise RuntimeError("could not save the session")
    log("session saved")

    # 5. one image per candidate
    def snap(vid, pos0):
        path = OUT / "snapshots" / f"{vid}.png"
        for _ in range(3):
            if goto(window(pos0), timeout=60.0) and file_menu(FILE_SAVE_PNG, "Save image", path) \
                    and path.stat().st_size > 8000:
                return path
        raise RuntimeError(f"no image for {vid}")
    digest = {}
    for vid, pos0, _, _ in cands:
        digest[vid] = hashlib.sha256(snap(vid, pos0).read_bytes()).hexdigest()
    for vid, pos0, _, _ in cands:           # an image filed twice means a navigation was missed: redo it
        if sum(1 for d in digest.values() if d == digest[vid]) > 1:
            log(f"image of {vid} is not unique, repeating")
            digest[vid] = hashlib.sha256(snap(vid, pos0).read_bytes()).hexdigest()
    if len(set(digest.values())) != len(digest):
        raise RuntimeError("images are not unique")
    m = re.search(r'locus="([^"]+)"', session.read_text())
    log(f"{len(digest)} images saved; session locus {m.group(1) if m else '?'}")


if __name__ == "__main__":
    main()
