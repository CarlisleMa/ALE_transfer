#!/usr/bin/env python3
"""Command-line front end to the desktop actions of the MCP server (server.py), for agents and
harnesses that only have a shell. Same display, same actions, same xdotool/scrot/wmctrl calls.

    desktop screenshot FILE.png [--region X,Y,W,H]   capture the desktop (1280x800), or a part of it, to a PNG file
    desktop click X Y [--button left|middle|right] [--clicks N]
    desktop move X Y
    desktop drag X1 Y1 X2 Y2                   press left at (X1,Y1), move, release at (X2,Y2)
    desktop type TEXT                          type text into the focused window
    desktop key COMBO [COMBO ...]              xdotool syntax: Return, Escape, ctrl+a, alt+F4, Down ...
    desktop scroll X Y CLICKS                  wheel at (X,Y); positive = down, negative = up
    desktop launch COMMAND [ARG ...]           start a GUI program without waiting for it
    desktop windows                            list top-level windows (id, desktop, x, y, w, h, host, title)
    desktop wait SECONDS

The virtual display (Xvfb :99 with a window manager) is started on first use.
"""
import argparse
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.realpath(__file__))
os.environ["DISPLAY"] = os.environ.get("DISPLAY", ":99")


def ensure(force=False):
    """Bring the display up if its X socket is not there yet (cheap check, no process per action)."""
    sock = "/tmp/.X11-unix/X" + os.environ["DISPLAY"].split(":")[-1].split(".")[0]
    if force or not os.path.exists(sock):
        subprocess.run(["bash", os.path.join(HERE, "start_desktop.sh")], capture_output=True, timeout=60)


def xdo(*args):
    ensure()
    r = subprocess.run(["xdotool", *map(str, args)], capture_output=True, text=True, timeout=60)
    if r.returncode != 0:           # stale socket of a display that is gone: start it and try once more
        ensure(force=True)
        r = subprocess.run(["xdotool", *map(str, args)], capture_output=True, text=True, timeout=60)
    out = (r.stdout + r.stderr).strip()
    if out:
        print(out)
    return r.returncode


def main(argv=None):
    ap = argparse.ArgumentParser(prog="desktop", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("screenshot")
    p.add_argument("file")
    p.add_argument("--region", default=None, help="X,Y,W,H in pixels (default: the whole desktop)")
    p = sub.add_parser("click")
    p.add_argument("x", type=int)
    p.add_argument("y", type=int)
    p.add_argument("--button", default="left", choices=["left", "middle", "right"])
    p.add_argument("--clicks", type=int, default=1)
    p = sub.add_parser("move")
    p.add_argument("x", type=int)
    p.add_argument("y", type=int)
    p = sub.add_parser("drag")
    for n in ("x1", "y1", "x2", "y2"):
        p.add_argument(n, type=int)
    p = sub.add_parser("type")
    p.add_argument("text", nargs="+")
    p = sub.add_parser("key")
    p.add_argument("combo", nargs="+")
    p = sub.add_parser("scroll")
    p.add_argument("x", type=int)
    p.add_argument("y", type=int)
    p.add_argument("clicks", type=int)
    p = sub.add_parser("launch")
    p.add_argument("command", nargs=argparse.REMAINDER)
    sub.add_parser("windows")
    p = sub.add_parser("wait")
    p.add_argument("seconds", type=float)
    a = ap.parse_args(argv)

    if a.cmd == "screenshot":
        ensure()
        out = os.path.abspath(a.file)
        os.makedirs(os.path.dirname(out), exist_ok=True)
        tmp = "/tmp/desktop_shot_%d.png" % os.getpid()
        cmd = ["scrot", "-o", "-z"]
        if a.region:
            try:
                x, y, w, h = (int(v) for v in a.region.split(","))
            except ValueError:
                ap.error("--region takes X,Y,W,H")
            cmd += ["-a", f"{x},{y},{w},{h}"]
        r = subprocess.run(cmd + [tmp], capture_output=True, text=True, timeout=60)
        if r.returncode != 0 or not os.path.exists(tmp):
            print((r.stdout + r.stderr).strip() or "screenshot failed", file=sys.stderr)
            return 1
        with open(tmp, "rb") as src, open(out, "wb") as dst:
            dst.write(src.read())
        os.unlink(tmp)
        print(out)
        return 0
    if a.cmd == "click":
        b = {"left": 1, "middle": 2, "right": 3}[a.button]
        return xdo("mousemove", a.x, a.y, "click", "--repeat", max(1, a.clicks), "--delay", 80, b)
    if a.cmd == "move":
        return xdo("mousemove", a.x, a.y)
    if a.cmd == "drag":
        return xdo("mousemove", a.x1, a.y1, "mousedown", 1, "mousemove", a.x2, a.y2, "mouseup", 1)
    if a.cmd == "type":
        return xdo("type", "--delay", 15, "--", " ".join(a.text))
    if a.cmd == "key":
        return xdo("key", "--", *a.combo)
    if a.cmd == "scroll":
        return xdo("mousemove", a.x, a.y, "click", "--repeat", abs(a.clicks) or 1, 5 if a.clicks > 0 else 4)
    if a.cmd == "launch":
        if not a.command:
            ap.error("launch needs a command")
        ensure()
        subprocess.Popen(a.command, stdout=open("/tmp/launch.log", "ab"), stderr=subprocess.STDOUT,
                         stdin=subprocess.DEVNULL, start_new_session=True)
        print("started")
        return 0
    if a.cmd == "windows":
        ensure()
        r = subprocess.run(["wmctrl", "-l", "-G"], capture_output=True, text=True, timeout=30)
        print(r.stdout.strip() or "(no windows)")
        return 0
    if a.cmd == "wait":
        time.sleep(max(0.0, min(300.0, a.seconds)))
        return 0
    return 2


if __name__ == "__main__":
    sys.exit(main())
