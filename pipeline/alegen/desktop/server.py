#!/usr/bin/env python3
"""GUI-as-tool: a small MCP (stdio) server that exposes the container's X desktop as ordinary
tools, the way ALE lifts CLI agents to generalist computer-use agents. The desktop is an Xvfb
display started on demand, so the same image serves terminal-only and GUI-capable agents.

Tools: screenshot, click, move, drag, type_text, key, scroll, launch, list_windows, wait.
"""
import os
import shlex
import subprocess
import time

from mcp.server.fastmcp import FastMCP, Image

DISPLAY = os.environ.get("DISPLAY", ":99")
os.environ["DISPLAY"] = DISPLAY
START = os.path.join(os.path.dirname(os.path.abspath(__file__)), "start_desktop.sh")
mcp = FastMCP("desktop")


def _ensure():
    subprocess.run(["bash", START], capture_output=True, timeout=60)


def _xdo(*args):
    _ensure()
    r = subprocess.run(["xdotool", *map(str, args)], capture_output=True, text=True, timeout=30)
    return (r.stdout + r.stderr).strip() or "ok"


@mcp.tool()
def screenshot() -> Image:
    """Capture the whole desktop (1280x800) as a PNG."""
    _ensure()
    path = "/tmp/desktop_shot.png"
    subprocess.run(["scrot", "-o", "-z", path], capture_output=True, timeout=30)
    with open(path, "rb") as fh:
        return Image(data=fh.read(), format="png")


@mcp.tool()
def click(x: int, y: int, button: str = "left", clicks: int = 1) -> str:
    """Click at pixel (x, y). button: left | middle | right. clicks=2 for a double click."""
    b = {"left": 1, "middle": 2, "right": 3}.get(button, 1)
    return _xdo("mousemove", x, y, "click", "--repeat", max(1, clicks), "--delay", 80, b)


@mcp.tool()
def move(x: int, y: int) -> str:
    """Move the pointer to pixel (x, y)."""
    return _xdo("mousemove", x, y)


@mcp.tool()
def drag(x1: int, y1: int, x2: int, y2: int) -> str:
    """Press the left button at (x1, y1), move to (x2, y2), release."""
    return _xdo("mousemove", x1, y1, "mousedown", 1, "mousemove", x2, y2, "mouseup", 1)


@mcp.tool()
def type_text(text: str) -> str:
    """Type text into the focused window."""
    return _xdo("type", "--delay", 15, "--", text)


@mcp.tool()
def key(combo: str) -> str:
    """Press a key or chord in xdotool syntax, e.g. Return, Escape, ctrl+s, alt+F4."""
    return _xdo("key", "--", combo)


@mcp.tool()
def scroll(x: int, y: int, clicks: int) -> str:
    """Scroll the wheel at (x, y); positive = down, negative = up."""
    return _xdo("mousemove", x, y, "click", "--repeat", abs(clicks) or 1, 5 if clicks > 0 else 4)


@mcp.tool()
def launch(command: str) -> str:
    """Start a GUI program on the desktop without waiting for it (e.g. /task/software/open_igv.sh)."""
    _ensure()
    subprocess.Popen(shlex.split(command), stdout=open("/tmp/launch.log", "ab"), stderr=subprocess.STDOUT,
                     start_new_session=True)
    return "started"


@mcp.tool()
def list_windows() -> str:
    """List visible top-level windows (id, geometry, title)."""
    _ensure()
    r = subprocess.run(["wmctrl", "-l", "-G"], capture_output=True, text=True, timeout=30)
    return r.stdout.strip() or "(no windows)"


@mcp.tool()
def wait(seconds: float) -> str:
    """Sleep for up to 30 seconds (let a window finish drawing)."""
    time.sleep(max(0.0, min(30.0, seconds)))
    return "ok"


if __name__ == "__main__":
    _ensure()
    mcp.run()
