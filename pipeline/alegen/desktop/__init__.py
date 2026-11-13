"""Desktop layer: adds an on-demand X display, a window manager and an MCP server with desktop
actions to an agent image, so a task can offer GUI professional software (IGV, Fiji, PyMOL, ...)
to agents that can use it while staying solvable from the terminal."""
from pathlib import Path

HERE = Path(__file__).parent

APT = ("RUN apt-get update && apt-get install -y --no-install-recommends "
       "xvfb openbox xdotool scrot wmctrl x11-utils fonts-dejavu-core "
       "libxtst6 libxrender1 libxi6 libfontconfig1 libfreetype6 "
       "&& rm -rf /var/lib/apt/lists/*\n"
       "ENV DISPLAY=:99\n")
PIP = ["mcp==1.30.0"]

FILES = {"desktop/server.py": HERE / "server.py", "desktop/start_desktop.sh": HERE / "start_desktop.sh"}

MCP_SERVER = {"name": "desktop", "transport": "stdio", "command": "/opt/tools/bin/python",
              "args": ["/task/software/desktop/server.py"]}

# ---- level 4 additions (the files above are shipped in level-3 tasks and stay byte-identical) ----
# `desktop` command: the same actions as the MCP server for shell-only agents and harnesses.
CLI_FILES = {**FILES, "desktop/desktop": HERE / "cli.py"}
CLI_LINK = "RUN ln -s /task/software/desktop/desktop /usr/local/bin/desktop\n"
