"""Read-only status of the benchmark runs: machine, queue, per-trial results, billed cost."""
import glob
import json
import os
import re
import shutil
import subprocess
import time
from pathlib import Path

S = Path(__file__).resolve().parent
LOGS = S / "logs"
JOBS = Path("/home/user/ale-lifesci-gym/results/jobs")
PORTS = {8770: "kimi-k3 (OpenRouter)", 8776: "kimi-k3 (level 5/6)", 8778: "gemini-3.8-flash (OpenRouter)", 8771: "deepseek-v4.1-flash", 8772: "qwen3.8-flash", 8773: "glm-5.3-flash", 8774: "minimax-m3.1-flash"}


def sh(*cmd):
    return subprocess.run(list(cmd), capture_output=True, text=True).stdout


print("time", time.strftime("%H:%M:%S", time.gmtime()), "| load", open("/proc/loadavg").read().split()[:3],
      "| disk free %.1f GB" % (shutil.disk_usage("/home/user").free / 1e9),
      "| mem avail %.1f GB" % (int(re.search(r"MemAvailable:\s+(\d+)", open("/proc/meminfo").read()).group(1)) / 1e6))
for name, pat in (("queue", "[f]lash_queue.py"), ("kimi driver", "[k]imi_or_driver.sh"), ("kimi guard", "[k]imi_or_guard.sh")):
    print(f"{name}: {'alive' if sh('pgrep', '-f', pat).strip() else 'not running'}", end=" | ")
print()
stop = S / "kimi_or.stop"
if stop.exists():
    print("KIMI STOP FILE:", stop.read_text().strip())
kl = LOGS / "kimi_or.log"
if kl.exists():
    print("kimi driver:", " / ".join(x for x in kl.read_text().splitlines() if re.search(r"start or|driver done", x)))
q = LOGS / "flash_queue.log"
if q.exists():
    lines = q.read_text().splitlines()
    for x in lines:
        if "DROPPED" in x or "queue done" in x or "stop file" in x:
            print("QUEUE:", x[:220])
    print("queue last line:", lines[-1][:210] if lines else "")

print("\n--- trials (jobs named kimi-k3-or* or <model>__<task>)")
rows = []
for T in sorted(glob.glob(str(JOBS / "*" / "*" / ""))):
    job = T.rstrip("/").split("/")[-2]
    if not os.path.isdir(T) or job == "_prepared" or not (job.startswith("kimi-k3-or") or "__" in job):
        continue
    trial = T.rstrip("/").split("/")[-1]
    task = trial.split("__")[0]
    model = job.split("__")[0] if "__" in job else "kimi-k3"
    r = os.path.join(T, "result.json")
    tj = os.path.join(T, "agent", "trajectory.json")
    turns = None
    if os.path.exists(tj):
        try:
            turns = sum(1 for s in json.load(open(tj))["steps"] if s.get("source") == "agent")
        except Exception:
            pass
    if os.path.exists(r):
        d = json.load(open(r))
        rw = (d.get("verifier_result") or {}).get("rewards")
        rw = rw.get("reward") if isinstance(rw, dict) else rw
        ex = (d.get("exception_info") or {}).get("exception_type")
        low = ""
        for f in glob.glob(os.path.join(T, "verifier", "**", "details.json"), recursive=True):
            dj = json.load(open(f))
            comps = dj.get("components") or []
            if comps:
                c = min(comps, key=lambda c: c["value"])
                low = "" if c["value"] == 1.0 else f"{c['name']}={c['value']:.2f}"
            elif dj.get("failure"):
                low = dj["failure"][:60]
        from datetime import datetime
        f_ = lambda x: datetime.fromisoformat(str(x).replace("Z", "+00:00"))
        mins = (f_(d["finished_at"]) - f_(d["started_at"])).total_seconds() / 60 if d.get("finished_at") else None
        rows.append((model, task, "DONE", f"reward={rw}", f"{mins:.0f} min" if mins is not None else "", f"turns={turns}", ex or "", low))
    else:
        started = time.strftime("%H:%M", time.gmtime(os.path.getctime(T)))
        rows.append((model, task, "running", f"since {started}", "", f"turns={turns}", "", ""))
for row in sorted(rows):
    print("  %-20s %-30s %-8s %-22s %-8s %-10s %s %s" % row)

print("\n--- billed (from the cost OpenRouter reports per call)")
tot = 0.0
for port, name in PORTS.items():
    f = LOGS / f"{port}.shim.log"
    if f.exists():
        s = f.read_text()
        c = sum(float(x) for x in re.findall(r" cost=([0-9.e-]+)", s))
        tot += c
        print("  %-22s ok %-5d failed %-3d $%.2f" % (name, s.count(" ok "), s.count("failed"), c))
print("  total $%.2f" % tot)
