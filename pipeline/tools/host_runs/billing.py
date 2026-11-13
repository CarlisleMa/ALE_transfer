"""Sum what OpenRouter billed, from the cost it reports with every call (logged by alegen/stream_shim.py).

  python3 billing.py <directory with <port>.shim.log files> <out.json>

One relay process (port) per model and route; the mapping below is the one used on the build host on
2026-11-12. Only calls that went through OpenRouter carry a cost; the direct Gemini and Kimi routes of the
morning are not in these logs.
"""
import json
import re
import sys
import time
from pathlib import Path

PORTS = {8770: "kimi-k3", 8776: "kimi-k3", 8778: "gemini-3.8-flash", 8771: "deepseek-v4.1-flash",
         8772: "qwen3.8-flash", 8773: "glm-5.3-flash", 8774: "minimax-m3.1-flash"}


def main(logs: Path, out: Path) -> int:
    per: dict[str, dict] = {}
    for port, model in PORTS.items():
        f = logs / f"{port}.shim.log"
        if not f.exists():
            continue
        s = f.read_text(errors="replace")
        m = per.setdefault(model, {"calls_ok": 0, "calls_failed": 0, "billed_usd": 0.0, "providers": {}})
        m["calls_ok"] += len(re.findall(r"(?m)^\S+ ok ", s))
        m["calls_failed"] += len(re.findall(r"(?m)^\S+ attempt \d+ failed", s))
        m["billed_usd"] += sum(float(x) for x in re.findall(r" cost=([0-9.e-]+)", s))
        for p in re.findall(r" provider='([^']*)'", s):
            m["providers"][p] = m["providers"].get(p, 0) + 1
    for m in per.values():
        m["billed_usd"] = round(m["billed_usd"], 2)
    doc = {"as_of_utc": time.strftime("%Y-%m-%d %H:%M", time.gmtime()), "source": "cost field of every OpenRouter response",
           "models": per, "total_usd": round(sum(m["billed_usd"] for m in per.values()), 2)}
    out.write_text(json.dumps(doc, indent=1) + "\n")
    print(json.dumps(doc, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main(Path(sys.argv[1]), Path(sys.argv[2])))
