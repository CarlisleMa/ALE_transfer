"""Queue runner for the extra flash models (through OpenRouter, one local shim per model).

One Harbor job per (model, task) so that the order and the number of concurrent trials are under control:
  - light tasks first, the two compute-heavy tasks (hcs, md) last;
  - limits are read from flash_caps.json every round (total trials of this queue, of which medium-or-heavy,
    of which heavy), one set while the kimi-k3 re-runs are still going and one for afterwards. Measured on
    this 4-core host: a light trial uses almost no CPU (it waits for the model), so only the compute-bound
    ones have to be rationed;
  - restart-safe: running trials are found by their process, finished ones by their job directory;
  - a model is dropped from the queue when its billed cost passes BUDGET or its provider keeps refusing.
State is printed to stdout (-> flash_queue.log). Touch flash_queue.stop to stop launching new trials.
"""
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

S = Path(__file__).resolve().parent
LOGS = S / "logs"
REPO = Path("/home/user/ale-lifesci-gym")
BUDGET = 12.0          # USD per model, from the cost OpenRouter reports per call
BASE_MODELS = ("deepseek-v4.1-flash", "qwen3.8-flash", "glm-5.3-flash", "minimax-m3.1-flash")   # run all ten TASKS
# further (model, task) pairs are read from flash_extra.json every round and go to the FRONT of the queue
MODELS = {             # key -> (port, OpenRouter id, model_info for LiteLLM)
    "deepseek-v4.1-flash": (8771, "deepseek/deepseek-v4.1-flash",
                            {"max_input_tokens": 1048576, "max_output_tokens": 943718, "input_cost_per_token": 3.0e-07,
                             "output_cost_per_token": 1.2e-06, "cache_read_input_token_cost": 6e-09}),
    "qwen3.8-flash": (8772, "qwen/qwen3.8-flash",
                      {"max_input_tokens": 1000000, "max_output_tokens": 131072, "input_cost_per_token": 1.5e-07,
                       "output_cost_per_token": 4.7e-07, "cache_read_input_token_cost": 1.6e-08}),
    "glm-5.3-flash": (8773, "z-ai/glm-5.3-flash",
                      {"max_input_tokens": 1048576, "max_output_tokens": 131072, "input_cost_per_token": 1.5e-07,
                       "output_cost_per_token": 5.0e-07, "cache_read_input_token_cost": 3e-08}),
    "minimax-m3.1-flash": (8774, "minimax/minimax-m3.1-flash",
                           {"max_input_tokens": 1048576, "max_output_tokens": 524288, "input_cost_per_token": 1.5e-07,
                            "output_cost_per_token": 6.0e-07, "cache_read_input_token_cost": 1.5e-08}),
}
MODELS.update({
    # Gemini ran on the user's direct Google key until its prepaid credits were used up (17:39 UTC, HTTP 402);
    # from then on through OpenRouter, provider pinned to Google AI Studio (same model, same list price).
    "gemini-3.8-flash": (8778, "google/gemini-3.8-flash",
                         {"max_input_tokens": 1048576, "max_output_tokens": 65536, "input_cost_per_token": 7.5e-07,
                          "output_cost_per_token": 3.75e-06, "cache_read_input_token_cost": 7.5e-08}),
    "kimi-k3": (8776, "moonshotai/kimi-k3",             # OpenRouter, provider pinned to Moonshot AI; own shim and log
                {"max_input_tokens": 1048576, "max_output_tokens": 1048576, "input_cost_per_token": 3e-06,
                 "output_cost_per_token": 1.5e-05, "cache_read_input_token_cost": 3e-07}),
    # ceiling measurement (2026-11-14): a non-flash model on the pack, OpenRouter, provider pinned to Anthropic
    "claude-opus-5.5": (8779, "anthropic/claude-opus-5.5",
                        {"max_input_tokens": 1000000, "max_output_tokens": 64000, "input_cost_per_token": 4e-06,
                         "output_cost_per_token": 2e-05, "cache_read_input_token_cost": 2e-07}),
})
TASKS = ["crispr-screen-s0001-l4", "flow-compensation-s0001-l4", "trio-liftover-s0001-l4", "xlms-model-selection-s0001-l4",
         "variant-review-s0001-l4", "rnaseq-confound-s0001-l4", "amplicon-16s-s0001-l3", "somatic-tn-s0001-l3",
         "hcs-dose-response-s0001-l4", "md-system-repair-s0001-l3"]
MODEL_BUDGET = {"gemini-3.8-flash": 70.0, "kimi-k3": 40.0, "claude-opus-5.5": 60.0}   # kimi: the pack runs the user asked for (2026-11-13)   # USD billed through OpenRouter (default: BUDGET)
MODEL_CAP = {"gemini-3.8-flash": 4, "claude-opus-5.5": 2}     # trials of one model at a time (provider rate limits on tokens per minute)
HEAVY = ("hcs-", "md-")
MEDIUM = ("xlms-", "rnaseq-", "amplicon-", "somatic-")
CAPS = {"kimi": {"total": 6, "busy": 1, "heavy": 1}, "after": {"total": 7, "busy": 4, "heavy": 2}}


def caps(k_alive: bool) -> dict:
    try:
        return json.loads((S / "flash_caps.json").read_text())["kimi" if k_alive else "after"]
    except Exception:
        return CAPS["kimi" if k_alive else "after"]


def live() -> set:
    """(model, task) of every modelrun process of this queue that is alive, whoever started it."""
    out = subprocess.run(["pgrep", "-af", "[a]legen.modelrun"], capture_output=True, text=True).stdout
    return {(m.group(1), re.sub(r"__t(\d+)$", r"@\1", re.sub(r"__r(\d+)(?=$|__t)", r"#\1", m.group(2))))
            for m in re.finditer(r"alegen\.modelrun (\S+?)__(\S+) ", out) if m.group(1) in MODELS}


def say(*a):
    print(time.strftime("%H:%M:%S", time.gmtime()), *a, flush=True)


def kimi_alive() -> bool:
    return subprocess.run(["pgrep", "-f", "[k]imi_or_driver.sh"], capture_output=True).returncode == 0


def billed(port: int) -> float:
    f = LOGS / f"{port}.shim.log"
    return sum(float(x) for x in re.findall(r" cost=([0-9.e-]+)", f.read_text())) if f.exists() else 0.0


def refusing(port: int) -> str | None:
    f = LOGS / f"{port}.shim.log"
    if not f.exists():
        return None
    tail = f.read_text().splitlines()[-8:]
    if len(tail) == 8 and all("failed" in x for x in tail):
        return tail[-1][:200]
    if any("(not retried)" in x and ("HTTP 401" in x or "HTTP 402" in x or "HTTP 403" in x or "HTTP 429" in x) for x in tail[-2:]):
        return tail[-1][:200]
    return None


def jobname(key: str, task: str) -> str:
    """Queue entries may be 'task#n' (n-th attempt of the same task): job '<model>__<task>__r<n>'."""
    return f"{key}__{task.replace('#', '__r').replace('@', '__t')}"      # 'task#n@turns': attempt n, own turn limit


def launch(key: str, task: str) -> subprocess.Popen:
    port, mid, info = MODELS[key]
    job = jobname(key, task)
    turns = task.split("@")[1] if "@" in task else ("450" if "project" in task else "150")
    task = re.split(r"[#@]", task)[0]
    cmd = [str(REPO / ".venv/bin/python"), "-m", "alegen.modelrun", job, f"openai/{mid}", f"../tasks/{task}", "--n", "1",
           "--ak", f"api_base=http://127.0.0.1:{port}/v1", "--ak", "record_terminal_session=false",
           "--ak", "max_turns=" + turns,   # default 150; a project is several tasks in one session (450)
           "--ak", "model_info=" + json.dumps(info, separators=(",", ":"))]
    env = dict(os.environ, OPENAI_API_KEY="proxy-injected", ALEGEN_PRUNE_BUILD_CACHE="1", PATH=os.path.expanduser("~/.local/bin") + ":" + os.environ["PATH"])
    out = open(LOGS / f"flash_{job}.log", "w")
    return subprocess.Popen(cmd, cwd=REPO / "pipeline", env=env, stdout=out, stderr=subprocess.STDOUT, start_new_session=True)


def result(key: str, task: str) -> str:
    for r in (REPO / "results" / "jobs" / jobname(key, task)).glob("*/result.json"):
        j = json.loads(r.read_text())
        rw = (j.get("verifier_result") or {}).get("rewards")
        ex = (j.get("exception_info") or {}).get("exception_type")
        return f"reward={rw.get('reward') if isinstance(rw, dict) else rw} exc={ex}"
    return "no result.json"


def main() -> int:
    jobs = REPO / "results" / "jobs"
    running = live()
    # a job directory without a finished trial and without a live process was cut off (machine restart):
    # move it aside so that the trial is run again
    aside = jobs / "_interrupted"
    for jd in sorted(jobs.glob("*__*")):
        m = re.match(r"(.+?)__(.+)$", jd.name)
        if not m or (m.group(1), re.sub(r"__t(\d+)$", r"@\1", re.sub(r"__r(\d+)(?=$|__t)", r"#\1", m.group(2)))) in running or not jd.is_dir():
            continue
        if not list(jd.glob("*/result.json")):
            aside.mkdir(exist_ok=True)
            jd.rename(aside / f"{jd.name}.{int(time.time())}")
            say("interrupted trial moved aside, will be run again:", jd.name)
    queue = [(k, t) for t in TASKS for k in BASE_MODELS if (k, t) not in running and not (jobs / jobname(k, t)).exists()]
    seen_extra: set = set()
    try:      # pairs added later (flash_extra*.json) that have not run yet: the loop below would end before reading them
        late = [tuple(x) for f in ("flash_extra.json", "flash_extra_back.json") for x in json.loads((S / f).read_text())]
    except Exception:
        late = []
    late = [kt for kt in dict.fromkeys(late) if kt[0] in MODELS and kt not in running and kt not in queue
            and not (jobs / jobname(*kt)).exists()]
    queue = late + queue
    seen_extra.update(late)
    say(f"queue (re)started: {len(running)} trials already running, {len(queue)} queued")
    dropped: dict[str, str] = {}
    procs: list[subprocess.Popen] = []
    while queue or running:
        for p in procs:
            p.poll()                 # reap finished children
        now = live()
        for kt in sorted(running - now):
            say("finished", kt[0], kt[1], result(*kt), f"| billed {kt[0]}: ${billed(MODELS[kt[0]][0]):.2f}")
            shutil.rmtree(jobs / "_prepared" / jobname(*kt), ignore_errors=True)
        running = now
        try:                                   # late additions
            extra = [tuple(x) for x in json.loads((S / "flash_extra.json").read_text())]
        except Exception:
            extra = []
        new = [kt for kt in extra if kt not in seen_extra and kt[0] in MODELS and kt not in running
               and kt not in queue and not (jobs / jobname(*kt)).exists()]
        seen_extra.update(extra)
        if new:
            queue = new + queue
            say("added to the front of the queue:", ", ".join(f"{k}/{t}" for k, t in new))
        try:                                   # low-priority additions go to the back
            back = [tuple(x) for x in json.loads((S / "flash_extra_back.json").read_text())]
        except Exception:
            back = []
        newb = [kt for kt in back if kt not in seen_extra and kt[0] in MODELS and kt not in running
                and kt not in queue and not (jobs / jobname(*kt)).exists()]
        seen_extra.update(back)
        if newb:
            queue = queue + newb
            say("added to the back of the queue:", len(newb), "trials")
        for k, (port, _, _) in MODELS.items():
            if k in dropped:
                continue
            lim = MODEL_BUDGET.get(k, BUDGET)
            why = f"budget ${lim} passed (${billed(port):.2f})" if billed(port) > lim else refusing(port)
            if why:
                dropped[k] = why
                say("DROPPED", k, "-", why)
                queue = [(a, b) for a, b in queue if a != k]
                subprocess.run(["pkill", "-INT", "-f", f"[h]arbor run .*--job-name {k}__"])
        if (S / "flash_queue.stop").exists() and queue:
            say("stop file present: not launching the remaining", len(queue), "trials")
            queue = []
        k_alive = kimi_alive()
        c = caps(k_alive)
        free_gb = shutil.disk_usage("/home/user").free / 1e9
        mem_gb = int(re.search(r"MemAvailable:\s+(\d+)", Path("/proc/meminfo").read_text()).group(1)) / 1e6
        if queue and len(running) < c["total"] and free_gb > 2.0 and mem_gb > 3.0:
            n_heavy = sum(1 for _, t in running if t.startswith(HEAVY))
            n_busy = sum(1 for _, t in running if t.startswith(HEAVY + MEDIUM))
            for i, (k, t) in enumerate(queue):
                if t.startswith(HEAVY) and n_heavy >= c["heavy"]:
                    continue
                if t.startswith(HEAVY + MEDIUM) and n_busy >= c["busy"]:
                    continue
                if sum(1 for rk, _ in running if rk == k) >= MODEL_CAP.get(k, 99):
                    continue
                procs.append(launch(k, t))
                running.add((k, t))
                del queue[i]
                say("started", k, t, f"| running {len(running)}/{c['total']} (busy {n_busy}, heavy {n_heavy}) | queued {len(queue)} "
                    f"| kimi re-runs alive: {k_alive} | disk {free_gb:.1f} GB, mem {mem_gb:.1f} GB free")
                time.sleep(20)          # let the new process show up before the next round counts it
                break
        time.sleep(10)
    say("queue done; dropped:", dropped or "none")
    for k, (port, _, _) in MODELS.items():
        say("billed", k, f"${billed(port):.2f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
