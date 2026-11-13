# Scripts used for the model runs on the build host (kept for the record, not part of the pipeline)

These are the scripts that drove the Harbor runs of 2026-11-12 on a 4-core cloud container whose task containers
cannot reach the internet. They contain absolute paths of that session and are not meant to run elsewhere unchanged;
on a normal machine `harbor run -p tasks -a <agent> -m <model>` is all that is needed (README section 12).

- `flash_queue.py` - queue runner: one Harbor job per (model, task) through `python -m alegen.modelrun`, one local
  streaming shim per model, limits on concurrent compute-bound trials, a spending cap per model from the cost
  OpenRouter reports, repeated attempts (`task#n`, job `<model>__<task>__r<n>`), an own turn limit (`task#n@turns`),
  and recovery of trials cut off by a machine restart.
- `recover.sh` - restarts Docker, the shims and the queue after a restart of the machine.
- `status.py`, `wait_re.py`, `wait_file.py` - read-only status and wait helpers.
- `billing.py` - sums the cost OpenRouter reported per call from the shim logs into `results/billing.json`
  (the source of the cost line in the main README).
- `kimi_or_driver.sh`, `kimi_or_guard.sh` - the four kimi-k3 tasks that were run again through OpenRouter, with a
  watchdog that stops on a spending cap or on repeated provider failures.

Lessons that went into them: a provider that answers "insufficient balance" must not be retried, a rate limit must;
on this host every image build leaves build cache of the size of the image (`docker builder prune -af`); a session
without a foreground or background task gets its machine restarted.
