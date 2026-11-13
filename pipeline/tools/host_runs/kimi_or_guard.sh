#!/bin/bash
# Watchdog for the OpenRouter runs: stop everything if the billed cost passes the budget or the provider keeps refusing.
LOG=$1; STOP=$2; BUDGET=$3
while pgrep -f "[k]imi_or_driver.sh" >/dev/null; do
  cost=$(grep -o ' cost=[0-9.e-]*' "$LOG" | cut -d= -f2 | awk '{s+=$1} END{printf "%.2f", s}')
  fails=$(tail -8 "$LOG" | grep -c 'failed')
  if awk "BEGIN{exit !($cost > $BUDGET)}"; then echo "$(date -u +%H:%M:%S) budget exceeded: \$$cost" > "$STOP"; fi
  if [ "$fails" -ge 8 ]; then echo "$(date -u +%H:%M:%S) eight consecutive upstream failures: $(tail -1 "$LOG" | cut -c1-200)" > "$STOP"; fi
  if [ -e "$STOP" ]; then
    pkill -INT -f "[h]arbor run .*kimi-k3-or"; sleep 45; pkill -TERM -f "[h]arbor run .*kimi-k3-or"; pkill -f "[k]imi_or_driver.sh"
    exit 0
  fi
  sleep 60
done
