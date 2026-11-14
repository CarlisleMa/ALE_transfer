#!/bin/bash
# Bring the benchmark machinery back after a machine restart (idempotent): Docker, shims, queue.
S=/tmp/claude-0/-home-user/5a931900-2990-511a-8e45-6939bdda2f59/scratchpad
L=$S/logs
PY=/root/.local/share/uv/tools/harbor/bin/python
SH=/home/user/ale-lifesci-gym/pipeline/alegen/stream_shim.py
OR=https://openrouter.ai/api/v1
if ! docker info >/dev/null 2>&1; then
  # after a machine restart the pid file of the old daemon is stale; if its number is taken by another process, dockerd refuses to start
  pgrep -x dockerd >/dev/null || rm -f /var/run/docker.pid
  (nohup dockerd >> $L/dockerd.log 2>&1 &)
  timeout 90 bash -c 'until docker info >/dev/null 2>&1; do sleep 2; done'
  if ! docker info >/dev/null 2>&1; then echo "DOCKER DID NOT START - queue not started; see $L/dockerd.log"; exit 1; fi
  docker container prune -f >/dev/null 2>&1; docker network prune -f >/dev/null 2>&1
  echo "docker restarted"
fi
shim() { pgrep -f "[s]tream_shim.py $1 " >/dev/null || (cd $L && setsid nohup $PY $SH "$@" > shim_$1.out 2>&1 < /dev/null &); }
shim 8771 $OR '{"provider":{"only":["fireworks","together"],"order":["fireworks","together"]}}'
shim 8772 $OR '{"provider":{"only":["alibaba"],"allow_fallbacks":false}}'
shim 8773 $OR '{"provider":{"only":["z-ai"],"allow_fallbacks":false}}'
shim 8774 $OR '{"provider":{"only":["minimax"],"allow_fallbacks":false}}'
shim 8778 $OR '{"provider":{"only":["google-ai-studio"],"allow_fallbacks":false}}'   # Gemini via OpenRouter (direct key out of credit)
shim 8776 $OR '{"provider":{"only":["moonshotai"],"allow_fallbacks":false}}'
shim 8779 $OR '{"provider":{"only":["anthropic"],"allow_fallbacks":false}}'   # ceiling model
pgrep -f "[f]lash_queue.py" >/dev/null || { (cd $S && setsid nohup python3 flash_queue.py >> $L/flash_queue.log 2>&1 < /dev/null &); echo "queue restarted"; }
sleep 3
echo "uptime: $(uptime | sed 's/,.*//') | shims: $(pgrep -f '[s]tream_shim.py' | wc -l) | queue: $(pgrep -f '[f]lash_queue.py' | wc -l) | disk free: $(df -h /home/user | awk 'NR==2{print $4}')"
