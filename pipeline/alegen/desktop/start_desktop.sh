#!/bin/bash
# Idempotent: bring up the virtual display and a window manager if they are not running yet.
export DISPLAY="${DISPLAY:-:99}"
if ! xdpyinfo -display "$DISPLAY" >/dev/null 2>&1; then
  nohup Xvfb "$DISPLAY" -screen 0 1280x800x24 -nolisten tcp >/tmp/xvfb.log 2>&1 &
  for _ in $(seq 1 100); do xdpyinfo -display "$DISPLAY" >/dev/null 2>&1 && break; sleep 0.1; done
  nohup openbox >/tmp/openbox.log 2>&1 &
  sleep 0.5
fi
