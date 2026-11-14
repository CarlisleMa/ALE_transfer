#!/bin/bash
# Verifier entrypoint (runs in the dedicated verifier image; the agent never sees /tests).
mkdir -p /logs/verifier
python3 /tests/score.py --output /task/output --reference /tests/reference --logs /logs/verifier
