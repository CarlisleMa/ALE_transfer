#!/bin/bash
# Run the four kimi-k3 tasks that have no valid score, through OpenRouter (provider pinned to Moonshot AI), two at a time.
cd /home/user/ale-lifesci-gym/pipeline || exit 1
export PATH="$HOME/.local/bin:$PATH" OPENAI_API_KEY=proxy-injected
MI='model_info={"max_input_tokens":1048576,"max_output_tokens":1048576,"input_cost_per_token":3e-06,"output_cost_per_token":1.5e-05,"cache_read_input_token_cost":3e-07}'
run() {  # job task task
  [ -e "$STOP" ] && return
  ../.venv/bin/python -m alegen.modelrun "$1" openai/moonshotai/kimi-k3 "../tasks/$2" "../tasks/$3" --n 2 \
     --ak api_base=http://127.0.0.1:8770/v1 --ak record_terminal_session=false --ak max_turns=150 --ak "$MI"
}
STOP=$1
date -u +"%H:%M:%S start or1"; run kimi-k3-or1 hcs-dose-response-s0001-l4 variant-review-s0001-l4
date -u +"%H:%M:%S start or2"; run kimi-k3-or2 md-system-repair-s0001-l3 rnaseq-confound-s0001-l4
date -u +"%H:%M:%S driver done"
