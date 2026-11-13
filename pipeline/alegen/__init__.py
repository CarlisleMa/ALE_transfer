"""alegen: generate ALE-style long-horizon life-science tasks as Harbor task directories.

Design rules every family follows (see README.md, "Pipeline"):
  1. truth-by-construction  - the hidden reference is planted by the generator, never looked up.
  2. goal explicit, obstacle discoverable, recipe withheld - traps live in the data, not the prompt.
  3. verifier is attacked before release - shortcut submissions must score low, two independent
     correct solutions must both score full.
"""
__version__ = "0.1.0"
