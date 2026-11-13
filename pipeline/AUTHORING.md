# How to author a task family for `alegen`

A *family* is a seeded generator that emits ALE-style long-horizon life-science tasks as Harbor
task directories. One family = one professional workflow; one `(seed, level)` = one task instance.
The reference implementation to imitate is `alegen/families/crispr_screen/`. Read it first
(`__init__.py`, `assets/score.py`, `assets/crispr_pipeline.py`, `assets/solve_alt.py`).

## The three non-negotiable design rules

1. **Truth by construction.** Decide the answer first, then generate the data from it (forward
   simulation on top of real reference material where that adds realism). The hidden reference is
   the *planted* quantity, never the output of the reference pipeline. If a tool-derived quantity
   is unavoidable, use the consensus of two independent tools and drop the disagreements from
   scoring (gray zone). Nothing the verifier needs may be retrievable from the public internet,
   from the task's inputs by a trivial transform, or from prior knowledge (blind real identifiers
   when they would leak the answer).
2. **Goal explicit, obstacle discoverable, recipe withheld.** The instruction states exactly what
   to deliver and defines every scored quantity unambiguously. It never names a trap, a threshold,
   a tolerance or the method. Every trap must be discoverable by a standard QC step an expert would
   run, or from domain jargon in the provided lab notes; write that evidence into `Trap.discoverable_by`.
   A requirement that cannot be discovered is a hidden requirement, which is a defect.
3. **Attack the verifier before release.** A task ships only if the shortcut submissions score low
   *and* two independent correct solutions score full (see "Battery").

## Files of a family

```
alegen/families/<family>/__init__.py     generate(seed:int, level:int, workdir:Path) -> TaskBuild
alegen/families/<family>/assets/score.py verifier (uses tests/_scoring.py, stdlib + numpy only unless justified)
alegen/families/<family>/assets/solve.sh oracle; `SOLUTION=2` must select an independent 2nd solution
alegen/families/<family>/assets/...      helper scripts
docs/families/<family>.md                design notes (see "Deliverables")
worklog.d/<family>.md                    dated log of what you tried, what failed, what you changed
```

Do **not** edit the core (`spec.py`, `harbor.py`, `validate.py`, `scoring.py`, `dockerlib.py`,
`cli.py`): other authors are working in parallel. If the core lacks something, write your own
Dockerfile string / helper inside your family and mention it in your report.

## Container contract (mirrors ALE's input/ software/ output/ reference/)

- Agent sees `/task/input` (read-only data + lab notes), `/task/software` (optional wrappers) and
  writes deliverables to `/task/output`. Tools live in one conda env at `/opt/tools` on `PATH`.
- Build the agent image with `dockerlib.agent_dockerfile(conda=[...pinned...], pip=[...])`.
  Pin versions (`tool=1.2.3`). Only install what a professional would have for this workflow.
- The verifier runs in its own image (`dockerlib.verifier_dockerfile(...)`) that bundles
  `tests/` (score.py + `reference/`). The agent image must not contain the reference or any file
  from which it can be derived. Only `/task/output` is transferred to the verifier.
- No network while the agent works (allowlist = model API hosts only) and none in the verifier.
  Everything the oracle needs must be inside the image.

## Budgets

| item | budget |
|---|---|
| `environment/input` total size | <= 25 MB (hard cap 40 MB) |
| generation time | <= 10 min |
| oracle runtime in the agent image | <= 20 min on 4 cores, <= 8 GB RAM |
| verifier runtime | <= 5 min |
| agent image | <= 3 GB |

## Levels

`level=1`: the long workflow without traps. `level=2`: some traps. `level=3`: all traps (the first
shipped difficulty). `level=4`: see "Level 4" below. Trap parameters and which sample/gene carries them must vary with the seed so
instances cannot be memorised. `task_id = f"{family-with-dashes}-s{seed:04d}-l{level}"`.

## Level 4 (second iteration, after the small-model pilot)

The pilot (README section 6) showed that a small model solved seven of the ten level-3 tasks in
2-15 minutes: every decision point that lab notes or a routine QC step reveal was handled. What it
failed on were (a) a final quantity that had to be *good enough* and needed self-checking and method
tuning, (b) a plausible shortcut that the generator had closed, (c) a process that had to be run to
a convergence criterion. Level 4 is built from those, not from more QC-discoverable traps.

`level=4` = the level-3 task plus **at least two** of these elements (deep families: three, one of
them H, C or T). Name them in `Trap.name` with the letter prefix (`Q_...`, `P_...`):

| code | element | what it must look like |
|---|---|---|
| Q | quality bar with an in-data control | A scored final quantity for which the default method gives a plausible but out-of-band value. The inputs contain a control an expert would use to measure their own error (mock sample, spike-in, reference compound, FMO tube, manually annotated fields, technical replicate). The SOP states the acceptance criterion on the control. Meeting it must be necessary and close to sufficient for the hidden scored set. |
| P | prior decoupled from truth | An authoritative-looking default (embedded metadata, a predecessor's parameter, a textbook rule of thumb, per-variant instead of per-haplotype logic) gives the wrong answer in this instance; the evidence for the right one is in the data. The generator must make rule and truth disagree on enough scored items to fail. |
| H | inherited project | The inputs include a predecessor's half-finished analysis: scripts, intermediate files, a draft result, a hand-over note. Some of it is right, some carries planted defects (seeded subset of a fixed checklist; include items that are *not* defective). The SOP contains the checklist. Deliverables: the corrected final result **and** `audit.tsv` (`item`, `verdict` = `ok`/`defective`), scored with MCC against the planted defects, plus where natural an erratum for the draft. Redoing everything from scratch is a legitimate strategy and must be able to score full; trusting the inherited intermediates must fail. Build defects as *mutations of the reference solution* (off-by-one, wrong sort order, wrong parameter, stale input, reversed contrast...) so the same operator can be reused across seeds and families. |
| C | dependent stages | At least three stages where a stage's output is the next stage's input and an early mistake only shows in the final deliverable. Final deliverable >= 0.65 of the weight. |
| T | workstation realism | A step for which the only software that can read the data is the professional tool (GUI or its own command line), as on a locked-down workstation. The instruction says so plainly. |
| I | run to convergence | A calibration/optimisation/equilibration that must be iterated until a stated criterion holds and needs real planning of the time budget (not idle waiting: <= 20 min of oracle compute). |

Fairness rules (level 4 is where unfair tasks are easiest to write):
- Every acceptance criterion is stated (instruction or SOP). Methods and thresholds *of the method*
  are not. No hidden requirement, no ambiguous definition: if two experts could read a definition
  differently, fix the definition.
- Two independent correct solutions still score >= 0.99 (`independent_solution` probe). A tolerance
  tighter than the gap between your own two solutions is a defect.
- Difficulty must not come from compute padding, obfuscated files, or contradicting the SOP.
- Mandatory extra attack `l3_reference`: the level-3 reference method applied to the level-4 inputs
  (it handles every level-3 decision point and none of the new elements), ceiling
  `pass_threshold - 0.15`. One "ignore exactly this element" attack per new element, same ceiling.
  The validator refuses a level-4 task without `l3_reference`.
- A level-4 reference solution should have `oracle_steps >= 20`; record it honestly.
- Levels 1-3 must keep generating byte-identical tasks (check `tasks/<family>-s0001-l3` with
  `diff -r` after your change, ignoring nothing).
- After release the task is piloted with a small model and every failure is read: a failure caused
  by an ambiguous or missing statement is a task defect and will come back to you.

### What the level-4 pilot showed (read before designing a level 5)

One attempt per task by a small model (README section 7, `results/pilot/FAILURE_AUDIT_l4.md`): six of
the seven level-4 tasks were passed in 2.5-17 minutes, although the level-3 method scores 0-0.65 on
every one of them. Element by element:

- **Q did not add difficulty.** With the criterion stated and the control in hand, the model iterated
  until the criterion held. A stated criterion plus a control is a checklist item. (The level-3 task
  that resisted, amplicon_16s, states no criterion.) Fairness and difficulty pull in opposite
  directions here.
- **P did not add difficulty.** None of the planted defaults was taken.
- **H did not add difficulty.** Re-running and comparing is what such a model does well, and the
  sound-but-odd decoy items were not flagged. H does expose ambiguous checklist wording: state in the
  SOP that every item is judged on its own and that a consequence belongs to the item that caused it,
  and check that no observation can be filed under two items.
- **T did not add difficulty.** The GUI route was used and took nine minutes.
- **What did resist** (both rounds): quantitative precision without a stated criterion; an artefact
  that needs one specific method; a process that takes real time to converge; and, at level 4, a
  step where the professional tool's default no longer answers the question, so the model wrote
  its own replacement and got a detail wrong, in a task with thirteen scored conclusions that all
  had to be right.

Design consequence: a harder level should be built from the last bullet (many independent
conclusions that must all be right; a necessary extension of a tool with no control to check
against; noise-limited precision; real compute with a time budget), not from more Q/P/H/T.

## Level 5 (third iteration: scaffolding removed)

Harbor runs with real flash models (README section 8) repeated the pilot's finding on a larger scale: every
task that states its acceptance criterion, or ships a checklist of what may be wrong, was solved; the two
tasks that defeated every model (somatic_tn, trio_liftover) state neither. Level 5 turns that observation
into an experiment: **the level-4 task of the same seed with the scaffolding taken away, and nothing else
changed.** It is not written by a family author; `alegen/descaffold.py` derives it from the level-4 build
(`python -m alegen make <family> --seed N --level 5`).

Scaffolding (removed):

| code | what | example |
|---|---|---|
| S1 | self-check criteria with thresholds | "accepted if the five medians lie within 50 units"; "at least 9 of the 10 reference genes" |
| S2 | procedures that say what to verify and how | "confirm the arm of every library from the data"; the marker-share formula of an identity control |
| S3 | audit checklists that enumerate the possible defects, and the audit deliverable | the eight-item hand-over checklists |
| S4 | statements of principle that warn about inherited material | "files left by a previous analyst are not raw data" |

What stays: the definition of every deliverable and of every term it uses (including the policy for a
case the analyst has to discover, e.g. "a library of mixed origin is excluded"); the facts in lab notes and
hand-over notes; the controls themselves (FMO tubes, reference genes, marker genes), present but unprompted.

Rules of the transformation (enforced in code):
- text only. Every edit must match exactly once (a changed generator is noticed); every input that is not
  an edited `.md` file is the same object as in the level-4 build; `tests/reference` passes through unchanged.
- if a deliverable disappears (the audit), its scorer component is removed and its weight is redistributed
  proportionally; attacks that only test that component are dropped. Nothing else in the scorer changes.
- all gates run again on the level-5 task. Ceilings that depended on the removed component are re-measured
  and the change is written into the task card (rnaseq: three "one defect kept" submissions score 0.79
  instead of at most 0.75, still 0.11 under the pass line).
- fairness no longer rests on a stated criterion. It rests on what the level-4 gates already established:
  the correct handling is standard practice, it is discoverable from the data with a standard QC step
  (`Trap.discoverable_by`), and two independent solutions agree with the planted truth.
- a family qualifies only if the truth stays unique without the criterion. hcs_dose_response does not:
  without the annotated fields and their tolerance a text-only agent has no way to tell whether a
  segmentation is right, which would be unfair rather than hard.

To compare level 4 and level 5 like with like, score the level-4 outputs again with the level-5 scorer
(`python -m alegen.rescore`).

## Level 7 (new silent decision points)

Levels 5 and 6 showed what stops a flash model and what does not: removing a stated acceptance criterion
does nothing; removing the sentence that tells the reader what to look at does (crispr_screen without the
"tubes were re-racked" remark: three models stopped checking sample identity). Where a model reads the
evidence unprompted (flow_compensation: voltages in the file headers, lots in the control table), text
removal is not enough, and a new decision point is needed. Its shape is fixed:

1. The routine workflow runs without an error and its result looks plausible.
2. It is wrong because one non-default professional step is missing.
3. No text names the problem or the step. The evidence is in the data or the file metadata, and a standard
   QC that practitioners run as a matter of course reveals it (`Trap.discoverable_by` names that QC and its
   measured effect size).
4. The existing definitions of the deliverables already determine the right answer; no sentence is added.

Mechanics: the family returns a level-7 build (`generate(seed, 7, workdir)`): level-4 content plus the new
decision point, with the level-4 texts unchanged. `python -m alegen make <family> --level 7` then applies the
level-6 text edits (`descaffold.apply(b, 7)`), so the family must not edit those texts. The validator
requires an attack `l6_reference` at level 7: the reference method of the level below on the level-7 inputs,
at most pass - 0.15. Specs are in `docs/specs/l7/`.

## G8: the difficulty gate

G1-G7 establish that an instance is correct. G8 measures whether it is hard for the model the benchmark is
aimed at (`python -m alegen.g8 <model>`, `alegen/g8.py`): an instance that passed G1-G7 is **kept** only if
the yardstick model fails every one of at least two attempts; it is **dropped** if the model passes every
attempt, and marked **unstable** if the outcomes are mixed. Repeated attempts are ordinary Harbor trials
whose job name ends in `__r<n>`.

- Every attempt, not one: the same flash model scored 1.0 and 0.60-0.66 on one task. Selecting on a single
  failure keeps unlucky draws, and the pass rate of the kept set drifts back up on the next run.
- G1-G7 first: without them a filter of this kind selects broken or ambiguous tasks, because those are
  failed most reliably. An instance that fails a gate is not released whatever the model does on it
  (crispr_screen level 6 seed 3: a weaker statistical test happens to pass on that seed -> rejected).
- The verdict is tied to the run conditions (agent, turn limit, time limit). A failure because the turn
  limit was reached is a failure under that limit; say so next to the verdict.
- Cost on this host: about one US dollar per attempt of gemini-3.8-flash, ten to twenty minutes each.

### What levels 5-7 and G8 showed (read before writing a new decision point)

- Removing a stated acceptance criterion does nothing (flow_compensation, crispr_screen at level 5).
- Removing the checklist item that names the right modelling choice does (rnaseq_confound level 5: the
  model found the exchanged samples and the mixed library without the SOP, and left the batch out of the design).
- Removing the incident report that prompts a check does (crispr_screen level 6: three models stopped
  checking sample identity; one of them validated its calls against the reference genes, which all agreed).
- Where a model reads the evidence unprompted from metadata (flow_compensation level 6: voltages in the file
  headers, lots in the control table), or the task is a fully specified computation (xlms_model_selection),
  or the standard workflow handles the peculiarity anyway (amplicon_16s level 6), text removal does nothing.
- A new silent decision point of the level-7 shape works where text removal failed (flow_compensation
  level 7: three models, none looked at the Time axis).
- It also works when the model sees the evidence (hcs_dose_response level 7, the plate turned on the stage:
  the yardstick noted vehicle wells without signal and readouts that rise with the concentration, explained
  them as an edge effect and a staining gradient, and reported every compound inactive, twice).
- It does not work where the standard workflow already looks at the place where it is planted
  (amplicon_16s level 7: the yardstick found the exchanged blank and specimen on both attempts while it
  inspected the blanks; another model took the roles from the sheet and lost exactly the affected components).
- A new decision point must leave alone what the level below relies on. hcs level 7 keeps the annotated well
  names, so the annotated fields became easy fields and the validation of the counting procedure lost its
  bite: the yardstick's counts met the criterion on them and failed on other wells. Such a side effect is
  reported next to the verdict (here the affected component alone would not have changed the outcome).
- Composing several instances into one task did not make it harder for the yardstick (project of three
  crispr_screen instances: 1.0).
- Not every family can be made harder this way: where the SOP clauses are the definition of the correct
  answer (md_system_repair), or the review tags are (variant_review), there is nothing to remove, and a new
  silent decision point would need a policy sentence, which is a pointer.

## Project tasks (several instances in one task)

`python -m alegen.project <family> --seed N --parts 3 --level 4` puts the level-4 (or level-5) instances of
seeds N..N+2 into one Harbor task: inputs under `/task/input/part_X/` with the part's own instruction as
`BRIEF.md`, deliverables under `/task/output/part_X/`. Every part is scored by its unchanged scorer on its
own reference; **the project score is the minimum**, and `details.json` keeps every part with all its
components (a dense reward for training, a strict pass for evaluation). The battery has the planted truth,
the within-tolerance submission and both solutions for all parts, the level-3 method on all parts, and
"all parts right except one" (must stay under the pass line: no passing on the average).

A project needs a per-family adapter (how to run the reference code on a part's directories). What it does
not have yet is coupling: the parts are independent instances. A cross-family project needs one tool image
for all families involved and a shared biological story in the generators.

## Instruction style (align with ALE task prompts)

One role sentence; `Task directory`; the inputs with paths and one-line descriptions; available
software; numbered task; **Required outputs** with exact filenames, column names, units, coordinate
conventions and definitions of every scored quantity; constraints (`Do not modify /task/input`,
`write only under /task/output`, `Do not ask for confirmation. Execute directly.`). Lab notes
(`*.md` in input) may state facts in domain jargon (as a real hand-over would) but never say what
to do about them.

## Scoring rules (`assets/score.py`)

- Use `_scoring.run(main, pass_threshold)`, `ctx.gate` for hard preconditions (score 0) and
  `ctx.component(name, value, weight)` with weights summing to 1. `CandidateError` = bad submission
  (score 0). Any other exception = verifier bug (no reward is written; that is intended).
- Never execute, import, unpickle or `eval` agent output. Cap sizes and row counts.
- Metrics must be robust to base rates and trivial answers (MCC or balanced measures instead of F1
  when one class dominates; penalise false positives; use the *minimum* over samples when losing a
  sample is a wrong deliverable).
- Every tolerance needs a measured justification: report the oracle's and the independent
  solution's raw metric values and the trap-ignoring values, and put the band between them with
  margin on both sides. Check at least two seeds.
- Give partial credit for real partial progress (earlier artefacts carry <= ~0.35 total weight),
  but a submission that ignores any single trap must stay below `pass_threshold`.
- Accept equivalent correct outputs (row order, extra columns, equivalent encodings).

## Battery (declared on the TaskBuild, run by `python -m alegen make`)

Probes (must score high):
- `planted_truth` - a submission written directly from the planted truth, floor 0.999 (when the
  deliverables can be written from truth).
- `within_tolerance` - truth perturbed inside the stated tolerance, floor 0.99.
- `independent_solution` - `script="SOLUTION=2 bash /solution/solve.sh"`, a genuinely different
  method/tool chain, floor 0.99 (lower only with a written reason, never below pass_threshold).

Attacks (must score low): the validator adds `empty_output`, `inputs_copied`, `zero_byte_files`,
`header_only`. You add:
- a prior-only / constant answer (ceiling <= 0.1),
- one "naive" variant **per trap** that ignores exactly that trap but still produces real outputs
  (ceiling <= pass_threshold - 0.15), plus the fully naive default pipeline,
- one perturbation just beyond tolerance.
A naive variant that crashes is not a test: make it produce the output a careless analyst would.
Use `Attack(script=..., script_files={...})`; helper files are mounted at `/alegen_script`.

## Running

```bash
cd /home/user/ale-lifesci-gym/pipeline
export PATH="$HOME/.local/bin:$PATH"
../.venv/bin/python -m alegen make <family> --seed 1 --level 3 --keep-work     # all gates incl. Harbor
../.venv/bin/python -m alegen make <family> --seed 2 --level 3 --no-harbor
../.venv/bin/python -m alegen make <family> --seed 1 --level 1 --no-harbor
```
Gates: G1 images build, G2 oracle >= 0.99, G3 probes, G4 attacks, G5 deterministic verifier,
G6 `harbor run -a oracle` = 1 and `-a nop` = 0, G7 the task's tools resolve in a login shell. Per-check details land in
`results/validation/<task_id>.json`; with `--keep-work`, each run's `details.json` is under
`results/work/<task_id>/`.

## Environment notes for this build host

- Docker daemon is already running and shared. Never restart it, never `docker system prune`,
  never remove images you did not create. Image builds must go through `python -m alegen make`
  or `dockerlib.build_image` (they add the proxy CA automatically). Containers run with
  `--network none`.
- Host downloads (reference genomes, structures) work through the proxy with plain `curl -fsSL`.
  Cache them under `/home/user/ale-lifesci-gym/.cache/<family>/` and reuse; record URL + sha256 in
  your design notes. Prefer small, stable sources (UCSC, Ensembl, RCSB, NCBI FTP).
- Generators run on the host venv (`../.venv/bin/python`, numpy available; add pure-python deps
  with `uv pip install --python ../.venv/bin/python <pkg>` and list them in your report). For
  steps that need domain tools at generation time use
  `tag = dockerlib.tool_image(conda=[...])` and `dockerlib.run_in_image(tag, script, mounts={hostdir: "/work"}, copy_out={...})`.
- Disk: on this host (containerd image store) every `docker build` leaves build cache of about the size of
  the image, and `docker builder prune -f` does not release it; `docker builder prune -af` does, without
  touching images. Run it after builds. A task whose `environment/Dockerfile` is byte-identical to that of
  a task whose image exists is built on top of that image by the validator (no network, no tool install).
- The machine has 4 cores and 16 GB RAM shared with other authors: use <= 2 threads for long
  steps and expect builds to be slow. Docker Hub pulls are rate limited: base every image on
  `python:3.12-slim` (already local) via the dockerlib helpers.

## Deliverables of a family author

1. The family code, passing the three `make` commands above (paste the gate lines in your report).
2. `docs/families/<family>.md`: scenario; truth model (what is planted, why it is correct); each
   trap with its discoverability evidence; output contract; scoring table with measured numbers
   (oracle, independent solution, each naive variant, two seeds); known limitations and anything
   you are unsure about. Be candid: a stated weakness is worth more than a hidden one.
3. `worklog.d/<family>.md`: dated entries (what was tried, what failed, what changed and why).
