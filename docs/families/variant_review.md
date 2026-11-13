# Family `variant_review` - manual review of somatic variant candidates (IGV)

Written by the maintainer; the GUI-capable family.

## Scenario
A caller proposed 48 somatic candidates for one tumour/normal pair. The analyst reviews each one
against the reads following the laboratory's review SOP (evidence definition + seven failure tags
with precedence) and delivers calls, tags, allele fractions and an IGV session.

## Truth model
BAM files are synthesised read by read from a 60 kb slice of hg38 chr22 (UCSC, sha256 in
`task.toml`), so the generator knows what every read carries. Each candidate gets exactly one
planted mechanism with parameters far from every SOP threshold:
true somatic SNV (clonal 0.25-0.45, subclonal 0.09-0.16), true somatic deletion (unique
representation), true deletion in a homopolymer, germline het, 1-2 supporting molecules,
2-3 molecules amplified into flagged duplicates, strand-specific artefact, read-end artefact,
low base quality, low mapping quality, paralog-like reads with 5-7 further mismatches.
`truth.json` = planted call and tag, plus the realised molecule-level allele fraction.

## Decision points and discoverability
| trap | effect if ignored | evidence |
|---|---|---|
| `artefact_classes` | counts-only review calls 17 artefacts somatic (MCC 0.47) | SOP tags; per-read attributes in IGV or pysam |
| `duplicate_only_support` | read-level counting sees 8-15 balanced reads and calls 4 false variants | flag 0x400; IGV hides duplicates by default; SOP counts molecules |
| `indel_representation` | position-exact counting finds 0 support for 4 true deletions | the gap is drawn at the right end of the repeat (see `docs/img/igv_repeat_deletion.png`); SOP: judge the sequence, not the gap position |

## GUI
The image ships IGV 2.19.8 and the desktop MCP server (`alegen/desktop/`). `/task/software/open_igv.sh`
starts IGV offline with the reference, both BAMs and the candidate track; verified by screenshot
and by driving IGV through its batch port and through the MCP tools. Grading is artefact-based, so
a terminal-only agent can solve it; the IGV session file is checked for structure only (weight 0.05).

## Scoring and measured values (level 3; seeds 1, 2 and 3 give identical numbers because the class
counts are fixed and the truth is exact)
| submission | calls MCC | tag accuracy | somatic with correct VAF | score |
|---|---|---|---|---|
| oracle (pysam, haplotype comparison) | 1.0 | 1.0 | 1.0 | 1.0 |
| independent solution (samtools-view text, own CIGAR walker, left-normalisation) | 1.0 | 1.0 | 1.0 | 1.0 |
| one wrong call + VAF noise 0.008 | 0.957 | 0.9667 | 1.0 | 0.990 |
| counts-only review | 0.4721 | 0.4333 | 1.0 | 0.25 |
| duplicates counted | 0.8421 | 0.8667 | 1.0 | 0.655 |
| position-exact indels | 0.8284 | 1.0 | 0.7778 | 0.634 |
| every candidate somatic | 0.0 | 0.0 | 1.0 | 0.25 |
| every candidate failed as LOW_SUPPORT | 0.0 | 0.2667 | 0.0 | 0.05 |

(In-container values from `results/work/variant-review-s0001-l3/*/logs/details.json`.)

## Levels
1: germline and low-support only (a counts-only review is correct: scores 1.0). 2: + five artefact
classes. 3: + duplicate-only support and repeat-context deletions.

## Small-model pilot (not a Harbor run)
One attempt by Claude Haiku through `docker exec`, terminal only (the desktop MCP server was not
connected to the pilot harness), answers sealed and audit clean (`results/pilot/AUDIT_batchB.txt`):
reward **1.0** in 4 min 6 s with 13 tool calls. It counted read pairs as molecules, excluded
duplicate-flagged reads and decided deletion support by comparing the junction sequence, which
handles the repeat-context deletions. Reading: with an explicit SOP this is a careful-execution
task that a small model completes; the prior (0.78, the highest of the ten) is wrong.

## Known limitations
- The SOP is explicit, so this is a protocol-execution task: it trains careful read-level work and
  rule precedence, not the judgement calls of real manual review, where evidence is borderline.
- Planted artefacts are textbook-clean (100% one strand, all at read ends). Real artefacts are
  mixed; a model trained only on this would be over-confident.
- Whether IGV was actually used cannot be verified from artefacts; the session check can be
  satisfied by writing XML by hand. It is worth 0.05 and cannot decide pass/fail.
- Mates never overlap in the simulated library, so molecule vs read counting only matters through
  duplicates.

## Level 4 - review on a locked-down workstation where IGV is the only reader (T + P)

Added 2026-11-12 after the pilot (level 3 was solved in 4 min from the terminal with pysam, IGV never
opened). Code: `alegen/families/variant_review/level4.py`, assets `sop_review.py`, `igv_gui.py`,
`review_cli.py`, `score_l4.py`, `solve_l4.sh`, `harden_igv.py`; desktop layer `alegen/desktop/cli.py`.
Levels 1-3 are untouched (`diff -r` of a regenerated level 3 seed 1 against `tasks/variant-review-s0001-l3`
is empty).

### What was added and why
| element | what | why |
|---|---|---|
| `T_igv_only_reader` | Alignments are CRAM 3.1 (reference-compressed, written with htslib 1.24). The image holds IGV 2.19.8, the desktop layer and Python with the standard library, numpy and the MCP package: no samtools, pysam, htslib, compiler. The Java runtime is re-linked with the 27 modules IGV needs; IGV's bundled libraries lose what could be used as a reader (see "Bypasses"). The instruction says in bold that IGV is the only software that can read the alignments and names its three interfaces (desktop application, batch port, igvtools). | the level-3 pilot never opened IGV |
| `P_mnv_same_molecule` | 4 pairs of adjacent substitutions on the same molecules (SOP tag 9 `MNV_PARTNER`: one event, the right-hand candidate fails) and 4 pairs on disjoint molecules with the same number of molecules on each strand (both stay somatic). `igvtools count` shows identical per-strand counts at the two positions in both situations. | per-variant instead of per-haplotype logic |
| `P_discordant_pair_support` | 4 candidates with 11-13 supporting molecules of which all but one (91-92%) are read pairs with a 1.2-1.9 kb insert (SOP tag 8 `DISCORDANT_PAIRS`, threshold 1,000 bp; ordinary pairs measure 232-477 bp). Counts, strands, qualities and positions in the read look like a clonal variant. | support confined to anomalous molecules |
| `P_normal_trace` | 4 of the 12 true SNVs have exactly one alternate molecule in the normal (1.6-3.0% of the normal depth; planted germline variants: 41-67%); the germline rule needs 2 molecules and 5%. | "any alternate read in the normal = germline" |
| mixed artefacts | every artefact class has one clean supporter (except strand bias) and carries other blemishes below their thresholds. Realised over both seeds: low base quality 89-92% of supporters (rule: 75%), low mapping quality 89-92%, read end 90-93%, multi-mismatch 5.7-6.2 further mismatches per read (rule: 4), strand bias 100% (rule: 90%); the same quantities in every other class with at least 4 supporters are <= 25%, <= 20%, <= 26%, <= 1.2 and <= 75%. True SNVs carry 1-6 low-quality (5-25%) and 0-4 poorly mapped (0-12%) supporters, have >= 8 supporters and <= 71% of them on one strand | level-3 artefacts were textbook-clean |
| documentation | `review_session.xml` must be written by IGV (track panels with `clazz=...AlignmentTrack`, `RenderOptions colorOption="READ_STRAND"` in both alignment tracks, view on the first candidate, <= 200 bp) and there is one IGV image per candidate (`snapshots/<id>.png`, >= 800x300, shows data, not filed twice). | the level-3 session check accepted five hand-typed lines |
| `desktop` command | `alegen/desktop/cli.py`, installed as `/usr/local/bin/desktop`: `screenshot <png> [--region X,Y,W,H]`, `click`, `move`, `drag`, `type`, `key`, `scroll`, `launch`, `windows`, `wait` - the actions of the MCP server for shell-only agents and harnesses. The MCP server (`server.py`, unchanged) still works in the level-4 image (tools listed, `list_windows` and `screenshot` called over stdio). | the pilot harness had no MCP client |
| `igvtools` launcher, `igv/README.md` | the bioconda IGV package has neither the `igvtools` launcher nor the igvtools readme, and there is no internet: a launcher and a reference of the batch commands and igvtools options are provided as software documentation (every listed item was run in the hardened image except `exit` and `--strands first`; `group base` answers "position is required" when the position, which the reference shows as optional, is left out) | fairness: tool documentation is not a recipe |
| `open_igv.sh` waits | IGV's batch port answers ~5 s before the start-up tracks are loaded (measured: port at 4.6 s, tracks at 9.6 s). The script now returns when they are loaded. | found when the first oracle run wrote empty images; it would have hit any agent |

The generator re-derives the SOP outcome of every candidate from the reads it wrote (with `sop_review.py`
on SAM text produced by pysam on the host) and refuses a seed if a call differs from the planted one or if any
quantity lies within its safety margin of a threshold (75% rules: must be <= 40% or >= 85%; mean further
mismatches <= 2 or >= 5; one-strand fraction <= 80% or >= 97%; shared-molecule fraction <= 0.2 or >= 0.98;
>= 6 supporters for anything that is not `LOW_SUPPORT`; normal fraction >= 12% for germline). Seeds 1-16 all
generate (seeds 3 and 5 only after the level-4 placement fallback to 380 bp spacing).

### Routes through IGV that were run
Both routes were run end to end in the hardened image, from a fresh container, and scored by the verifier.

1. **Command-line route = oracle** (`review_cli.py`). `igvtools count -w 1 --bases --strands read` on each
   CRAM (duplicates are excluded by default) and once more on the tumour with `--minMapQuality 20`: depth,
   alternate molecules per strand, deletions (the DEL column summed over every equivalent placement of the
   gap), the normal, the mapping-quality rule. What counts cannot give is read from SVG images of a 221-bp
   view around each candidate written through the batch port (`preference SAM.COMPUTE_ISIZES false`,
   `maxPanelHeight 20000`, `goto`, `snapshot <id>.svg`): each read is a polygon (arrow side = strand), each
   mismatch a rectangle in the base colour whose fill opacity is IGV's base-quality shading (opaque from
   quality 20 = the SOP threshold), a deletion is a line between two blocks, a pair with an insert above
   1,000 bp is filled `rgb(200,0,0)`. For all 68 candidates of both seeds the number of supporting reads
   found in the SVG equals the igvtools count. Session and images: `colorBy READ_STRAND`, `goto`,
   `saveSession`, `snapshot`.
2. **GUI route = independent solution** (`igv_gui.py` + `sop_review.py`), only `desktop` actions: type a locus
   into the locus box, right-click the alignment track, "Export Alignments..." (IGV's file dialog, full path
   typed into "Enter file name") for three overlapping 20-kb views and both tracks, apply the SOP to the
   exported SAM text, "Color alignments by > read strand" in both tracks, File > Save Session..., then for each
   candidate locus box + File > Save PNG Image .... Each step is checked on its result (window titles, a
   column of screen pixels through both tracks, the written file) and repeated if the result is missing.
   Evidence: `docs/img/variant_review_l4_gui_track_menu.png` (the track menu with "Export Alignments..."),
   `docs/img/variant_review_l4_gui_after_export.png` (status bar "180 reads written" at a discordant-pair
   candidate; the red read carries the alternate base).

A third way was not automated but is what a vision-capable agent would do: look at the images.

### Bypasses (anything that reads the alignments without IGV's own interfaces)
Each line was tried in the unmodified image (bioconda `igv=2.19.8` + conda openjdk 21, see worklog) and again in the
shipped image: `pipeline/alegen/families/variant_review/dev/image_checks.sh`, output of the final seed-1 image in
`results/variant_review_l4_image_checks.log`. The same script checks that IGV, igvtools, the batch commands and
igvtools options listed in `igv/README.md`, the `desktop` command and the MCP server work.

| bypass | unmodified image | shipped image |
|---|---|---|
| samtools / pysam / htslib / compilers | not installed | not installed (`perl` of the base image is present, without any sequence library) |
| pure-Python BAM reader (BAM is gzip blocks + fixed binary records: zlib and struct suffice) | n/a | avoided by shipping CRAM 3.1 instead of BAM |
| `java Foo.java` / `javac` / `jshell` against htsjdk | `javac`, `jshell` and the `jdk.compiler` module are there (a full JDK); not run | closed: runtime re-linked without `jdk.compiler`, `jdk.jshell`; `java Foo.java` and `java -m jdk.jshell/...` fail (`Module ... not found`) |
| `htsjdk.samtools.example.PrintReadsExample` and the other htsjdk example programs | present in htsjdk-4.3.0; my one try started and stopped at a host-name lookup (container without network), not at reading | closed: classes removed |
| Nashorn JavaScript shell bundled with IGV (`org.openjdk.nashorn.tools.Shell`): script with full access to htsjdk | works (printed reads) | closed: nashorn jar removed, `jdk.dynalink` and `java.scripting` not in the runtime |
| `jrunscript` / `com.sun.tools.script.shell.Main` with the Nashorn engine (commons-jexl registers a script engine the same way; not tried) | works (printed reads) | closed: `java.scripting` not in the runtime |
| `org.apache.commons.jexl2.Main` (expression evaluator) | runs (evaluated `1+1`); reaching htsjdk through it was not tried | closed: class removed |
| `igvtools sort` with the CRAM under a `.bam` name and `JAVA_TOOL_OPTIONS=-Dsamjdk.reference_fasta=...` -> SAM text | works (15 MB of SAM) | closed: `BAMSorter`/`SAMSorter` classes removed (README says so) |
| any other class with a `main` in IGV's jars (e.g. `org.broad.igv.tools.PairedEndStats`) | present | closed: `main(String[])` renamed in 107 class files, kept only in `org.broad.igv.ui.Main` and `org.broad.igv.tools.IgvTools` |
| a class file assembled by hand (or a copy of a jar with the renaming undone) that calls htsjdk | not tried | **open**: `java` and htsjdk must stay for IGV; needs a class-file writer in Python |
| pure-Python CRAM 3.1 decoder written from memory (rANS Nx16, adaptive arithmetic coding, name tokeniser) | not tried | **open** |
| installing a reader | no network, package cache cleaned | same |

Not claimed: that no other entry point exists. What was verified is the list above.

### Measured values
In-container values from `results/validation/variant-review-s000{1,2}-l4.json` and
`results/work/.../logs/details.json`. Columns: MCC of calls / fraction of planted failures with the planted tag /
fraction of true variants called with VAF within 0.05 / valid images / score. Pass mark 0.9; careless variants must
stay <= 0.75.

| submission | seed 1 | seed 2 |
|---|---|---|
| oracle (command-line route: igvtools count + batch port) | 1.0 / 1.0 / 1.0 / 68 / **1.0** | 1.0 / 1.0 / 1.0 / 68 / **1.0** |
| independent solution (GUI route) | 1.0 / 1.0 / 1.0 / 68 / **1.0** | 1.0 / 1.0 / 1.0 / 68 / **1.0** |
| `planted_truth` (table from truth, stand-in session and images) | 1.0 / 1.0 / 1.0 / 68 / 1.0 | 1.0 / 1.0 / 1.0 / 68 / 1.0 |
| `within_tolerance` (one wrong call, VAF noise 0.008) | 0.9707 / 0.9737 / 1.0 / 68 / 1.0 | 0.9705 / 1.0 / 0.9667 / 68 / 1.0 |
| `review_credit_without_documentation` = `handwritten_session_blank_snapshots` | 1.0 / 1.0 / 1.0 / 0 / 0.85 | 1.0 / 1.0 / 1.0 / 0 / 0.85 |
| `l3_reference` (level-3 `review.py`, pysam) | no output / 0.0 | no output / 0.0 |
| `l3_logic_on_igv_export` (level-3 SOP on exported reads; ignores both new tags) | 0.7895 / 0.7895 / 1.0 / 68 / 0.30 | 0.7895 / 0.7895 / 1.0 / 68 / 0.30 |
| `igvtools_counts_only` (real igvtools, no switches, no documentation) | 0.432 / 0.3421 / 1.0 / 0 / 0.15 | 0.432 / 0.3421 / 1.0 / 0 / 0.15 |
| `counts_only_review` | 0.432 / 0.3421 / 1.0 / 68 / 0.30 | 0.432 / 0.3421 / 1.0 / 68 / 0.30 |
| `ignore_mnv_phase` | 0.8885 / 0.8947 / 1.0 / 68 / 0.5205 | 0.8885 / 0.8947 / 1.0 / 68 / 0.5205 |
| `merge_every_adjacent_pair` | 0.8855 / 1.0 / 0.8667 / 68 / 0.6097 | 0.8855 / 1.0 / 0.8667 / 68 / 0.6097 |
| `ignore_discordant_pairs` | 0.8885 / 0.8947 / 1.0 / 68 / 0.5205 | 0.8885 / 0.8947 / 1.0 / 68 / 0.5205 |
| `any_alt_read_in_normal_is_germline` | 0.8855 / 1.0 / 0.8667 / 68 / 0.6097 | 0.8855 / 1.0 / 0.8667 / 68 / 0.6097 |
| `duplicates_counted` | 0.8885 / 0.8947 / 1.0 / 68 / 0.5205 | 0.8885 / 0.8947 / 1.0 / 68 / 0.5205 |
| `position_exact_indels` | 0.8855 / 1.0 / 0.8667 / 68 / 0.6097 | 0.8855 / 1.0 / 0.8667 / 68 / 0.6097 |
| `trust_the_caller` (all somatic) | 0.0 / 0.0 / 1.0 / 68 / 0.30 | 0.0 / 0.0 / 1.0 / 68 / 0.30 |
| `fail_everything` | 0.0 / 0.2105 / 0.0 / 68 / 0.15 | 0.0 / 0.2105 / 0.0 / 68 / 0.15 |
| `two_wrong_calls` | 0.9404 / 0.9737 / 0.9667 / 68 / 0.89 | 0.9404 / 0.9737 / 0.9667 / 68 / 0.89 |
| `vaf_beyond_tolerance` (+0.07) | 1.0 / 1.0 / 0.0 / 68 / 0.85 | 1.0 / 1.0 / 0.0 / 68 / 0.85 |
| generic: empty, inputs copied, zero-byte files, header only | 0.0 | 0.0 |

The two seeds give the same numbers for every careless variant because the class counts are fixed and the truth is
exact; they differ only where the battery draws a random candidate (`within_tolerance`). The review tables of the two
routes are identical row by row (call, tag, all four counts) for both seeds.

`l3_reference` is the unmodified level-3 `review.py`: it stops at `import pysam` and writes nothing.
`igvtools_counts_only` runs the real `igvtools count --bases` in the image without `--strands` and
`--minMapQuality` and writes no session or images; its table is identical to `counts_only_review` (computed on the
host from the planted reads) for all 68 rows, i.e. igvtools counts and molecule counts agree everywhere.
`review_credit_without_documentation` (probe, floor 0.85) and `handwritten_session_blank_snapshots` (attack,
ceiling 0.85) are the same submission: a correct table with a hand-typed session and blank images scores exactly 0.85.

Bands (`score_l4.py`): calls 0.40 (MCC 0.86 -> 0.965), failure tags 0.30 (0.85 -> 0.97), allele fraction 0.15
(0.80 -> 0.96), session 0.07, images 0.08. If the review components reach < 0.9 of their weight the total is capped at
0.89 (`two_wrong_calls`: review 0.8895, documentation full, capped to 0.89). One wrong call still gives 1.0
(MCC 0.9707, tags 0.9737), two wrong calls fail.

### Runtime, size, reliability
| quantity | measured |
|---|---|
| agent image (`docker images`) | 1.84 GB; verifier image 178 MB |
| `environment/input` | 412 KB (tumour CRAM 217 KB, normal CRAM 100 KB) |
| generation on the host | 10-11 s per seed |
| oracle (command-line route), `runtime_min` of gate G2 | 2.21 / 2.16 min in the first complete validations (host load average 10-14 from other authors' builds), 1.52 / 0.89 min in the final ones; 45 s / 40 s when re-run alone (`results/work/<task>/oracle_stdout_rerun.log`) |
| GUI route | 12 min 10 s in the first complete run (load 13-14, before the number of process starts per action was reduced); 6 min 43 s alone from a fresh IGV start: 63 s until the session is saved, then 5.0 s per image; in the first validation probes 4.5 s (seed 1) and 4.2 s (seed 2) per image, from the file times |
| `python -m alegen make` | seed 1 with Harbor 13.6 min and 10.2 min (final); seed 2 without Harbor 12.5 min and 6.7 min (final) |

Reliability of the automation, counted over everything that was run (level 4 only):

| route | before the fix named in the worklog | after it |
|---|---|---|
| command-line route (oracle) | 1 of 2: correct in a container where IGV was already running, **0.0 in the first `make`** (fresh container, batch commands sent before IGV had loaded the tracks) | 7 of 7 scored 1.0: 1 alone (seed 1), gate G2 twice and Harbor oracle twice on seed 1, gate G2 twice on seed 2; plus 2 unscored re-runs (one per seed) that completed with the planted tag distribution |
| GUI route (independent solution) | 0 of 1: the first attempt stopped at the first export (`xdotool click --repeat 3` selected one word of the locus box, IGV raised a message box) | 6 of 6 complete runs scored 1.0: seed 1 four times (two alone, two as validation probe), seed 2 twice (validation probe). The two runs whose output I kept log no repeated step; the validator does not keep the output of probes, so for the other four only the result is known |

Each seed went through two complete validations; `results/validation/*.json` are the final ones (made after the last
change to the code; regenerating both tasks with the final code gives directories identical to `tasks/variant-review-s000{1,2}-l4`).

### Deviations from the spec
- **Oracle and second solution are swapped relative to my first plan**: the oracle is the command-line route
  (no pixel positions, about 2 min), the GUI route is the independent solution. Both are run by every validation.
- **Python is not "standard library only"**: gate G7 of the validator imports numpy in a login shell, so numpy
  stays; the `mcp` package and its dependencies stay for the MCP server. None of them reads alignments. scipy and
  pandas are gone (own Dockerfile string instead of `dockerlib.agent_dockerfile`).
- **Soft clipping was not planted**; "support confined to anomalous molecules" is implemented with insert size
  only. Reasons: one exactly defined quantity (TLEN > 1,000) that both routes see (red reads at IGV's default
  1,000-bp threshold; TLEN in the export), and soft-clipped reads would have made the SOP's "end of the read"
  ambiguous for `READ_END`.
- **"Variant at low level in the normal" is a decoy, not a failure class**: one alternate molecule, which the
  SOP's germline rule (2 molecules and 5%) does not catch; a third value between "absent" and "germline" would have
  needed a new threshold close to both.
- **`review.tsv` has no new columns**; the new classes are expressed by two new tags.
- **igvtools "duplicate switch"**: `igvtools count` excludes duplicates by default (`--includeDuplicates` adds
  them), so the counts-only attack is "no strand, no mapping-quality, no read-level data", not "duplicates counted";
  `duplicates_counted` is a separate attack.
- `oracle_steps` is 14 + 4 x candidates = 286 commands of the command-line route (3 igvtools runs, batch
  commands), recorded as counted.

### Known limitations
- **The verifier cannot prove that IGV wrote the session and the images.** The `planted_truth` probe passes
  with a session string modelled on IGV's output and generator-drawn PNGs. An image of the wrong locus, or of a
  view without alignments but with >= 40 distinct pixel rows, also counts. Documentation is worth 0.15 and
  cannot lift a failing review; the review table cannot be produced without reading the alignments.
- **The two open bypasses above.** A model that writes a JVM class file or a CRAM 3.1 decoder by hand solves the
  task without IGV's interfaces.
- **The GUI reference solution is tied to pixel positions** of IGV 2.19.8 at 1280x800 as opened by
  `open_igv.sh` (track menu, File menu, file dialog). It checks every step and retries, but a different window
  geometry would need new positions. This affects the reference solution, not what an agent has to do.
- The SOP is still explicit: level 4 adds tool use and per-molecule reasoning, not judgement under ambiguity.
- Multi-nucleotide events are pairs only; mates still never overlap; the discordant pairs differ from the
  library in insert size only (same orientation). Found while measuring: the 8 duplicate-flagged copies of
  discordant-pair reads per seed carry the proper-pair flag (the copy step does not carry the flag over). They are
  duplicates and add no evidence, and the SOP uses the insert size, not the flag; not fixed because it would change
  the validated bytes.
- `SAM.DOWNSAMPLE_READS=false` is set by `open_igv.sh`. By the numbers IGV's default (100 reads per 50 bp) should
  not drop reads either (at most 71 and 76 non-duplicate reads start in any 50-bp window of the tumour for seeds 1
  and 2; the generator refuses seeds above 95), but a run with down-sampling switched on was not made.
- The read-level part of the command-line route depends on how IGV 2.19.8 draws SVG (polygons, opacity for base
  quality, `rgb(200,0,0)` for long inserts). That is how the reference solution reads them, not an interface promise.
- Agent pilots were not run (by instruction): whether the task is hard for a small model is not measured. The
  oracle shows that a shell-only path exists (igvtools counts plus SVG images; `review_cli.py` has 373 lines); a
  vision-capable agent has the shorter path through the track menu's alignment export.
- The generator can still refuse a seed (margin checks, placement). Seeds 1-16 generate; nothing beyond 16 was tried.
