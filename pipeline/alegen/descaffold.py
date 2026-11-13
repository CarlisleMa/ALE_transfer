"""Level 5 = a level-4 task with the scaffolding taken away. Data and planted truth are not touched.

Why: every task that stated its acceptance criterion, or shipped a checklist of what could be wrong, was
solved by the flash models; the tasks that defeated all of them (somatic_tn, trio_liftover) state neither.
Level 5 measures how much that scaffolding is worth, on the same inputs and the same truth.

What counts as scaffolding (removed):
  S1  self-check criteria with thresholds ("accepted if ... within 50 units", "at least 9 of 10 reference genes")
  S2  procedures that tell the analyst what to verify and how ("confirm the arm of every library from the data",
      the marker-share formula of an identity control)
  S3  audit checklists that enumerate the possible defects, and the audit deliverable that goes with them
  S4  statements of principle that warn about inherited files
What stays: the definition of every deliverable and of every term it uses, the facts in the lab notes and
hand-over notes, the controls themselves (FMO tubes, reference genes, marker genes) - present, but unprompted.

Level 6 = level 5 with the incident reports taken out of the lab notes as well (S5): the sentences in which
somebody wrote down that something went wrong or was changed ("the tubes were re-racked", "the voltage was
raised after sample X", "the labels had become illegible"). The flash models that passed level 5 did so by
reading exactly these sentences and then checking. At level 6 the evidence is only where it physically is:
in the data and in the file metadata (detector voltages in the FCS headers, reagent lots in the control
table, the sequencing run in the sample sheet, the sgRNA profiles and marker genes themselves). Protocol
descriptions stay; so do the predecessor's hand-over notes, which are part of the material to be judged.

The transformation is text-only and checked: every edit must match exactly once, every input file that is
not edited must be the same object as in the level-4 build, and tests/reference is passed through unchanged.
The scorer loses the audit component (its weight is redistributed proportionally) and nothing else.

  python -m alegen make <family> --seed 1 --level 5|6    (cli.py builds level 4 and applies this module)
"""
from __future__ import annotations

import re
from pathlib import Path

from .spec import TaskBuild

FAMILIES = ("flow_compensation", "crispr_screen", "rnaseq_confound")     # level 5 and 6
# Level the transformation starts from (the hardest level a family author wrote). Families that never had
# stated criteria or checklists have no level 5; their level 6 removes the pointer sentences directly.
TOP_LEVEL = {"flow_compensation": 4, "crispr_screen": 4, "rnaseq_confound": 4, "xlms_model_selection": 4,
             "amplicon_16s": 3, "hcs_dose_response": 4}


class Edit:
    """Text of one file with exact, counted replacements."""

    def __init__(self, text: str, what: str):
        self.t, self.what = text, what

    def sub(self, old: str, new: str, count: int = 1) -> "Edit":
        n = self.t.count(old)
        if n != count:
            raise ValueError(f"descaffold: {self.what}: expected {count} occurrence(s) of {old[:60]!r}, found {n}")
        self.t = self.t.replace(old, new)
        return self

    def cut(self, start: str, end: str | None) -> "Edit":
        """Remove from the line that starts with `start` up to (not including) the line that starts with `end`
        (or to the end of the text)."""
        i = self.t.find("\n" + start)
        if i < 0 or self.t.count("\n" + start) != 1:
            raise ValueError(f"descaffold: {self.what}: section start {start!r} not found exactly once")
        if end is None:
            self.t = self.t[:i].rstrip("\n") + "\n"
            return self
        j = self.t.find("\n" + end, i + 1)
        if j < 0:
            raise ValueError(f"descaffold: {self.what}: section end {end!r} not found after {start!r}")
        self.t = self.t[:i + 1] + self.t[j + 1:]
        return self

    def drop_line(self, pattern: str, repl: str = "") -> "Edit":
        """Remove (or replace by `repl`) the one line that matches `pattern` from its start."""
        self.t, n = re.subn(r"(?m)^" + pattern + r".*\n", repl, self.t)
        if n != 1:
            raise ValueError(f"descaffold: {self.what}: line pattern {pattern[:60]!r} matched {n} time(s), expected 1")
        return self

    def resub(self, pattern: str, repl: str, count: int = 1) -> "Edit":
        self.t, n = re.subn(pattern, repl, self.t, flags=re.S)
        if n != count:
            raise ValueError(f"descaffold: {self.what}: pattern {pattern[:60]!r} matched {n} time(s), expected {count}")
        return self


def _text(src) -> str:
    if isinstance(src, Path):
        return src.read_text()
    return src.decode() if isinstance(src, bytes) else src


def _drop_audit_from_scorer(score_py: str, block_start: str, block_end: str, weights: dict[str, tuple[str, str]],
                            what: str) -> str:
    """Remove the audit block of a scorer and rescale the remaining component weights (given as exact
    old -> new source snippets, so that a changed scorer is noticed)."""
    e = Edit(score_py, what)
    e.cut(block_start, block_end)
    for name, (old, new) in weights.items():
        e.sub(old, new)
    if 'ctx.table("audit.tsv"' in e.t:
        raise ValueError(f"descaffold: {what}: the scorer still reads audit.tsv")
    return e.t


# ----------------------------------------------------------------------------------------- flow_compensation
def _flow(b: TaskBuild) -> list[str]:
    sop = Edit(_text(b.input_files["gating_sop.md"]), "gating_sop.md")
    sop.cut("## 3. Compensation acceptance for the HLA-DR read-out", None)
    b.input_files["gating_sop.md"] = sop.t
    ins = Edit(b.instruction_md, "instruction")
    ins.sub("definitions,\n  the control-based HLA-DR gate and the compensation acceptance check",
            "definitions\n  and the control-based HLA-DR gate")
    b.instruction_md = ins.t
    return ["S1: the compensation acceptance criterion of the gating SOP (section 3: medians of five lymphocyte "
            "populations in the FMO tube within 50 detector units)"]


# --------------------------------------------------------------------------------------------- crispr_screen
def _crispr(b: TaskBuild) -> list[str]:
    sop = Edit(_text(b.input_files["screen_sop.md"]), "screen_sop.md")
    sop.resub(r"\n- The sample sheet is what the wet lab wrote down\. Before any contrast is computed.*?recorded in `qc\.json`\.\n",
              "\n")
    sop.sub("## 3. Counting and sample identity", "## 3. Counting")
    sop.cut("## 5. Release criterion for modifier calls", "## 6. Hand-over audit checklist")
    sop.cut("## 6. Hand-over audit checklist", "## 7. Erratum")
    sop.sub("## 7. Erratum", "## 5. Erratum")
    sop.sub("This\nSOP states what a released analysis has to satisfy. It does not prescribe software.",
            "This\nSOP defines the quantities of a released analysis. It does not prescribe software.")
    b.input_files["screen_sop.md"] = sop.t

    ho = Edit(_text(b.input_files["handover/HANDOVER.md"]), "handover/HANDOVER.md")
    ho.sub("Not done: the SOP checks, the reference genes, anything",
           "Not done: the reference genes, anything")
    b.input_files["handover/HANDOVER.md"] = ho.t

    ins = Edit(b.instruction_md, "instruction")
    ins.sub("You\ndeliver the final results, an audit of the inherited project and an erratum for the mailed list.",
            "You\ndeliver the final results and an erratum for the mailed list.")
    ins.sub("the facility SOP for drug-modifier screens: definitions, acceptance criteria,\n  the hand-over audit checklist",
            "the facility SOP for drug-modifier screens: definitions of the reported\n  quantities and of the erratum")
    ins.sub("3. Audit the inherited project against the checklist in section 6 of the SOP.\n4. Write the erratum",
            "3. Write the erratum")
    ins.sub(" The calls must satisfy the release criterion of SOP section 5.", "")
    ins.sub("  - Libraries are grouped by the treatment they actually received (SOP section 3).\n", "")
    ins.resub(r"- `audit\.tsv` - tab-separated, columns `item`, `verdict`; one row per checklist item of SOP\n  section 6, `verdict` is `ok` or `defective`\.\n",
              "")
    if "audit" in ins.t or "section 5" in ins.t or "section 6" in ins.t or "section 3" in ins.t:
        raise ValueError("descaffold: crispr instruction still refers to removed material")
    b.instruction_md = ins.t

    b.tests_files["score.py"] = _drop_audit_from_scorer(
        _text(b.tests_files["score.py"]), "    # ---------------- audit.tsv", "    # ---------------- erratum.tsv",
        {"guide_counts": ('ctx.component("guide_counts", ramp(min(per_sample), 0.50, 0.90), 0.15,',
                          'ctx.component("guide_counts", ramp(min(per_sample), 0.50, 0.90), 0.1875,'),
         "gene_lfc": ('ctx.component("gene_lfc", ramp(worst, 0.60, 0.25), 0.20,',
                      'ctx.component("gene_lfc", ramp(worst, 0.60, 0.25), 0.25,'),
         "drug_class": ('ctx.component("drug_class", ramp(m, 0.75, 0.90), 0.30,',
                        'ctx.component("drug_class", ramp(m, 0.75, 0.90), 0.375,'),
         "erratum": ('ctx.component("erratum", ramp(bal, 0.60, 0.85), 0.15,',
                     'ctx.component("erratum", ramp(bal, 0.60, 0.85), 0.1875,')},
        "crispr score.py")
    return ["S2: the instruction to confirm the arm of every day-14 library from the data, and the sentence that "
            "libraries are grouped by the treatment they actually received (SOP section 3, instruction)",
            "S1: the release criterion on the reference genes (SOP section 5); reference_genes.tsv stays in the inputs",
            "S3: the eight-item hand-over audit checklist (SOP section 6) and the deliverable audit.tsv; "
            "its weight (0.20) is redistributed proportionally over the other four components"]


# ------------------------------------------------------------------------------------------- rnaseq_confound
def _rnaseq(b: TaskBuild) -> list[str]:
    if "SOP_takeover.md" not in b.input_files:
        raise ValueError("descaffold: rnaseq: SOP_takeover.md not among the inputs")
    del b.input_files["SOP_takeover.md"]

    ho = Edit(_text(b.input_files["previous_analysis/HANDOVER.md"]), "previous_analysis/HANDOVER.md")
    ho.sub("- Samples were used as listed in the sample sheet. I did not get to the identity control of the SOP.\n",
           "- Samples were used as listed in the sample sheet.\n")
    b.input_files["previous_analysis/HANDOVER.md"] = ho.t

    ins = Edit(b.instruction_md, "instruction")
    ins.sub("You\nrelease the final result: gene-level counts, a per-sample QC table, the genes that respond to the\n"
            "compound, an audit of the previous analysis and an erratum for the mailed list.",
            "You\nrelease the final result: gene-level counts, a per-sample QC table, the genes that respond to the\n"
            "compound and an erratum for the mailed list.")
    ins.resub(r"- `SOP_takeover\.md` - the group's SOP for taking over an analysis \(identity control of libraries,\n"
              r"  audit checklist, erratum\)\. It is binding for this task\.\n", "")
    ins.sub("2. Check every sample (SOP section 2) and decide how it enters the statistical analysis.",
            "2. Decide how every sample enters the statistical analysis.")
    ins.sub("4. Audit the previous analysis against the checklist of SOP section 3.\n"
            "5. Tell the wet lab which of the mailed genes hold up (SOP section 4).",
            "4. Tell the wet lab which of the mailed genes hold up.")
    ins.sub("`mixed` for a library that fails the\n    identity control of the SOP.",
            "`mixed` for a library that contains\n    material of more than one donor.")
    ins.sub("or `excluded` if you left the sample out.",
            "or `excluded` if you left the sample out.\n"
            "    A library that contains material of more than one donor cannot be assigned to a record of the\n"
            "    sample sheet and is excluded; every other library enters the model under the condition its\n"
            "    cells actually received.")
    ins.resub(r"- `audit\.tsv` - tab-separated with at least the columns `item`, `verdict`; one row per checklist item\n"
              r"  of SOP section 3, `item` exactly as written there; `verdict` = `ok` or `defective`\.\n", "")
    if "audit" in ins.t or "SOP" in ins.t:
        raise ValueError("descaffold: rnaseq instruction still refers to removed material")
    b.instruction_md = ins.t

    b.tests_files["score.py"] = _drop_audit_from_scorer(
        _text(b.tests_files["score.py"]), "    # ---------------- audit.tsv", "    # ---------------- erratum.tsv",
        {"W": ('W = {"gene_counts": 0.08, "sample_sex": 0.02, "sample_labels": 0.07, "de_false_discoveries": 0.19,\n'
               '     "de_sensitivity": 0.16, "de_log2fc": 0.13, "audit": 0.20, "erratum": 0.15}',
               'W = {"gene_counts": 0.10, "sample_sex": 0.025, "sample_labels": 0.0875, "de_false_discoveries": 0.2375,\n'
               '     "de_sensitivity": 0.20, "de_log2fc": 0.1625, "erratum": 0.1875}')},
        "rnaseq score.py")
    # Without the audit component the other components weigh 1.25 x more, so a submission that keeps one planted
    # defect of the predecessor (and had no audit credit at level 4 either) scores higher than at level 4:
    # measured 0.7875 for the three below, 0.11 under the pass line. Their level-4 ceiling of 0.75 is raised to 0.80.
    for a in b.attacks:
        if a.name == "trust_inherited" or a.name.startswith("defect_kept_"):
            a.max_score = max(a.max_score, 0.80)
    return ["S4 + S2 + S1: the take-over SOP is not shipped: its principle about inherited files, the identity "
            "control with the marker-share formula and the 0.02 limit, and the sentence in the hand-over note that "
            "refers to it. The policy for a mixed library (reported as `mixed`, excluded) stays, as a definition of "
            "the deliverable; the marker genes stay in the gene panel",
            "S3: the eight-item audit checklist and the deliverable audit.tsv; its weight (0.20) is redistributed "
            "proportionally over the other seven components"]


def _flow6(b: TaskBuild) -> list[str]:
    n = Edit(_text(b.input_files["acquisition_notes.md"]), "acquisition_notes.md")
    n.drop_line(r"- Experiment template: created on ")
    n.drop_line(r"- After PBMC_\w+ the HLA-DR signal of the T cells was judged too close to the axis")
    n.drop_line(r"- PBMC_\w+: large events were piling up in the top FSC channel")
    n.sub(", stained the day before from the vials that were open at the time", "")
    n.sub(", taken from the vial that was used for the cocktail", "")
    b.input_files["acquisition_notes.md"] = n.t
    return ["S5: acquisition notes without the three incident reports (the template compensation that was never "
            "re-run; the detector voltage raised in mid-session and which controls were recorded again; the sample "
            "re-acquired with a lower FSC voltage) and without the remarks on which vials the controls were stained "
            "from. Voltages are in the FCS headers, reagent lots in controls.tsv and panel.tsv"]


def _crispr6(b: TaskBuild) -> list[str]:
    n = Edit(_text(b.input_files["screen_design.md"]), "screen_design.md")
    n.drop_line(r"- The day-14 gDNA tubes were re-racked")
    n.drop_line(r"- LIB\d+ failed QC on run A\.", "- Sequencing run B: NovaSeq, paired-end lane shared with other users.\n")
    b.input_files["screen_design.md"] = n.t
    ins = Edit(b.instruction_md, "instruction")
    ins.sub('"mapped_reads": <int>, "arm": "<T0|vehicle|drug>"}, ...}}`', '"mapped_reads": <int>}, ...}}`')
    ins.resub(r"assigned to a library sgRNA and `arm` the group in which you analysed the library\.",
              "assigned to a library sgRNA.")
    b.instruction_md = ins.t
    return ["S5: wet-lab notes without the two incident reports (the day-14 tubes re-racked before the sample sheet "
            "was filled in; the library that failed QC and came back as a Read-2 file). The sequencing run of every "
            "library is in the sample sheet",
            "the `arm` field of qc.json is no longer asked for (it was never scored; its description told the "
            "reader that the analysed group can differ from the sample sheet)"]


def _rnaseq6(b: TaskBuild) -> list[str]:
    n = Edit(_text(b.input_files["lab_notes.md"]), "lab_notes.md")
    n.drop_line(r"- Incident \(recorded for completeness\)")
    b.input_files["lab_notes.md"] = n.t
    return ["S5: wet-lab notes without the incident report (tube labels illegible, tubes re-labelled from the box map)"]


def _xlms6(b: TaskBuild) -> list[str]:
    n = Edit(_text(b.input_files["xlms_notes.md"]), "xlms_notes.md")
    n.resub(r"- The models come from the group's structure-prediction pipeline\. The prediction input were the expression "
            r"constructs.*?Both copies of a protein were predicted from the same construct\.\n",
            "- The models come from the group's structure-prediction pipeline.\n")
    n.sub(" They were pooled from several runs that were set up by different people.", "")
    b.input_files["xlms_notes.md"] = n.t
    ins = Edit(b.instruction_md, "instruction")
    ins.sub("1. Work out which residues of the models each cross-linked position corresponds to.\n2. For every model",
            "1. For every model")
    ins.sub("\n3. Rank the models by the number", "\n2. Rank the models by the number")
    ins.resub(r"- `mapping\.tsv` - tab-separated with the columns `accession`, `position`, `model_chain`,.*?present in the model\.\n",
              "")
    if "mapping" in ins.t:
        raise ValueError("descaffold: xlms instruction still refers to the mapping deliverable")
    b.instruction_md = ins.t
    sc = Edit(_text(b.tests_files["score.py"]), "xlms score.py")
    sc.cut("    # ---------------- mapping.tsv", "    # ---------------- components")
    sc.resub(r"\n    ok = 0\n    for tr in truth\[\"mapping\"\]:.*?otherwise 0\.6 x fraction\"\)\n", "\n")
    for old, new in (
            ('ctx.component("restraint_counts", 1.0 if exact == len(names) else 0.6 * frac, 0.30,',
             'ctx.component("restraint_counts", 1.0 if exact == len(names) else 0.6 * frac, 0.352941,'),
            ('ctx.component("link_status", ramp(mcc, 0.95, 0.995), 0.15,',
             'ctx.component("link_status", ramp(mcc, 0.95, 0.995), 0.176471,'),
            ('ctx.component("best_model", 1.0 if models[best]["near_native"] else 0.0, 0.30,',
             'ctx.component("best_model", 1.0 if models[best]["near_native"] else 0.0, 0.352941,'),
            ('ctx.component("ranking", ramp(auc, 0.70, 0.95), 0.10,',
             'ctx.component("ranking", ramp(auc, 0.70, 0.95), 0.117647,')):
        sc.sub(old, new)
    if "mapping.tsv" in sc.t or "mp.get" in sc.t:
        raise ValueError("descaffold: xlms scorer still uses mapping.tsv")
    b.tests_files["score.py"] = sc.t
    return ["S5: hand-over notes without the two remarks on where the models come from (predicted from trimmed "
            "constructs; pooled from several runs). That model numbering and chain letters differ from the "
            "full-length positions is visible in the model files and the FASTA sequences",
            "the instruction no longer lists 'work out which residues correspond' as a step, and the deliverable "
            "mapping.tsv is gone; its weight (0.15) is redistributed proportionally over the other four components"]


def _amplicon6(b: TaskBuild) -> list[str]:
    n = Edit(_text(b.input_files["library_prep_notes.md"]), "library_prep_notes.md")
    n.resub(r" The primer stocks are equimolar mixes of four phased versions per primer.*?\(heterogeneity spacers\)\.", "")
    n.resub(r"- Library construction: PCR products were cleaned up, end-repaired/A-tailed and the indexed Illumina adapters "
            r"were added by ligation \(PCR-free ligation kit\), one index pair per library\.\n",
            "- Library construction: indexed Illumina libraries, one index pair per library.\n")
    n.resub(r" The two \w+ strains of the standard differ at \d+ position\(s\) of the V4 amplicon\.", "")
    b.input_files["library_prep_notes.md"] = n.t
    return ["S5: wet-lab notes without the three sentences that describe what makes the reads unusual (phased primers "
            "with heterogeneity spacers; adapters added by ligation, i.e. inserts in both orientations; two strains of "
            "the mock standard that differ at one or two positions). All three are visible in the reads themselves and "
            "in the mock library"]


def _hcs6(b: TaskBuild) -> list[str]:
    n = Edit(_text(b.input_files["assay_notes.md"]), "assay_notes.md")
    n.drop_line(r"- Incubation: ", "- Incubation: 37 C / 5% CO2, standard polystyrene lid.\n")
    n.drop_line(r"- Staining: ", "- Staining: Alexa Fluor 488 conjugate against the marker protein, two PBS washes, "
                                 "wells left in 100 uL PBS.\n")
    n.drop_line(r"- Focus: ", "- Focus: hardware autofocus on the plate bottom.\n")
    n.drop_line(r"- Wells were not washed between treatment and fixation")
    b.input_files["assay_notes.md"] = n.t
    return ["S5: assay notes without the four remarks that explain where the plate artefacts come from (dry humidity "
            "tray and shelf position; stain added column by column and removed at once; one fixed focus offset and no "
            "review at the microscope; unwashed wells with floating material). The artefacts are visible in the vehicle "
            "wells of the plate map and in the images. The SOP, the annotated fields and the scorer are unchanged"]


_BY_FAMILY = {"flow_compensation": _flow, "crispr_screen": _crispr, "rnaseq_confound": _rnaseq}
_BY_FAMILY_6 = {"hcs_dose_response": _hcs6, "flow_compensation": _flow6, "crispr_screen": _crispr6, "rnaseq_confound": _rnaseq6,
                "xlms_model_selection": _xlms6, "amplicon_16s": _amplicon6}


def apply(b: TaskBuild, level: int = 5) -> TaskBuild:
    """Turn a level-4 build into the level-5 or level-6 build of the same seed (in place; returns it)."""
    if level not in (5, 6, 7):
        raise ValueError("descaffold produces level 5, 6 or 7")
    # level 7: a family author's build with new silent decision points and the level-4 texts; it gets the level-6 edits
    src = 7 if level == 7 else TOP_LEVEL.get(b.family)
    if b.level != src or b.family not in TOP_LEVEL or (level == 5 and b.family not in _BY_FAMILY):
        raise ValueError(f"level {level} of {b.family} is derived from its level-{src} build; got level {b.level} "
                         f"(level 5 exists for {FAMILIES}, level 6 for {tuple(TOP_LEVEL)})")
    before = dict(b.input_files)
    reference = {k: v for k, v in b.tests_files.items() if k != "score.py"}
    removed = _BY_FAMILY[b.family](b) if b.family in _BY_FAMILY else []
    if level >= 6:
        removed += _BY_FAMILY_6[b.family](b)
    # data and truth untouched: everything that is not a text edit is the very same object as before
    edited = sorted(k for k in before if k not in b.input_files or b.input_files[k] is not before[k])
    if any(not k.endswith(".md") for k in edited):
        raise AssertionError(f"descaffold changed a non-text input: {edited}")
    if any(b.tests_files[k] is not v for k, v in reference.items()) or set(b.tests_files) != set(reference) | {"score.py"}:
        raise AssertionError("descaffold changed the hidden reference")
    l4_id = b.task_id
    b.level = level
    b.task_id = re.sub(r"-l\d+$", f"-l{level}", l4_id)
    b.title = b.title + {5: " (level 5: scaffolding removed)", 6: " (level 6: pointers removed)",
                         7: " (level 7: pointers removed, new silent decision points)"}[level]
    b.summary = b.summary.rstrip() + f" Level {level}: same inputs and planted truth as level {src}; " \
                                     + ("without the stated acceptance criteria, verification procedures and audit checklist."
                                        if level == 5 else "the text no longer says what to look out for (see the task card).")
    b.attacks = [a for a in b.attacks if not a.name.startswith("audit_")]
    b.extra_metadata = dict(b.extra_metadata, derived_from=l4_id,
                            source_level=(TOP_LEVEL[b.family] if level == 7 else src), scaffolding_removed=removed,
                            text_inputs_edited=edited)
    b.judge_explanation_md = (
        f"**Level {level}.** Inputs and planted truth are those of `" + l4_id + "`; only text was taken away:\n\n"
        + "\n".join(f"- {r}" for r in removed)
        + "\n\nFairness no longer rests on a stated criterion but on three things that are checked for the level-4 "
          "task this one is derived from: the correct handling is standard practice of the field, it is discoverable "
          "from the data with a standard QC step (see the decision points above), and two independent solutions "
          "agree with the planted truth.\n\n" + b.judge_explanation_md)
    return b
