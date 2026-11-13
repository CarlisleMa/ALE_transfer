"""Family: choosing among predicted models of a two-protein complex with XL-MS restraints.

Truth-by-construction: a real heterodimer from RCSB is the hidden native. Every shipped model is a
rigid-body variant of it (one subunit moved relative to the other), so its interface RMSD / fnat /
DockQ to the native are known exactly. Cross-links are *sampled from the native geometry* (plus
planted false positives), so which model the restraints support is decided by the generator and
then measured on the coordinates that are actually written to disk. The native itself, its PDB id,
chain ids and numbering never reach the agent image; accessions are blinded.

Decision points:
  level 2  numbering_offset_and_gap   cross-link positions are full-length positions, the models
                                      are numbered 1..n per chain over trimmed constructs
           chain_identity_varies      which chain letter is which protein differs per model
  level 3  misleading_confidence      the pipeline's top-ranked model is wrong
           false_positive_crosslinks  10-15% of the reported residue pairs are false
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np

from ... import dockerlib
from ...spec import Attack, Probe, TaskBuild, Trap

ASSETS = Path(__file__).parent / "assets"
CACHE = Path(__file__).resolve().parents[4] / ".cache" / "xlms_model_selection"
RCSB = "https://files.rcsb.org/download/{pdb}.pdb"

# Hidden natives: (PDB id, the two chains used, sha256 of the cached download).
ENTRIES = [
    {"pdb": "1FIN", "chains": ("A", "B"),
     "sha256": "12f7920d265a08c0145152e6f960d942a17dd947585fab0470869111e8136fca"},
    {"pdb": "1EFU", "chains": ("A", "B"),
     "sha256": "f621d20b05abfc7550a5250a971023cbfa11a49fefaf53946153674d11b687c0"},
    {"pdb": "1KXP", "chains": ("A", "D"),
     "sha256": "5d68ab0d9bec2230c7240a77988cd00a4cadff3dd97d814c60a85e6182d76418"},
    {"pdb": "1IBR", "chains": ("A", "B"),
     "sha256": "181e03e9c5b5122ed91b91600927f017286a0cabb43575cb0a36550aa6a65c14"},
    {"pdb": "1I2M", "chains": ("A", "B"),
     "sha256": "d11d60e277912fc62f2379ae4d2a64fb5d7970dbc9742197a11ea6f4cf4a096c"},
]
# cross-linker, full name, spacer arm (A), the group's CA-CA cutoff (A)
LINKERS = [("DSSO", "disuccinimidyl sulfoxide", 10.3, 30.0),
           ("BS3", "bis(sulfosuccinimidyl) suberate", 11.4, 30.0),
           ("DSBU", "disuccinimidyl dibutyric urea", 12.5, 32.0)]
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
       "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
       "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "MSE": "M"}
AA1 = {v: k for k, v in AA3.items() if k != "MSE"}
BACKBONE = ("N", "CA", "C", "O")
MIN_ISLAND = 12        # modelled segments shorter than this are trimmed away (keeps sequence mapping unambiguous)
MARGIN = 0.6           # every evaluated distance is at least this far from the cutoff (requirement: 0.5 A)
COUNT_GAP = 3          # every near-native model satisfies >= this many more inter-links than any wrong model
NEAR_IRMSD = 2.0       # near-native: interface RMSD below this
WRONG_IRMSD = 6.0      # wrong: interface RMSD above this (and DockQ < 0.23)


# ----------------------------------------------------------------------------- native structure

def _fetch(entry: dict) -> Path:
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"{entry['pdb']}.pdb"
    if not path.is_file():
        tmp = path.with_suffix(".part")
        subprocess.run(["curl", "-fsSL", "--max-time", "120", "-o", str(tmp), RCSB.format(pdb=entry["pdb"])],
                       check=True)
        tmp.rename(path)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != entry["sha256"]:
        raise RuntimeError(f"{path} has sha256 {digest}, expected {entry['sha256']}")
    return path


def _load_chain(path: Path, chain: str) -> dict:
    """One chain of the native: SEQRES sequence and the resolved residues (complete backbone only)."""
    seqres: list[str] = []
    res: dict[int, dict] = {}
    for line in path.read_text().splitlines():
        rec = line[:6]
        if rec == "SEQRES" and line[11] == chain:
            seqres += line[19:70].split()
        elif rec == "ENDMDL":
            break
        elif rec in ("ATOM  ", "HETATM") and line[21] == chain:
            resn = line[17:20]
            if resn not in AA3 or (rec == "HETATM" and resn != "MSE"):
                continue
            if line[26] != " ":
                raise RuntimeError("insertion codes are not supported")
            name, elem = line[12:16].strip(), line[76:78].strip().upper()
            if elem in ("H", "D") or name.startswith("H"):
                continue
            if resn == "MSE" and name == "SE":
                name, elem = "SD", "S"
            num = int(line[22:26])
            r = res.setdefault(num, {"aa": AA3[resn], "atoms": {}, "b": []})
            if name in r["atoms"] or name == "OXT":      # first alternate location only
                continue
            r["atoms"][name] = (elem or name[0], [float(line[30:38]), float(line[38:46]), float(line[46:54])])
            r["b"].append(float(line[60:66]))
    seq = "".join(AA3.get(x, "X") for x in seqres)
    nums = sorted(res)
    best = None
    for c in range(-max(nums) - 1, len(seq)):
        if all(0 <= n + c < len(seq) and seq[n + c] == res[n]["aa"] for n in nums):
            best = c
            break
    if best is None:
        raise RuntimeError(f"{path.name} chain {chain}: residue numbering is not a constant shift of SEQRES")
    out = []
    for n in nums:
        r = res[n]
        if not all(b in r["atoms"] for b in BACKBONE):
            continue
        names = list(r["atoms"])
        out.append({"idx": n + best, "aa": r["aa"], "names": names,
                    "elems": [r["atoms"][a][0] for a in names],
                    "xyz": np.array([r["atoms"][a][1] for a in names]), "b": float(np.mean(r["b"]))})
    return {"seqres": seq, "res": out}


# ----------------------------------------------------------------------------- geometry helpers

def _rot(axis: np.ndarray, ang: float) -> np.ndarray:
    a = axis / np.linalg.norm(axis)
    k = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + np.sin(ang) * k + (1 - np.cos(ang)) * (k @ k)


def _unit(rng) -> np.ndarray:
    v = rng.normal(size=3)
    return v / np.linalg.norm(v)


def _superpose(p: np.ndarray, q: np.ndarray):
    """Rotation r and translation t that put p onto q (least squares): p @ r.T + t ~ q."""
    pc, qc = p.mean(0), q.mean(0)
    u, _, vt = np.linalg.svd((p - pc).T @ (q - qc))
    d = np.sign(np.linalg.det(vt.T @ u.T))
    r = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    return r, qc - pc @ r.T


def _rmsd(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(((a - b) ** 2).sum() / len(a)))


def _fit_rmsd(p: np.ndarray, q: np.ndarray) -> float:
    r, t = _superpose(p, q)
    return _rmsd(p @ r.T + t, q)


def _resmin(xa: list[np.ndarray], xb: list[np.ndarray]) -> np.ndarray:
    """Minimum heavy-atom distance for every residue pair (len(xa) x len(xb))."""
    allb = np.vstack(xb)
    starts = np.cumsum([0] + [len(x) for x in xb[:-1]])
    out = np.empty((len(xa), len(xb)))
    for i, a in enumerate(xa):
        d = np.sqrt(((a[:, None, :] - allb[None, :, :]) ** 2).sum(-1)).min(0)
        out[i] = np.minimum.reduceat(d, starts)
    return out


# ----------------------------------------------------------------------------- sequence layout

_DISORDER = "SGPEKADTQNRSGPEASTQDGKPLVE"   # composition of the synthetic N-terminal regions


def _nterm(rng, n: int) -> str:
    return "M" + "".join(_DISORDER[i] for i in rng.integers(0, len(_DISORDER), size=n - 1))


def _segments(idx: list[int]) -> list[list[int]]:
    segs, cur = [], [idx[0]]
    for a in idx[1:]:
        if a == cur[-1] + 1:
            cur.append(a)
        else:
            segs.append(cur)
            cur = [a]
    segs.append(cur)
    return segs


def _ambiguous(full: str, modelled: set[int]) -> set[int]:
    """0-based positions whose assignment to 'modelled' vs 'missing' is not unique because an internal
    gap could be slid along a repeat. Lysines there are never used for cross-links."""
    bad: set[int] = set()
    order = sorted(modelled)
    for a0, b0 in zip(order, order[1:]):
        if b0 == a0 + 1:
            continue
        a, b = a0 + 1, b0 - 1            # the gap, inclusive
        for s in range(1, 40):
            if a - s >= 0 and full[a - s:a] == full[b - s + 1:b + 1]:
                bad.update(range(a - s, a))
                bad.update(range(b - s + 1, b + 1))
            if b + s < len(full) and full[a:a + s] == full[b + 1:b + 1 + s]:
                bad.update(range(a, a + s))
                bad.update(range(b + 1, b + 1 + s))
    return bad


def _plant_deletion(rng, prot: dict, partner_ca: np.ndarray) -> tuple[int, int]:
    """Pick an internal segment (surface, away from the interface, with a lysine if possible) that is
    removed from the modelled construct. Returns inclusive residue-list indices."""
    res, full = prot["res"], prot["full"]
    ca = np.array([r["xyz"][r["names"].index("CA")] for r in res])
    nb = (np.linalg.norm(ca[:, None] - ca[None], axis=2) < 10.0).sum(1)
    far = np.linalg.norm(ca[:, None] - partner_ca[None], axis=2).min(1) > 14.0
    n = len(res)
    cands = []
    for ln in range(6, 15):
        for s in range(MIN_ISLAND, n - ln - MIN_ISLAND):
            e = s + ln - 1
            fp = [res[i]["pos"] for i in range(s - MIN_ISLAND, e + MIN_ISLAND + 1)]
            if fp[-1] - fp[0] != len(fp) - 1:           # needs >= MIN_ISLAND contiguous residues both sides
                continue
            if not 0.15 <= s / n <= 0.60 or not far[s:e + 1].all():
                continue
            a, b = res[s]["pos"] - 1, res[e]["pos"] - 1  # 0-based in full
            if full[a - 1] == full[b] or full[a] == full[b + 1]:
                continue
            nk = sum(1 for i in range(s, e + 1) if res[i]["aa"] == "K")
            cands.append((nk > 0, -float(nb[s:e + 1].mean()), s, e))
    if not cands:
        raise RuntimeError("no place for the trimmed loop")
    with_k = [c for c in cands if c[0]] or cands
    with_k.sort(key=lambda c: c[1], reverse=True)          # most exposed first
    top = with_k[:max(1, len(with_k) // 4)]
    _, _, s, e = top[int(rng.integers(0, len(top)))]
    return s, e


# ----------------------------------------------------------------------------- PDB writer

def _pdb_text(chains: list[tuple[str, list[dict], np.ndarray, list[int], np.ndarray]]) -> str:
    """chains: (label, residues, coordinates of all atoms, residue numbers, per-residue pLDDT)."""
    out, serial = [], 1
    for label, residues, xyz, numbers, plddt in chains:
        k = 0
        for r, num, pl in zip(residues, numbers, plddt):
            for name, elem in zip(r["names"], r["elems"]):
                nm = name if len(name) == 4 else " " + name.ljust(3)
                x, y, z = xyz[k]
                out.append(f"ATOM  {serial:5d} {nm} {AA1[r['aa']]} {label}{num:4d}    "
                           f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{pl:6.2f}          {elem:>2s}")
                serial += 1
                k += 1
        last = residues[-1]
        out.append(f"TER   {serial:5d}      {AA1[last['aa']]} {label}{numbers[-1]:4d}")
        serial += 1
    out.append("END")
    return "\n".join(out) + "\n"


def _read_ca(text: str) -> dict[tuple[str, int], np.ndarray]:
    ca = {}
    for line in text.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            ca[(line[21], int(line[22:26]))] = np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return ca


# ----------------------------------------------------------------------------- instruction / notes

INSTRUCTION = """\
You are the integrative structural biologist of a protein-complex group. A structure-prediction run
produced a set of candidate models of a two-protein complex, and the mass-spectrometry facility
delivered lysine-lysine cross-links (XL-MS) measured on the same complex in solution. The group
needs every model evaluated against the cross-links and wants to know which model's arrangement of
the two proteins the data support best.

Task directory:
- `/task`

Inputs (read-only) under `/task/input`:
- `models/model_NN.pdb` - {n_models} predicted models of the complex (`model_01` ... `model_{n_models:02d}`),
  PDB format, two protein chains each
{ranking_line}- `crosslinks.csv` - the identified cross-links, one row per unique cross-linked residue pair:
  `xl_id`, `protein_a`, `position_a`, `protein_b`, `position_b`, `n_csm`, `score`. Proteins are
  given by accession, positions are residue positions in the full-length protein sequence.
- `sequences.fasta` - full-length sequences of the two proteins (the accession is the first word of
  the header). Protein identities were blinded.
- `xlms_notes.md` - hand-over notes of the mass-spectrometry facility and the modelling group

Available software (on `PATH`): Python 3.12 with numpy, scipy, pandas and Biopython. There is no
internet access.

Your task:
1. Work out which residue of the models each cross-linked position corresponds to.
2. For every model, count the cross-links it satisfies, separately for inter-protein and
   intra-protein cross-links.
3. Rank the models by how well the cross-linking data support the arrangement of the two proteins
   relative to each other, and select the best-supported model.

Required outputs under `/task/output`:
- `restraints.tsv` - tab-separated with at least the columns `model`, `inter_satisfied`,
  `inter_total`, `intra_satisfied`, `intra_total`; one row per model; `model` is the file name
  without `.pdb`; the four counts are integers. Definitions:
  - Every row of `crosslinks.csv` is one cross-link. It is inter-protein if `protein_a` and
    `protein_b` are different proteins and intra-protein if they are the same protein.
  - A cross-link is evaluable in a model if both cross-linked residues are present in that model.
    `inter_total` / `intra_total` are the numbers of evaluable cross-links of each class. A
    cross-link with a residue that is not present in the model is counted neither as satisfied nor
    in the total.
  - An evaluable cross-link is satisfied if the two cross-linked residues meet the distance
    criterion that `xlms_notes.md` gives for this experiment.
- `selection.json` - `{{"best_model": "model_NN", "ranking": ["model_NN", ...]}}`. `ranking` lists
  every model exactly once, from the model whose inter-protein arrangement is best supported by the
  cross-links to the one that is least supported; its first element is `best_model`.
- `mapping.tsv` - tab-separated with the columns `accession`, `position`, `model_chain`,
  `model_resnum`; one row for every distinct cross-linked position (accession + position) that
  occurs in `crosslinks.csv`. `model_chain` and `model_resnum` are the chain identifier and the
  residue number of that residue in `models/model_01.pdb`; write `NA` in both columns if the
  residue is not present in the model.

Constraints:
- Do not modify anything under `/task/input`; write deliverables only under `/task/output`.
- Every model must be reported.
- Do not ask for confirmation. Execute directly.
"""

JUDGE = """\
**Where the truth comes from.** The hidden native is a real heterodimer from RCSB. Every shipped
model is the native with one subunit moved as a rigid body, so interface RMSD, fnat and DockQ of
each model are known from the transform itself. Cross-links were sampled from lysine pairs that are
within range in the native (false positives: lysine pairs that are far apart in the native). The
per-model satisfied counts in `reference/truth.json` are computed with the planted residue mapping
on the coordinates exactly as written to the model files; the generator refuses to emit an instance
unless every evaluated distance is >= 0.6 A away from the cutoff, every near-native model satisfies
at least 3 more inter-protein cross-links than any wrong model, and there is no model between the
two quality classes (near-native: interface RMSD < 2 A; wrong: interface RMSD > 6 A and DockQ < 0.23).

**Hard gates (score 0).** Missing/unparseable file; `restraints.tsv` not listing every model exactly
once with integer counts; `selection.json` without a `best_model` that is a model, without a
`ranking` that is a permutation of all models, or with `ranking[0] != best_model`; duplicated
positions in `mapping.tsv`.

**Components.**

| weight | component | metric | full credit | zero |
|---|---|---|---|---|
| 0.30 | restraint_counts | fraction of models whose four counts are all exact | all models exact (otherwise 0.6 x fraction) | no model exact |
| 0.35 | best_model | hidden interface RMSD of `best_model` | near-native (< 2 A) | any other model |
| 0.20 | ranking | ROC AUC of the submitted order for near-native versus wrong models | >= 0.95 | <= 0.70 |
| 0.15 | residue_mapping | fraction of cross-linked positions with exact chain + residue number in model_01 (or NA) | every position exact (otherwise 0.6 x fraction) | no position exact |

**Why these tolerances.** Counts are exact because the definition is exact and no distance is
within 0.6 A of the cutoff (PDB coordinates are rounded to 0.001 A), so any correct implementation
reproduces them; a single wrong model already costs 0.30 x (1 - 0.6 x (N-1)/N) > 0.12 and fails the
task. The mapping is exact for the same reason (lysines whose assignment would depend on where an
aligner opens a gap are never cross-linked); a missing or useless mapping table costs 0.15 and
fails the task. Any near-native model is accepted as `best_model` because near-native models can
tie on the satisfied count. The ranking is scored only on near-native versus wrong pairs: the evidence cannot
order wrong models among themselves, and tie-breaking among equally supported models must not
matter. A count-based ranking has AUC 1.0 by construction; the prediction pipeline's own confidence
order is planted with AUC <= 0.60.

**What is deliberately not scored.** The order among wrong models and among near-native models; any
distance table, plot or report.
"""


def _notes(level: int, linker: tuple, accs: list[str], fdr_pct: int) -> str:
    name, longname, spacer, cutoff = linker
    lines = [
        "# XL-MS hand-over notes", "",
        "## Sample and cross-linking",
        f"- Sample: heterodimer of {accs[0]} and {accs[1]} (identities blinded for this hand-over), both "
        "full-length, co-expressed, purified by affinity and size-exclusion chromatography; 1:1 by SEC-MALS.",
        f"- Cross-linker: {name} ({longname}), amine-reactive, spacer arm {spacer:.1f} A. Only Lys-Lys "
        "cross-links are reported (protein N-termini, Ser/Thr/Tyr side reactions and mono-links were not "
        "considered).",
        "- Digest: trypsin, SCX enrichment of cross-linked peptides, LC-MS/MS on an Orbitrap instrument.",
        "",
        "## Identification",
        "- Database: the two sequences in `sequences.fasta` plus common contaminants. Positions in "
        "`crosslinks.csv` are positions in these sequences.",
        "- The table is collapsed to unique residue pairs. `n_csm` = number of cross-link spectrum matches "
        "supporting the pair, `score` = best match score (higher is better).",
    ]
    if level >= 3:
        lines.append("- Error control: target-decoy, list cut at 5% FDR on the spectrum-match level, which the "
                     f"facility estimates to correspond to about {fdr_pct}% FDR on the residue-pair level. Decoy "
                     "hits have been removed from the table; there was no manual validation.")
    else:
        lines.append("- Error control: target-decoy at 1% residue-pair FDR, followed by manual validation of "
                     "every spectrum; pairs that did not pass manual validation were removed.")
    lines += [
        "",
        "## Distance criterion used by the group",
        f"- For {name} the group counts a cross-link as satisfied by a structure if the CA-CA distance of the "
        f"two linked lysines is at most {cutoff:.1f} A, and as violated otherwise.",
        "- The sample is a 1:1 heterodimer: a pair of lysines of the same protein is treated as an "
        "intra-protein cross-link, a pair of lysines of the two different proteins as an inter-protein "
        "cross-link.",
        "- A cross-link to a residue that is absent from a structure cannot be evaluated for that structure; "
        "the group leaves it out of both the satisfied count and the total.",
        "",
        "## Models",
        "- The models come from the group's structure-prediction pipeline. The prediction input were the "
        "expression constructs used for crystallisation trials, which are shorter than the full-length "
        "proteins (disordered terminal regions and flexible internal loops had been trimmed).",
    ]
    if level >= 2:
        lines.append("- The model files are the unmodified output of the pipeline. They were pooled from "
                     "several runs that were set up by different people.")
    else:
        lines.append("- The modelling group renumbered every model to the full-length sequence numbering and "
                     "gave the two proteins the same chain identifiers in all models.")
    if level >= 3:
        lines.append("- `ranking.json` holds the pipeline's ranking confidence (0.8 x ipTM + 0.2 x pTM); the "
                     "model files are named in the order of that confidence.")
    else:
        lines.append("- The models are in no particular order.")
    return "\n".join(lines) + "\n"


# ----------------------------------------------------------------------------- generator

def generate(seed: int, level: int, workdir: Path) -> TaskBuild:
    if level >= 4:      # level 4 lives in its own module; levels 1-3 below are unchanged
        from . import l4
        return l4.generate(seed, level, workdir)
    rng = np.random.default_rng(7000 + seed)
    entry = ENTRIES[(seed - 1) % len(ENTRIES)]
    linker = LINKERS[int(rng.integers(0, len(LINKERS)))]
    cutoff = linker[3]
    path = _fetch(entry)

    # ---- proteins: full-length sequence = synthetic N-terminal region + SEQRES of the native chain
    alnum = "ABCDEFGHJKLMNPQRSTUVWXYZ23456789"
    accs = []
    while len(accs) < 2:
        a = "XL" + "".join(alnum[i] for i in rng.integers(0, len(alnum), size=4))
        if a not in accs:
            accs.append(a)
    prots = []
    for k, ch in enumerate(entry["chains"]):
        c = _load_chain(path, ch)
        ext = _nterm(rng, int(rng.integers(14, 95)))
        res = [r for seg in _segments([r["idx"] for r in c["res"]]) if len(seg) >= MIN_ISLAND
               for r in c["res"] if r["idx"] in set(seg)]
        for r in res:
            r["pos"] = len(ext) + r["idx"] + 1                      # 1-based full-length position
        prots.append({"acc": accs[k], "full": ext + c["seqres"], "res": res})
    ca_of = [np.array([r["xyz"][r["names"].index("CA")] for r in p["res"]]) for p in prots]
    for k, p in enumerate(prots):
        s, e = _plant_deletion(rng, p, ca_of[1 - k])
        p["deleted"] = set(range(s, e + 1))
        p["model_idx"] = [i for i in range(len(p["res"])) if i not in p["deleted"]]
        p["modelled_pos0"] = {p["res"][i]["pos"] - 1 for i in p["model_idx"]}
        p["ambiguous"] = _ambiguous(p["full"], p["modelled_pos0"])
        # model residue number of every modelled residue
        if level >= 2:
            p["resnum"] = {p["res"][i]["pos"]: n + 1 for n, i in enumerate(p["model_idx"])}
        else:
            p["resnum"] = {p["res"][i]["pos"]: p["res"][i]["pos"] for i in p["model_idx"]}

    # ---- lysines
    # status: "model" (in the models), "native" (resolved in the native but trimmed from the models),
    #         "none" (no coordinates anywhere: N-terminal region or unresolved in the native)
    lys = []     # (protein index, full position, status, CA or None)
    for k, p in enumerate(prots):
        by_pos = {r["pos"]: (i, r) for i, r in enumerate(p["res"])}
        for j, aa in enumerate(p["full"]):
            if aa != "K" or j in p["ambiguous"]:
                continue
            pos = j + 1
            if pos in by_pos:
                i, r = by_pos[pos]
                lys.append((k, pos, "native" if i in p["deleted"] else "model", r["xyz"][r["names"].index("CA")]))
            else:
                lys.append((k, pos, "none", None))
    L = len(lys)
    nat_d = np.full((L, L), np.nan)
    for a in range(L):
        for b in range(a + 1, L):
            if lys[a][3] is not None and lys[b][3] is not None:
                nat_d[a, b] = nat_d[b, a] = float(np.linalg.norm(lys[a][3] - lys[b][3]))

    def pairs(cond):
        out = [(a, b) for a in range(L) for b in range(a + 1, L) if cond(a, b)]
        return [out[i] for i in rng.permutation(len(out))]

    used: dict[int, int] = {}
    links: list[dict] = []
    taken: set[tuple[int, int]] = set()

    def take(cands, n, kind, cap=4):
        got = 0
        for a, b in cands:
            if got >= n:
                break
            if (a, b) in taken or used.get(a, 0) >= cap or used.get(b, 0) >= cap:
                continue
            taken.add((a, b))
            used[a] = used.get(a, 0) + 1
            used[b] = used.get(b, 0) + 1
            links.append({"a": a, "b": b, "kind": kind})
            got += 1
        return got

    both_model = lambda a, b: lys[a][2] == "model" and lys[b][2] == "model"          # noqa: E731
    inter = lambda a, b: lys[a][0] != lys[b][0]                                      # noqa: E731
    n_inter = take(pairs(lambda a, b: both_model(a, b) and inter(a, b) and nat_d[a, b] <= cutoff - 4.0),
                   int(rng.integers(13, 18)), "true_inter", cap=3)
    if n_inter < 11:
        raise RuntimeError("not enough inter-protein lysine pairs in range")
    take(pairs(lambda a, b: both_model(a, b) and not inter(a, b) and abs(lys[a][1] - lys[b][1]) >= 2
               and nat_d[a, b] <= cutoff - 1.5), int(rng.integers(24, 33)), "true_intra")
    take(pairs(lambda a, b: {lys[a][2], lys[b][2]} == {"model", "native"} and nat_d[a, b] <= cutoff - 4.0),
         3, "true_unresolved", cap=5)
    n_true = len(links)
    # one reported pair reaches into a region that has no coordinates at all (N-terminal region or a
    # segment that is disordered in the native); it can never be evaluated, whatever its origin
    take(pairs(lambda a, b: "none" in (lys[a][2], lys[b][2]) and "model" in (lys[a][2], lys[b][2])),
         1, "false_unresolved" if level >= 3 else "disordered_region", cap=5)
    if level >= 3:
        n_fp = max(5, int(round(n_true * 0.14)))
        n_fp_inter = max(2, int(round(0.4 * n_fp)))
        take(pairs(lambda a, b: both_model(a, b) and inter(a, b) and nat_d[a, b] >= cutoff + 8.0),
             n_fp_inter, "false_inter", cap=5)
        take(pairs(lambda a, b: both_model(a, b) and not inter(a, b) and nat_d[a, b] >= cutoff + 3.0),
             n_fp - n_fp_inter - 1, "false_intra", cap=5)
    else:
        # a manually validated list can still hit a lysine that is missing from the models
        take(pairs(lambda a, b: {lys[a][2], lys[b][2]} == {"model", "native"} and nat_d[a, b] <= cutoff - 4.0),
             1, "true_unresolved", cap=6)
    fdr = sum(1 for x in links if x["kind"].startswith("false")) / len(links)

    # ---- rigid-body poses of protein 2 (protein 1 stays in the native frame)
    bb = [np.vstack([r["xyz"][[r["names"].index(a) for a in BACKBONE]] for r in p["res"]]) for p in prots]
    heavy = [[r["xyz"] for r in p["res"]] for p in prots]
    resmin = _resmin(heavy[0], heavy[1])
    if_res = [np.where((resmin < 10.0).any(1))[0], np.where((resmin < 10.0).any(0))[0]]
    if_bb = [np.vstack([bb[k][4 * i:4 * i + 4] for i in if_res[k]]) for k in (0, 1)]
    contacts = np.argwhere(resmin < 5.0)
    nat_if = np.vstack(if_bb)
    c1, c2 = ca_of[0].mean(0), ca_of[1].mean(0)
    u0 = (c2 - c1) / np.linalg.norm(c2 - c1)
    if_centre = if_bb[1].mean(0)
    ev_inter = [x for x in links if x["kind"] in ("true_inter", "false_inter")]
    k1 = np.array([lys[x["a"] if lys[x["a"]][0] == 0 else x["b"]][3] for x in ev_inter])
    k2 = np.array([lys[x["b"] if lys[x["a"]][0] == 0 else x["a"]][3] for x in ev_inter])

    def assess(r, t):
        """Pose = protein 2 coordinates x @ r.T + t. Returns (iRMSD, satisfied inter count, margin)."""
        mod_if = np.vstack([if_bb[0], if_bb[1] @ r.T + t])
        d = np.linalg.norm(k1 - (k2 @ r.T + t), axis=1)
        return _fit_rmsd(mod_if, nat_if), int((d <= cutoff).sum()), float(np.abs(d - cutoff).min())

    near = []
    for _ in range(600):
        r = _rot(_unit(rng), np.deg2rad(rng.uniform(1.0, 9.0)))
        t = if_centre - if_centre @ r.T + _unit(rng) * rng.uniform(0.0, 1.6)
        irmsd, cnt, mar = assess(r, t)
        if 0.5 <= irmsd <= NEAR_IRMSD - 0.1 and mar >= MARGIN:
            near.append((irmsd, cnt, r, t))
    n_models = int(rng.integers(22, 29))
    n_near = 3
    near.sort(key=lambda x: x[0])
    targets = np.linspace(0.75, 1.75, n_near)
    chosen_near = []
    for tg in targets:
        best = min((x for x in near if not any(x is y for y in chosen_near)), key=lambda x: abs(x[0] - tg))
        chosen_near.append(best)
    min_near = min(x[1] for x in chosen_near)

    # Wrong poses: protein 2 in a random orientation, pushed against a different part of the surface
    # of protein 1 until the two chains share about as many CA-CA contacts (< 8 A) as the native
    # interface does, so that contact number or compactness do not give the near-native models away.
    wrong = []
    ca_mod = [np.array([p["res"][i]["xyz"][p["res"][i]["names"].index("CA")] for i in p["model_idx"]])
              for p in prots]
    nat_contacts = int((np.linalg.norm(ca_mod[0][:, None] - ca_mod[1][None], axis=2) < 8.0).sum())
    ca2c = ca_mod[1] - c2
    for _ in range(1600):
        perp = np.cross(u0, _unit(rng))
        theta = np.deg2rad(rng.uniform(15, 55) if rng.random() < 0.35 else rng.uniform(55, 180))
        u = _rot(perp, theta) @ u0
        q = _rot(_unit(rng), np.deg2rad(rng.uniform(10, 180)))
        b = ca2c @ q.T
        d = (ca_mod[0] - c1)[:, None, :] - b[None, :, :]
        p_ = d @ u
        q2 = (d ** 2).sum(-1) - p_ ** 2
        m = q2 < 64.0
        reach = np.sort(p_[m] + np.sqrt(64.0 - q2[m]))[::-1]     # slide at which each CA pair comes within 8 A
        want = int(rng.uniform(0.6, 1.6) * nat_contacts)
        if len(reach) < want:
            continue
        s = float(reach[want - 1])
        t = c1 + s * u - c2 @ q.T
        irmsd, cnt, mar = assess(q, t)
        if irmsd >= WRONG_IRMSD and mar >= MARGIN and cnt <= min_near - COUNT_GAP:
            wrong.append((irmsd, cnt, q, t))
    n_wrong = n_models - n_near
    bins: dict[int, list] = {}
    for w in wrong:
        bins.setdefault(w[1], []).append(w)
    chosen_wrong = []
    while len(chosen_wrong) < n_wrong:
        progressed = False
        for cnt in sorted(bins, reverse=True):
            if bins[cnt] and len(chosen_wrong) < n_wrong:
                chosen_wrong.append(bins[cnt].pop(int(rng.integers(0, len(bins[cnt])))))
                progressed = True
        if not progressed:
            raise RuntimeError("not enough wrong poses")

    def quality(r, t):
        mod_if = np.vstack([if_bb[0], if_bb[1] @ r.T + t])
        irmsd = _fit_rmsd(mod_if, nat_if)
        kept = sum(1 for i, j in contacts
                   if np.sqrt(((heavy[0][i][:, None, :] - (heavy[1][j] @ r.T + t)[None, :, :]) ** 2).sum(-1)).min() < 5.0)
        fnat = kept / len(contacts)
        rec, lig = (0, 1) if len(prots[0]["res"]) >= len(prots[1]["res"]) else (1, 0)
        mod = [bb[0], bb[1] @ r.T + t]
        rr, tt = _superpose(mod[rec], bb[rec])
        lrmsd = _rmsd(mod[lig] @ rr.T + tt, bb[lig])
        dockq = (fnat + 1 / (1 + (irmsd / 1.5) ** 2) + 1 / (1 + (lrmsd / 8.5) ** 2)) / 3
        return {"irmsd": round(irmsd, 3), "fnat": round(fnat, 3), "lrmsd": round(lrmsd, 3), "dockq": round(dockq, 3)}

    poses = [{"near": True, "r": x[2], "t": x[3], **quality(x[2], x[3])} for x in chosen_near]
    for x in chosen_wrong:
        qd = quality(x[2], x[3])
        if qd["dockq"] >= 0.23 or qd["fnat"] >= 0.3:
            raise RuntimeError("a wrong pose is not clearly wrong")
        poses.append({"near": False, "r": x[2], "t": x[3], "pre_count": x[1], **qd})
    for p in poses[:n_near]:
        if p["irmsd"] >= NEAR_IRMSD or p["dockq"] < 0.49:
            raise RuntimeError("a near-native pose is not near-native")

    # ---- order of the models = file names. The name order never favours the near-native models:
    # at level 3 it is the pipeline's confidence order, below it is arbitrary.
    N = len(poses)
    pre = [(p.get("pre_count", 99)) for p in poses]
    while True:
        order = [int(i) for i in rng.permutation(N)]
        if level >= 3:      # the most convincing wrong model gets the top confidence
            top = max(range(n_near, N), key=lambda i: pre[i])
            order.remove(top)
            order.insert(0, top)
        rank = {m: i for i, m in enumerate(order)}
        auc = np.mean([rank[a] < rank[b] for a in range(n_near) for b in range(n_near, N)])
        if order[0] >= n_near and 0.30 <= auc <= 0.60:
            break
    poses = [poses[i] for i in order]
    names = [f"model_{i + 1:02d}" for i in range(N)]
    conf = np.sort(rng.uniform(0.36, 0.89, size=N))[::-1]
    conf[0] = max(conf[0], conf[1] + 0.02)

    # ---- chain labelling per model
    if level >= 2:
        schemes = [((0, 1), ("A", "B")), ((1, 0), ("A", "B")), ((0, 1), ("B", "C")), ((1, 0), ("B", "C"))]
        first = schemes[int(rng.integers(0, 2))]
        lab01 = dict(zip(first[0], first[1]))
        others = [s for s in schemes if dict(zip(s[0], s[1])) != lab01]
        while True:
            per_model = [first]
            for p in poses[1:]:
                if p["near"] or rng.random() < 0.55:
                    per_model.append(others[int(rng.integers(0, len(others)))])
                else:
                    per_model.append(first)
            if len({s for s in per_model}) >= 3 and sum(s == first for s in per_model) >= 5:
                break
    else:
        per_model = [((0, 1), ("A", "B"))] * N

    # ---- write the models, then measure everything on the written coordinates
    mdir = workdir / "models"
    mdir.mkdir(parents=True, exist_ok=True)
    for f in mdir.glob("*.pdb"):
        f.unlink()
    xyz0 = [np.vstack([p["res"][i]["xyz"] for i in p["model_idx"]]) for p in prots]
    res_m = [[p["res"][i] for i in p["model_idx"]] for p in prots]
    nums = [[p["resnum"][r["pos"]] for r in res_m[k]] for k, p in enumerate(prots)]
    bfac = [np.array([r["b"] for r in res_m[k]]) for k in (0, 1)]
    pl_base = [96.0 - 34.0 * (np.argsort(np.argsort(b)) / max(len(b) - 1, 1)) ** 2 for b in bfac]
    ev = [x for x in links if lys[x["a"]][2] == "model" and lys[x["b"]][2] == "model"]
    truth_models = {}
    files = {}
    min_margin = 9.0
    for name, pose, (order2, labels) in zip(names, poses, per_model):
        g = _rot(_unit(rng), rng.uniform(0, np.pi))
        coords = [xyz0[0], xyz0[1] @ pose["r"].T + pose["t"]]
        centre = np.vstack(coords).mean(0)
        shift = rng.uniform(-25, 25, size=3)
        coords = [(c - centre) @ g.T + shift for c in coords]
        lab = dict(zip(order2, labels))
        text = _pdb_text([(lab[k], res_m[k], coords[k], nums[k],
                           np.clip(pl_base[k] + rng.normal(0, 1.2, size=len(nums[k])), 30, 98.5))
                          for k in order2])
        (mdir / f"{name}.pdb").write_text(text)
        files[f"models/{name}.pdb"] = mdir / f"{name}.pdb"
        ca = _read_ca(text)
        by_chain = [np.array([v for (c, _), v in ca.items() if c == lab[k]]) for k in (0, 1)]
        dmat = np.linalg.norm(by_chain[0][:, None] - by_chain[1][None], axis=2)
        cnt = {"inter_satisfied": 0, "inter_total": 0, "intra_satisfied": 0, "intra_total": 0}
        false_sat = 0
        for x in ev:
            (ka, pa, _, _), (kb, pb, _, _) = lys[x["a"]], lys[x["b"]]
            d = float(np.linalg.norm(ca[(lab[ka], prots[ka]["resnum"][pa])] - ca[(lab[kb], prots[kb]["resnum"][pb])]))
            min_margin = min(min_margin, abs(d - cutoff))
            cls = "inter" if ka != kb else "intra"
            cnt[cls + "_total"] += 1
            cnt[cls + "_satisfied"] += int(d <= cutoff)
            false_sat += int(d <= cutoff and x["kind"] == "false_inter")
        truth_models[name] = {**cnt, "near_native": bool(pose["near"]),
                              **{k: pose[k] for k in ("irmsd", "fnat", "lrmsd", "dockq")},
                              "false_inter_satisfied": false_sat, "ca_contacts_8A": int((dmat < 8.0).sum()), "ca_clashes_4A": int((dmat < 4.0).sum()),
                              "chains": {prots[k]["acc"]: lab[k] for k in (0, 1)}}

    # ---- release conditions measured on the written files
    near_names = [n for n in names if truth_models[n]["near_native"]]
    wrong_names = [n for n in names if not truth_models[n]["near_native"]]
    inter_counts = [truth_models[n]["inter_satisfied"] for n in names]
    gap = min(truth_models[n]["inter_satisfied"] for n in near_names) - \
        max(truth_models[n]["inter_satisfied"] for n in wrong_names)
    distinct = len(set(inter_counts))
    biggest_tie = max(inter_counts.count(v) for v in set(inter_counts))
    if min_margin < 0.5:
        raise RuntimeError(f"a distance is only {min_margin:.3f} A from the cutoff")
    if gap < COUNT_GAP:
        raise RuntimeError(f"near-native models lead by only {gap} satisfied inter-links")
    if distinct < 6 or biggest_tie > 0.4 * N:
        raise RuntimeError(f"degenerate satisfied counts: {sorted(inter_counts)}")
    if len({truth_models[n]["intra_satisfied"] for n in names}) != 1:
        raise RuntimeError("intra-protein counts must be identical in rigid-body models")

    # ---- cross-link table, FASTA, mapping truth
    rows = []
    for x in links:
        a, b = (x["a"], x["b"]) if rng.random() < 0.5 else (x["b"], x["a"])
        false = x["kind"].startswith("false")
        n_csm = 1 + int(rng.poisson(0.7 if false else 3.5))
        score = float(np.clip(rng.normal(78, 16) if false else rng.normal(112, 24), 41, 240))
        rows.append((prots[lys[a][0]]["acc"], lys[a][1], prots[lys[b][0]]["acc"], lys[b][1], n_csm, score, x["kind"]))
    rows = [rows[i] for i in rng.permutation(len(rows))]
    xl_csv = "xl_id,protein_a,position_a,protein_b,position_b,n_csm,score\n" + "".join(
        f"XL{i + 1:03d},{r[0]},{r[1]},{r[2]},{r[3]},{r[4]},{r[5]:.1f}\n" for i, r in enumerate(rows))
    fasta_order = [0, 1] if rng.random() < 0.5 else [1, 0]
    fasta = "".join(f">{prots[k]['acc']} blinded protein, full length, {len(prots[k]['full'])} aa\n"
                    + "\n".join(prots[k]["full"][i:i + 60] for i in range(0, len(prots[k]["full"]), 60)) + "\n"
                    for k in fasta_order)
    lab01 = dict(zip(per_model[0][0], per_model[0][1]))
    positions = sorted({(lys[i][0], lys[i][1]) for x in links for i in (x["a"], x["b"])})
    mapping = [{"accession": prots[k]["acc"], "position": pos,
                "chain": lab01[k] if pos in prots[k]["resnum"] else None,
                "resnum": prots[k]["resnum"].get(pos)} for k, pos in positions]
    n_unres_pos = sum(1 for m in mapping if m["chain"] is None)
    if n_unres_pos < 2:
        raise RuntimeError("fewer than two cross-linked positions are missing from the models")
    truth = {"cutoff": cutoff, "models": truth_models, "mapping": mapping,
             "crosslinks": [{"xl_id": f"XL{i + 1:03d}", "kind": r[6]} for i, r in enumerate(rows)],
             "source": {"pdb": entry["pdb"], "chains": list(entry["chains"])},
             "measured": {"min_margin_A": round(min_margin, 3), "count_gap": gap, "distinct_inter_counts": distinct,
                          "fdr": round(fdr, 3), "name_order_auc": round(float(auc), 3),
                          "unresolved_positions": n_unres_pos, "n_positions": len(mapping)}}
    ranking_json = json.dumps({"description": "ranking confidence = 0.8*ipTM + 0.2*pTM",
                               "order": names,
                               "confidence": {n: round(float(c), 3) for n, c in zip(names, conf)}}, indent=1)

    # ---- truth-derived probes / attacks
    def restraints_tsv(tm, order=None, extra=False):
        hdr = "model\tinter_satisfied\tinter_total\tintra_satisfied\tintra_total" + ("\tinter_fraction" if extra else "")
        out = [hdr]
        for n in (order or names):
            m = tm[n]
            line = f"{n}\t{m['inter_satisfied']}\t{m['inter_total']}\t{m['intra_satisfied']}\t{m['intra_total']}"
            if extra:
                line += f"\t{m['inter_satisfied'] / max(m['inter_total'], 1):.3f}"
            out.append(line)
        return "\n".join(out) + "\n"

    def selection_json(ranking):
        return json.dumps({"best_model": ranking[0], "ranking": ranking}, indent=1)

    def mapping_tsv(entries, na="NA"):
        return "accession\tposition\tmodel_chain\tmodel_resnum\n" + "".join(
            f"{m['accession']}\t{m['position']}\t{m['chain'] if m['chain'] else na}\t"
            f"{m['resnum'] if m['resnum'] is not None else na}\n" for m in entries)

    by_support = sorted(names, key=lambda n: (-truth_models[n]["inter_satisfied"], -truth_models[n]["dockq"]))
    prng = np.random.default_rng(seed + 31)
    near_sorted = sorted(near_names, key=lambda n: truth_models[n]["inter_satisfied"])   # weakest first
    wrong_shuf = [wrong_names[i] for i in prng.permutation(len(wrong_names))]
    tol_rank = near_sorted[:-1] + [wrong_shuf[0], near_sorted[-1]] + wrong_shuf[1:]
    shuffled = [names[i] for i in prng.permutation(N)]
    one_off = json.loads(json.dumps(truth_models))
    one_off[wrong_names[len(wrong_names) // 2]]["inter_satisfied"] += 1
    top_wrong = max(wrong_names, key=lambda n: truth_models[n]["inter_satisfied"])
    n_rows_inter = sum(1 for r in rows if r[0] != r[2])
    constant = {n: {"inter_satisfied": 0, "inter_total": n_rows_inter, "intra_satisfied": 0,
                    "intra_total": len(rows) - n_rows_inter} for n in names}
    acc_first = prots[fasta_order[0]]["acc"]
    identity_map = [{"accession": m["accession"], "position": m["position"],
                     "chain": "A" if m["accession"] == acc_first else "B", "resnum": m["position"]} for m in mapping]
    truth_map = mapping_tsv(mapping)
    pipe = {"xl_pipeline.py": ASSETS / "xl_pipeline.py"}

    probes = [
        Probe("planted_truth", "submission written directly from the planted truth", 0.999,
              files={"restraints.tsv": restraints_tsv(truth_models), "selection.json": selection_json(by_support),
                     "mapping.tsv": truth_map}),
        Probe("within_tolerance", "rows shuffled, extra column, lower-case NA, the weakest near-native model chosen "
              "as best, one wrong model ranked above one near-native model, wrong models in random order", 0.99,
              files={"restraints.tsv": restraints_tsv(truth_models, shuffled, extra=True),
                     "selection.json": selection_json(tol_rank),
                     "mapping.tsv": mapping_tsv([mapping[i] for i in prng.permutation(len(mapping))], na="na")}),
        Probe("independent_solution", "second implementation: no Biopython, own PDB reader, mapping by exact "
              "matching blocks (difflib), ranking by satisfied fraction with distance tie-break", 0.99,
              script="SOLUTION=2 bash /solution/solve.sh"),
    ]
    attacks = [
        Attack("prior_only_constant", "no analysis: zero satisfied links, totals = rows of the table, best = "
               "model_01, ranking = file order, positions used as residue numbers on chains A/B", 0.1,
               files={"restraints.tsv": restraints_tsv(constant), "selection.json": selection_json(names),
                      "mapping.tsv": mapping_tsv(identity_map)}),
        Attack("one_count_off", "everything correct except one satisfied count of one model (+1)", 0.89,
               files={"restraints.tsv": restraints_tsv(one_off), "selection.json": selection_json(by_support),
                      "mapping.tsv": truth_map}),
        Attack("best_is_best_wrong_model", "correct counts and mapping, but the best-scoring wrong model is "
               "put first", 0.70,
               files={"restraints.tsv": restraints_tsv(truth_models),
                      "selection.json": selection_json([top_wrong] + [n for n in by_support if n != top_wrong]),
                      "mapping.tsv": truth_map}),
        Attack("naive_pooled_counts", "inter- and intra-protein cross-links not separated: pooled numbers "
               "reported as inter, intra left at zero", 0.75,
               script="python3 /alegen_script/xl_pipeline.py --pooled", script_files=pipe),
    ]
    if level >= 2:
        attacks += [
            Attack("naive_positions_as_resnums", "cross-link positions used directly as model residue numbers",
                   0.75, script="python3 /alegen_script/xl_pipeline.py --numbering raw", script_files=pipe),
            Attack("naive_single_offset", "one constant offset per protein (start of the construct), the "
                   "trimmed internal loop ignored", 0.75,
                   script="python3 /alegen_script/xl_pipeline.py --numbering offset", script_files=pipe),
            Attack("naive_constant_chains", "chain letters of model_01 assumed for every model", 0.75,
                   script="python3 /alegen_script/xl_pipeline.py --chains first", script_files=pipe),
        ]
    if level >= 3:
        attacks += [
            Attack("naive_top_confidence", "correct counts and mapping, but selection and ranking follow the "
                   "pipeline confidence", 0.75,
                   script="python3 /alegen_script/xl_pipeline.py --select confidence", script_files=pipe),
            Attack("naive_prune_unsatisfied", "cross-links that no model satisfies are discarded as false "
                   "positives before counting", 0.75,
                   script="python3 /alegen_script/xl_pipeline.py --prune", script_files=pipe),
            Attack("naive_default_pipeline", "positions as residue numbers, chain letters of model_01, "
                   "selection by confidence", 0.60,
                   script="python3 /alegen_script/xl_pipeline.py --numbering raw --chains first --select confidence",
                   script_files=pipe),
        ]

    traps = []
    if level >= 2:
        g = [sorted(p["deleted"]) for p in prots]
        traps += [
            Trap("numbering_offset_and_gap",
                 "Cross-link positions are full-length positions; the models are numbered 1..n over trimmed "
                 f"constructs (N-terminal regions of {prots[0]['res'][0]['pos'] - 1} and "
                 f"{prots[1]['res'][0]['pos'] - 1} residues missing, internal loops of {len(g[0])} and {len(g[1])} "
                 "residues missing), so the shift changes along each chain.",
                 "the residue at a cross-linked position is not a lysine in the model (or does not exist); "
                 "chain sequences are shorter than the FASTA sequences; xlms_notes.md says the prediction "
                 "input were trimmed constructs", "Understanding / domain knowledge gap"),
            Trap("chain_identity_varies",
                 "Which chain letter (and which position in the file) is which protein differs between models.",
                 "chain lengths/sequences per chain letter differ between model files; xlms_notes.md says the "
                 "models were pooled from several runs", "Execution / unchecked assumption"),
        ]
    if level >= 3:
        traps += [
            Trap("misleading_confidence",
                 "The model with the highest pipeline confidence (model_01) is wrong; the confidence order "
                 f"separates near-native from wrong models with AUC {auc:.2f}.",
                 "the satisfied inter-protein counts disagree with ranking.json", "Approach / wrong strategy"),
            Trap("false_positive_crosslinks",
                 f"{fdr:.0%} of the reported residue pairs are false; the best model violates several of them, "
                 "and intra-protein false positives are violated by every model.",
                 "xlms_notes.md states the residue-pair FDR and that there was no manual validation; no model "
                 "satisfies every cross-link", "Understanding / domain knowledge gap"),
        ]

    inputs = {**files, "crosslinks.csv": xl_csv, "sequences.fasta": fasta,
              "xlms_notes.md": _notes(level, linker, [prots[k]["acc"] for k in fasta_order],
                                      int(round(fdr * 100)))}
    if level >= 3:
        inputs["ranking.json"] = ranking_json
    ranking_line = ("- `ranking.json` - the prediction pipeline's own confidence for every model\n"
                    if level >= 3 else "")
    return TaskBuild(
        task_id=f"xlms-model-selection-s{seed:04d}-l{level}", family="xlms_model_selection", seed=seed, level=level,
        title="Selecting a predicted complex model with XL-MS restraints",
        summary="Evaluate predicted models of a two-protein complex against lysine-lysine cross-links reported "
                "in full-length numbering and select the model whose inter-protein arrangement the data support.",
        subdomain="3.1 Biomolecular Structure & Design",
        software=["Python", "Biopython", "numpy/scipy/pandas"],
        instruction_md=INSTRUCTION.format(n_models=N, ranking_line=ranking_line),
        agent_dockerfile=dockerlib.agent_dockerfile(conda=["biopython=1.85"]),
        verifier_dockerfile=dockerlib.verifier_dockerfile(pip=[]),
        input_files=inputs,
        solution_files={"solve.sh": ASSETS / "solve.sh", "xl_pipeline.py": ASSETS / "xl_pipeline.py",
                        "solve_alt.py": ASSETS / "solve_alt.py"},
        tests_files={"score.py": ASSETS / "score.py", "reference/truth.json": json.dumps(truth, indent=1)},
        traps=traps, attacks=attacks, probes=probes,
        capabilities=["mapping experimental residue numbering onto model numbering",
                      "per-file sanity checks instead of assuming a constant layout",
                      "weighing experimental evidence against a tool's self-reported confidence",
                      "exact output contracts"],
        judge_explanation_md=JUDGE, pass_threshold=0.9,
        oracle_steps=6, oracle_runtime_min=0.1, cpus=2, memory_mb=4000, storage_mb=10000,
        extra_metadata={"n_models": N, "n_crosslinks": len(rows)},
    )
