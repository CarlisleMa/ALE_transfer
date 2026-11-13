"""Level 4 of the XL-MS model selection family (imported lazily by `generate` for level >= 4).

Level 3 plus two elements taken from the small-model pilot (README section 6, AUTHORING "Level 4"):

  P_surface_distance  the group's distance criterion is a surface distance (shortest path through
                      solvent on a lattice, written out in xlms_notes.md), not the straight line.
                      The generator picks the cross-links so that the straight-line count ranks a
                      wrong model first.
  P_homo_oligomer     the complex is a 2:2 assembly. A cross-link is satisfied if any assignment of
                      its two positions to chains is within range; a link between the same position
                      of one protein can only run between the two copies. "Same protein = same
                      chain" ranks a wrong model first.

The hidden native is a real pair of heterodimers from one RCSB entry. Models keep heterodimer 1 and
move heterodimer 2 as a rigid body. Surface distances are a computed quantity: the generator
computes them with a numpy bucket search (assets/l4/sasd_np.py), uses a cross-link only if its
surface distance is at least GRAY angstrom away from the cutoff in every shipped model, and the
release gate (G2) requires the reference solution, a separate implementation on scipy, to reproduce
every count and status. The generator itself does not run the scipy implementation.
"""
from __future__ import annotations

import importlib.util
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from ... import dockerlib
from ...spec import Attack, Probe, TaskBuild, Trap
from . import (ASSETS, BACKBONE, ENTRIES, MIN_ISLAND, NEAR_IRMSD, WRONG_IRMSD, _ambiguous, _fetch, _fit_rmsd,
               _load_chain, _nterm, _pdb_text, _plant_deletion, _resmin, _rmsd, _rot, _segments, _unit)

A4 = ASSETS / "l4"
# hidden natives: two heterodimers of one entry (body 1 is kept, body 2 is moved)
L4_ENTRIES = [{"pdb": "1FIN", "bodies": (("A", "B"), ("C", "D")), "oracle_min": 2.0},
              {"pdb": "1EFU", "bodies": (("A", "B"), ("C", "D")), "oracle_min": 3.0}]
# cross-linker, full name, spacer arm (A), the group's surface-distance cutoff (A)
L4_LINKERS = [("DSSO", "disuccinimidyl sulfoxide", 10.3, 32.0),
              ("BS3", "bis(sulfosuccinimidyl) suberate", 11.4, 34.0),
              ("DSBU", "disuccinimidyl dibutyric urea", 12.5, 36.0)]
R_BLOCK = 3.0     # a lattice point closer than this to a heavy atom is blocked
TETHER = 7.0      # longest straight segment between a CA atom and the first/last lattice point
GRAY = 3.0        # every scored surface distance is at least this far from the cutoff
COUNT_GAP = 3     # every near-native model satisfies >= this many more cross-links than any wrong model
CONDA = ["biopython=1.85"]
POOL_WRONG = 36   # wrong poses that get a full distance table; the models are chosen among them


def _sasd():
    spec = importlib.util.spec_from_file_location("alegen_xl4_sasd_np", A4 / "sasd_np.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _table_worker(args):
    """Surface distances between all lysine CA atoms of one structure (numpy implementation)."""
    atoms, anchors, limit = args
    s = _sasd()
    lat = s.Lattice(atoms, R_BLOCK)
    return s.pairwise(lat, anchors, TETHER, limit)


def _read_atoms(text: str):
    """All atoms of a model file as written: (xyz, [(chain, resnum, atom name)])."""
    xyz, meta = [], []
    for line in text.splitlines():
        if line.startswith("ATOM"):
            xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            meta.append((line[21], int(line[22:26]), line[12:16].strip()))
    return np.array(xyz), meta


def _has_exact_tie(xyz: np.ndarray, ca: np.ndarray) -> bool:
    """True if a lattice point is exactly R_BLOCK from an atom or exactly TETHER from a lysine CA
    (integer arithmetic on the coordinates as written, 0.001 A)."""
    for pts, rad in ((xyz, R_BLOCK), (ca, TETHER)):
        q = np.rint(pts * 1000).astype(np.int64)
        base = np.rint(q / 1000.0).astype(np.int64)
        f = q - 1000 * base
        k = int(rad) + 1
        r2 = int(round(rad * 1000)) ** 2
        for dx in range(-k, k + 1):
            for dy in range(-k, k + 1):
                for dz in range(-k, k + 1):
                    if ((f[:, 0] - 1000 * dx) ** 2 + (f[:, 1] - 1000 * dy) ** 2 + (f[:, 2] - 1000 * dz) ** 2 == r2).any():
                        return True
    return False


class Retry(RuntimeError):
    """The sampled instance does not meet a release condition; the generator draws another one."""


def generate(seed: int, level: int, workdir: Path) -> TaskBuild:
    last = None
    for attempt in range(6):
        try:
            return _generate(seed, level, workdir, attempt)
        except Retry as e:
            last = e
            print(f"  xlms l4 seed {seed}: attempt {attempt} rejected: {e}")
    raise RuntimeError(f"no level-4 instance for seed {seed}: {last}")


def _build(seed: int, attempt: int, workdir: Path) -> dict:
    """Proteins, poses, model files and the distance tables of one candidate instance."""
    rng = np.random.default_rng([74000 + seed, attempt])
    e4 = L4_ENTRIES[(seed - 1) % len(L4_ENTRIES)]
    entry = next(e for e in ENTRIES if e["pdb"] == e4["pdb"])
    linker = L4_LINKERS[int(rng.integers(0, len(L4_LINKERS)))]
    cutoff = linker[3]
    path = _fetch(entry)

    alnum = "ABCDEFGHJKLMNPQRSTUVWXYZ23456789"
    accs = []
    while len(accs) < 2:
        a = "XL" + "".join(alnum[i] for i in rng.integers(0, len(alnum), size=4))
        if a not in accs:
            accs.append(a)

    # ---- proteins: two copies each; the construct is what is resolved in both copies
    prots = []
    for k in (0, 1):
        cc = [_load_chain(path, e4["bodies"][b][k]) for b in (0, 1)]
        if cc[0]["seqres"] != cc[1]["seqres"]:
            raise RuntimeError("the two copies differ in sequence")
        common = {r["idx"] for r in cc[0]["res"]} & {r["idx"] for r in cc[1]["res"]}
        keep = {i for seg in _segments(sorted(common)) if len(seg) >= MIN_ISLAND for i in seg}
        ext = _nterm(rng, int(rng.integers(14, 95)))
        copies = []
        for c in cc:
            res = [r for r in c["res"] if r["idx"] in keep]
            for r in res:
                r["pos"] = len(ext) + r["idx"] + 1
            copies.append(res)
        prots.append({"acc": accs[k], "full": ext + cc[0]["seqres"], "res": copies[0], "copies": copies})
    ca_chain = [[np.array([r["xyz"][r["names"].index("CA")] for r in p["copies"][b]]) for b in (0, 1)] for p in prots]
    for k, p in enumerate(prots):
        others = np.vstack([ca_chain[1 - k][0], ca_chain[1 - k][1], ca_chain[k][1]])
        s, e = _plant_deletion(rng, p, others)
        p["deleted"] = set(range(s, e + 1))
        p["model_idx"] = [i for i in range(len(p["res"])) if i not in p["deleted"]]
        p["modelled_pos0"] = {p["res"][i]["pos"] - 1 for i in p["model_idx"]}
        p["ambiguous"] = _ambiguous(p["full"], p["modelled_pos0"])
        p["resnum"] = {p["res"][i]["pos"]: n + 1 for n, i in enumerate(p["model_idx"])}

    # ---- lysine positions. status "model" / "native" (trimmed loop) / "none" (no coordinates)
    lys = []
    for k, p in enumerate(prots):
        by_pos = {r["pos"]: i for i, r in enumerate(p["res"])}
        for j, aa in enumerate(p["full"]):
            if aa != "K" or j in p["ambiguous"]:
                continue
            pos = j + 1
            if pos in by_pos:
                lys.append((k, pos, "native" if by_pos[pos] in p["deleted"] else "model"))
            else:
                lys.append((k, pos, "none"))
    mod_lys = [(k, pos) for k, pos, st in lys if st == "model"]
    U = len(mod_lys)

    # ---- rigid bodies: body b = copy b of both proteins
    def body_stack(fn):
        return [np.vstack([fn(r) for k in (0, 1) for r in prots[k]["copies"][b]]) for b in (0, 1)]

    bb = body_stack(lambda r: r["xyz"][[r["names"].index(a) for a in BACKBONE]])
    heavy = [[r["xyz"] for k in (0, 1) for r in prots[k]["copies"][b]] for b in (0, 1)]
    resmin = _resmin(heavy[0], heavy[1])
    if_res = [np.where((resmin < 10.0).any(1))[0], np.where((resmin < 10.0).any(0))[0]]
    if len(if_res[0]) < 8 or len(if_res[1]) < 8:
        raise RuntimeError("the two heterodimers do not touch in this entry")
    if_bb = [np.vstack([bb[b][4 * i:4 * i + 4] for i in if_res[b]]) for b in (0, 1)]
    contacts = np.argwhere(resmin < 5.0)
    nat_if = np.vstack(if_bb)
    ca_mod = [np.vstack([np.array([prots[k]["copies"][b][i]["xyz"][prots[k]["copies"][b][i]["names"].index("CA")]
                                   for i in prots[k]["model_idx"]]) for k in (0, 1)]) for b in (0, 1)]
    c1, c2 = ca_mod[0].mean(0), ca_mod[1].mean(0)
    u0 = (c2 - c1) / np.linalg.norm(c2 - c1)
    if_centre = if_bb[1].mean(0)
    nat_contacts = int((np.linalg.norm(ca_mod[0][:, None] - ca_mod[1][None], axis=2) < 8.0).sum())

    def quality(r, t):
        mod_if = np.vstack([if_bb[0], if_bb[1] @ r.T + t])
        irmsd = _fit_rmsd(mod_if, nat_if)
        kept = sum(1 for i, j in contacts
                   if np.sqrt(((heavy[0][i][:, None, :] - (heavy[1][j] @ r.T + t)[None, :, :]) ** 2).sum(-1)).min() < 5.0)
        fnat = kept / len(contacts)
        lrmsd = _rmsd(bb[1] @ r.T + t, bb[1])            # body 1 is the fixed receptor
        dockq = (fnat + 1 / (1 + (irmsd / 1.5) ** 2) + 1 / (1 + (lrmsd / 8.5) ** 2)) / 3
        return {"irmsd": round(irmsd, 3), "fnat": round(fnat, 3), "lrmsd": round(lrmsd, 3), "dockq": round(dockq, 3)}

    def irmsd_of(r, t):
        return _fit_rmsd(np.vstack([if_bb[0], if_bb[1] @ r.T + t]), nat_if)

    # lysine CA atoms in the native frame: index a = 2 * u + copy
    def lys_ca(b):
        out = []
        for k, pos in mod_lys:
            r = next(x for x in prots[k]["copies"][b] if x["pos"] == pos)
            out.append(r["xyz"][r["names"].index("CA")])
        return np.array(out)

    kca = [lys_ca(0), lys_ca(1)]

    def euclid_between(r, t):
        """U x U matrix: smaller CA-CA distance of the two between-body assignments."""
        m = kca[1] @ r.T + t
        d01 = np.linalg.norm(kca[0][:, None] - m[None], axis=2)      # u in body 1, v in body 2
        return np.minimum(d01, d01.T)

    nat_eb = euclid_between(np.eye(3), np.zeros(3))
    close = np.triu(nat_eb <= 0.72 * cutoff)                         # pairs that are close across the bodies

    near = []
    for _ in range(900):
        r = _rot(_unit(rng), np.deg2rad(rng.uniform(1.0, 9.0)))
        t = if_centre - if_centre @ r.T + _unit(rng) * rng.uniform(0.0, 1.6)
        ir = irmsd_of(r, t)
        if 0.5 <= ir <= NEAR_IRMSD - 0.1:
            near.append((ir, r, t))
    n_models = int(rng.integers(15, 19))
    pool = []
    for tg in (0.75, 1.25, 1.75):                      # two candidates per target; one of each is kept
        for _ in range(2):
            x = min((x for x in near if not any(x[1] is p["r"] for p in pool)), key=lambda x: abs(x[0] - tg))
            qd = quality(x[1], x[2])
            if qd["irmsd"] >= NEAR_IRMSD or qd["dockq"] < 0.49:
                raise Retry("a near-native pose is not near-native")
            pool.append({"near": True, "target": tg, "r": x[1], "t": x[2], **qd})

    # wrong poses: body 2 in another orientation, pushed against body 1 until the two bodies share
    # about as many CA-CA contacts as in the native. Part of them stays near the native site
    # (rotated in place), which is where straight-line and surface distance disagree most.
    wrong = []
    ca2c = ca_mod[1] - c2
    for _ in range(1400):
        cls = rng.random()
        perp = np.cross(u0, _unit(rng))
        if cls < 0.35:
            theta, spin = np.deg2rad(rng.uniform(3, 25)), np.deg2rad(rng.uniform(20, 75))
        elif cls < 0.6:
            theta, spin = np.deg2rad(rng.uniform(15, 55)), np.deg2rad(rng.uniform(10, 180))
        else:
            theta, spin = np.deg2rad(rng.uniform(55, 180)), np.deg2rad(rng.uniform(10, 180))
        u = _rot(perp, theta) @ u0
        q = _rot(_unit(rng), spin)
        b = ca2c @ q.T
        d = (ca_mod[0] - c1)[:, None, :] - b[None, :, :]
        p_ = d @ u
        q2 = (d ** 2).sum(-1) - p_ ** 2
        m = q2 < 64.0
        reach = np.sort(p_[m] + np.sqrt(64.0 - q2[m]))[::-1]
        want = int(rng.uniform(0.6, 1.6) * nat_contacts)
        if len(reach) < want:
            continue
        t = c1 + float(reach[want - 1]) * u - c2 @ q.T
        if irmsd_of(q, t) < WRONG_IRMSD:
            continue
        eb = euclid_between(q, t)
        wrong.append((float((eb[close] <= cutoff - GRAY).mean()), q, t))
    if len(wrong) < 120:
        raise Retry("not enough wrong poses")
    wrong.sort(key=lambda x: -x[0])
    # the poses that keep most of the close pairs within straight-line range, then a random spread
    cand = wrong[:10] + [wrong[int(i)] for i in 10 + rng.permutation(len(wrong) - 10)]
    for x in cand:
        if sum(1 for p in pool if not p["near"]) >= POOL_WRONG:
            break
        qd = quality(x[1], x[2])
        if qd["dockq"] < 0.23 and qd["fnat"] < 0.3:
            pool.append({"near": False, "r": x[1], "t": x[2], "close_frac": round(x[0], 3), **qd})
    if sum(1 for p in pool if not p["near"]) < POOL_WRONG:
        raise Retry("not enough clearly wrong poses")

    # ---- every pool pose gets its own frame; coordinates are rounded as they will be written
    res_m = {(k, b): [prots[k]["copies"][b][i] for i in prots[k]["model_idx"]] for k in (0, 1) for b in (0, 1)}
    xyz0 = {kb: np.vstack([r["xyz"] for r in rs]) for kb, rs in res_m.items()}
    nums = {kb: [prots[kb[0]]["resnum"][r["pos"]] for r in rs] for kb, rs in res_m.items()}
    kbs = [(0, 0), (1, 0), (0, 1), (1, 1)]
    ca_row = {}                                              # row of the CA atom of a modelled lysine in a chain block
    for kb in kbs:
        row, k = 0, kb[0]
        for r in res_m[kb]:
            if r["aa"] == "K":
                ca_row[(kb, r["pos"])] = row + r["names"].index("CA")
            row += len(r["names"])
    jobs = []
    for pose in pool:
        for _try in range(50):
            g = _rot(_unit(rng), rng.uniform(0, np.pi))
            raw = {kb: (x if kb[1] == 0 else x @ pose["r"].T + pose["t"]) for kb, x in xyz0.items()}
            centre = np.vstack([raw[kb] for kb in kbs]).mean(0)
            shift = rng.uniform(-25, 25, size=3)
            coords = {kb: np.array([[float(f"{v:.3f}") for v in row] for row in (raw[kb] - centre) @ g.T + shift])
                      for kb in kbs}
            xyz = np.vstack([coords[kb] for kb in kbs])
            anchors = np.array([coords[(k, b)][ca_row[((k, b), pos)]] for k, pos in mod_lys for b in (0, 1)])
            if not _has_exact_tie(xyz, anchors):
                break
        else:
            raise Retry("could not place a model without an exact lattice tie")
        pose["coords"] = coords
        cab = [np.vstack([coords[(k, b)][[i for i, r in enumerate(_atom_res(res_m[(k, b)])) if r]] for k in (0, 1)])
               for b in (0, 1)]
        dmat = np.linalg.norm(cab[0][:, None] - cab[1][None], axis=2)
        pose["ca_contacts_8A"], pose["ca_clashes_4A"] = int((dmat < 8.0).sum()), int((dmat < 4.0).sum())
        jobs.append((xyz, anchors, cutoff + GRAY))
    # the native complex (constructs as modelled, native frame) decides which pairs are true cross-links
    nat_xyz = np.vstack([xyz0[kb] for kb in kbs])
    nat_anch = np.array([kca[b][u] for u in range(U) for b in (0, 1)])
    jobs.append((nat_xyz, nat_anch, cutoff + 8.5))
    with ProcessPoolExecutor(max_workers=2) as ex:
        tabs = list(ex.map(_table_worker, jobs))
    S = np.array(tabs[:-1])                                  # poses x 2U x 2U surface distances (inf beyond cutoff + GRAY)
    E = np.array([np.linalg.norm(j[1][:, None] - j[1][None], axis=2) for j in jobs[:-1]])
    return dict(rng=rng, e4=e4, linker=linker, cutoff=cutoff, prots=prots, lys=lys, mod_lys=mod_lys, pool=pool,
                S=S, S_nat=tabs[-1], E=E, res_m=res_m, nums=nums, kbs=kbs, kca=kca, n_models=n_models)


def _atom_res(residues):
    """Per atom of a chain block: True where the atom is a CA."""
    return [n == "CA" for r in residues for n in r["names"]]


def _reduce(T: np.ndarray):
    """(.., 2U, 2U) atom tables -> position tables: within one body (same copy index) and between bodies."""
    a00, a11 = T[..., 0::2, 0::2], T[..., 1::2, 1::2]
    a01, a10 = T[..., 0::2, 1::2], T[..., 1::2, 0::2]
    return np.minimum(a00, a11), np.minimum(a01, a10)


def _select(B: dict, rng) -> dict:
    """Choose the models and the cross-links from the pool tables so that every release condition holds."""
    c, pool, mod_lys, prots = B["cutoff"], B["pool"], B["mod_lys"], B["prots"]
    U = len(mod_lys)
    Sw, Sb = _reduce(B["S"])
    Ew, Eb = _reduce(B["E"])
    Nw, Nb = _reduce(B["S_nat"])
    kk = np.array([k for k, _ in mod_lys])
    pp = np.array([p for _, p in mod_lys])
    het = kk[:, None] != kk[None]
    diag = np.eye(U, dtype=bool)
    up = np.triu(np.ones((U, U), dtype=bool))
    Sm = np.where(diag, Sb, np.minimum(Sw, Sb))            # the group's rule: any assignment (same position: between copies)
    Em = np.where(diag, Eb, np.minimum(Ew, Eb))            # the same with straight-line distances
    nat = np.where(diag, Nb, np.minimum(Nw, Nb))
    sat, vio = Sm <= c - GRAY, ~(Sm < c + GRAY)
    clean = sat | vio
    esat = Em <= c
    n2 = np.where(het, sat, np.where(diag, True, Sw <= c))  # "two positions of one protein = one chain"
    true_any = up & (nat <= c - 4.0)
    within = true_any & ~diag & (Nw <= c - 4.0)
    between = true_any & (diag | ~(Nw < c + GRAY))
    far = up & ~(nat <= c + 8.0)
    near_by_t = {}
    for i, p in enumerate(pool):
        if p["near"]:
            near_by_t.setdefault(p["target"], []).append(i)
    wrong_idx = [i for i, p in enumerate(pool) if not p["near"]]
    best = None
    for a0 in near_by_t[0.75]:
        for a1 in near_by_t[1.25]:
            for a2 in near_by_t[1.75]:
                n = int((between & sat[[a0, a1, a2]].all(0)).sum())
                if best is None or n > best[0]:
                    best = (n, [a0, a1, a2])
    Nr = best[1]
    near_sat = sat[Nr].all(0)
    TBc = between & near_sat
    Fc = far & vio[Nr].all(0) & ~esat[Nr].any(0)           # false candidates: out of range in the near-native models

    def pick(mask, n, prefer=None):
        idx = np.argwhere(mask)
        idx = idx[rng.permutation(len(idx))]
        if prefer is not None:
            idx = np.array(sorted(idx.tolist(), key=lambda uv: -int(prefer[uv[0], uv[1]]))).reshape(-1, 2)
        return [tuple(int(v) for v in uv) for uv in idx[:n]]

    def as_mask(pairs):
        m = np.zeros((U, U), dtype=bool)
        for u, v in pairs:
            m[u, v] = True
        return m

    plans, why = [], []
    for h in wrong_idx:
        hv = TBc & ~het & vio[h]                            # same-protein between-copy links that pose h violates
        fh = Fc & het & sat[h]                              # false hetero pairs that pose h satisfies
        hs = TBc & het & sat[h]
        if hv.sum() >= 5 and fh.sum() >= 1:
            plans.append((min(int(hs.sum()), 8) + 0.1 * min(int(hv.sum()), 10), h))
    plans.sort(reverse=True)
    n_wrong = B["n_models"] - 3
    for _, h in plans[:10]:
        hv = TBc & ~het & vio[h]
        n_fh = 2 if (hv.sum() >= 6 and (Fc & het & sat[h]).sum() >= 2) else 1
        hv_l = pick(hv, 8, prefer=esat[h])
        tb = pick(TBc & het & sat[h], 8) + hv_l + pick(TBc & ~het & sat[h], 3)
        f_h = pick(Fc & het & sat[h], n_fh)
        FH = as_mask(f_h)

        def compatible(core, n_true):
            return [p for p in wrong_idx if p != h and clean[p][core].all()
                    and int(sat[p][core].sum()) <= n_true - COUNT_GAP]

        # drop the links that sit in the gray zone of many pool poses until enough wrong models remain
        while True:
            ok = compatible(as_mask(tb) | FH, len(tb))
            if len(ok) >= n_wrong - 1:
                break
            core = as_mask(tb) | FH
            bad = [p for p in wrong_idx if p != h and not clean[p][core].all()]
            keep_hv = sum(1 for x in tb if x in hv_l) <= n_fh + COUNT_GAP + 1
            score = {x: sum(1 for p in bad if not clean[p][x]) for x in tb if not (keep_hv and x in hv_l)}
            if not score or max(score.values()) == 0 or len(tb) <= 8:
                break
            tb.remove(max(score, key=lambda x: score[x]))
        if len(ok) < n_wrong - 1:
            why.append(f"h={h}: only {len(ok)} of {n_wrong - 1} further wrong models are compatible")
            continue
        TB = as_mask(tb)
        # straight-line look-alike: a wrong pose that keeps the chosen links within straight-line range
        options = []
        for e in [h] + ok:
            if not clean[e][TB | FH].all() or int(sat[e][TB | FH].sum()) > len(tb) - COUNT_GAP:
                continue
            rest = [p for p in ok if p != e]
            fe_c = Fc & esat[e] & vio[e] & vio[h] & ~FH
            fe_c &= np.array([clean[p] for p in rest]).sum(0) >= len(rest) - max(0, len(rest) - (n_wrong - len({h, e})))
            short = int((~esat[e][TB]).sum()) - int(esat[e][FH].sum())   # links the look-alike is behind by straight line
            need = max(0, short + 2)                                     # aim at a lead of 2, accept 1
            if need > min(3, int(fe_c.sum())):
                need = max(0, short + 1)
            if need <= min(4, int(fe_c.sum())):
                options.append((need, -pool[e]["close_frac"], e))
        if not options:
            why.append(f"h={h}: no straight-line look-alike")
            continue
        need, _, e = min(options)
        rest = [p for p in ok if p != e]
        fe_c = Fc & esat[e] & vio[e] & vio[h] & ~FH
        f_e = pick(fe_c, need, prefer=np.array([clean[p] for p in rest]).sum(0))
        FE = as_mask(f_e)
        core = TB | FH | FE
        ok = [p for p in rest if clean[p][core].all() and int(sat[p][core].sum()) <= len(tb) - COUNT_GAP]
        n_other = n_wrong - len({h, e})
        if len(ok) < n_other:
            why.append(f"h={h} e={e}: only {len(ok)} of {n_other} further wrong models are compatible")
            continue
        bins = {}
        for p in ok:
            bins.setdefault(int(sat[p][core].sum()), []).append(p)
        others = []
        while len(others) < n_other:
            for cnt in sorted(bins, reverse=True):
                if bins[cnt] and len(others) < n_other:
                    others.append(bins[cnt].pop(int(rng.integers(0, len(bins[cnt])))))
        M = Nr + sorted({h, e}, key=lambda x: (x != e, x)) + others      # pool indices of the models; e right after the near ones
        cl = clean[M].all(0) & ~core
        used = np.zeros(U, dtype=int)
        for u, v in tb + f_h + f_e:
            used[u] += 1
            used[v] += 1

        def take(mask, n, cap=4):
            got = []
            for u, v in pick(mask, 10 ** 6):
                if len(got) >= n:
                    break
                if used[u] >= cap or used[v] >= cap:
                    continue
                used[u] += 1
                used[v] += 1
                got.append((u, v))
            return got

        tw_het = take(up & het & within & cl & near_sat, int(rng.integers(10, 15)))
        tw_homo = take(up & ~het & ~diag & within & cl & near_sat & (np.abs(pp[:, None] - pp[None]) >= 2),
                       int(rng.integers(20, 27)))
        # false pairs that every model has within straight-line range but far apart along the surface
        # (the path would have to go through the protein), and one that is far by both measures
        thru = far & cl & vio[M].all(0) & esat[M].all(0) & ~diag
        f_thru = take(thru & het, 1, cap=5) + take(thru & ~het, 1, cap=5)
        f_thru += take(thru, 3 - len(f_thru), cap=5)
        f_far = take(far & cl & vio[M].all(0) & ~esat[M].any(0), 1, cap=5)
        links = ([(u, v, "true_between") for u, v in tb] + [(u, v, "true_within") for u, v in tw_het + tw_homo]
                 + [(u, v, "false_satisfied_by_a_wrong_model") for u, v in f_h]
                 + [(u, v, "false_straight_line_only_in_a_wrong_model") for u, v in f_e]
                 + [(u, v, "false_through_the_protein") for u, v in f_thru] + [(u, v, "false_far") for u, v in f_far])
        L = as_mask([(u, v) for u, v, _ in links])
        tot = lambda T: np.array([int(T[m][L].sum()) for m in M])          # noqa: E731
        truth, eu, nn = tot(sat), tot(esat), tot(n2)
        gap = int(truth[:3].min() - truth[3:].max())
        gap_e = int(eu[3:].max() - eu[:3].max())
        gap_n = int(nn[3:].max() - nn[:3].max())
        distinct = len(set(truth.tolist()))
        tie = max(truth.tolist().count(v) for v in set(truth.tolist()))
        if gap < COUNT_GAP or gap_e < 1 or gap_n < 1 or distinct < 5 or tie > 0.6 * len(M) \
                or len(tw_het) < 6 or len(tw_homo) < 12 or len(f_thru) < 2:
            why.append(f"h={h} e={e}: gap {gap}, straight-line lead {gap_e}, one-chain lead {gap_n}, distinct "
                       f"{distinct}, tie {tie}, within links {len(tw_het)}+{len(tw_homo)}")
            continue
        return dict(M=M, links=links, sat=sat, esat=esat, n2=n2, Sm=Sm, Em=Em, het=het, diag=diag, h=h, e=e,
                    measured={"count_gap": gap, "straight_line_lead_of_a_wrong_model": gap_e,
                              "one_chain_rule_lead_of_a_wrong_model": gap_n, "distinct_totals": distinct,
                              "truth_totals": truth.tolist(), "straight_line_totals": eu.tolist(),
                              "one_chain_rule_totals": nn.tolist(),
                              "n_between": len(tb), "n_between_hetero": int((as_mask(tb) & het).sum()),
                              "n_between_same_position": int((as_mask(tb) & diag).sum()),
                              "n_within": len(tw_het) + len(tw_homo),
                              "n_false_satisfied_by_a_wrong_model": len(f_h),
                              "n_false_straight_line_only_in_a_wrong_model": len(f_e),
                              "n_false_through_the_protein": len(f_thru), "n_false_far": len(f_far)})
    raise Retry("no combination of models and cross-links meets the release conditions: "
                + (" | ".join(why) if why else "no wrong pose violates >= 5 same-protein between-copy links while "
                                               "satisfying a false inter-protein pair"))


# ----------------------------------------------------------------------------- instruction / notes

INSTRUCTION = """\
You are the integrative structural biologist of a protein-complex group. A structure-prediction run
produced a set of candidate models of a protein complex that contains two copies of each of two
proteins, and the mass-spectrometry facility delivered lysine-lysine cross-links (XL-MS) measured on
the same complex in solution. The group needs every model evaluated against the cross-links under
the group's own evaluation rules and wants to know which model the data support best.

Task directory:
- `/task`

Inputs (read-only) under `/task/input`:
- `models/model_NN.pdb` - {n_models} predicted models of the complex (`model_01` ... `model_{n_models:02d}`),
  PDB format, four protein chains each (two copies of each of the two proteins)
- `ranking.json` - the prediction pipeline's own confidence for every model
- `crosslinks.csv` - the identified cross-links, one row per unique cross-linked residue pair:
  `xl_id`, `protein_a`, `position_a`, `protein_b`, `position_b`, `n_csm`, `score`. Proteins are
  given by accession, positions are residue positions in the full-length protein sequence. The
  table does not say which copy of a protein a residue belongs to.
- `sequences.fasta` - full-length sequences of the two proteins (the accession is the first word of
  the header). Protein identities were blinded.
- `xlms_notes.md` - hand-over notes of the mass-spectrometry facility and the modelling group,
  including the group's rules for deciding whether a structure satisfies a cross-link

Available software (on `PATH`): Python 3.12 with numpy, scipy, pandas and Biopython. There is no
internet access. No program for the group's distance measure is installed; `xlms_notes.md` defines
the measure completely.

Your task:
1. Work out which residues of the models each cross-linked position corresponds to.
2. For every model and every cross-link, decide whether the model satisfies the cross-link under
   the rules of `xlms_notes.md`, and count the satisfied cross-links per model, separately for
   inter-protein and intra-protein cross-links.
3. Rank the models by the number of cross-links they satisfy and select the best-supported model.

Required outputs under `/task/output`:
- `restraints.tsv` - tab-separated with at least the columns `model`, `inter_satisfied`,
  `inter_total`, `intra_satisfied`, `intra_total`; one row per model; `model` is the file name
  without `.pdb`; the four counts are integers. Definitions:
  - Every row of `crosslinks.csv` is one cross-link. It is inter-protein if `protein_a` and
    `protein_b` are different proteins and intra-protein if they are the same protein (no matter
    whether its two residues lie in one copy or in the two copies of that protein).
  - A cross-link is evaluable in a model if both cross-linked positions are present in that model.
    `inter_total` / `intra_total` are the numbers of evaluable cross-links of each class. A
    cross-link with a position that is not present in the model is counted neither as satisfied nor
    in the total.
  - An evaluable cross-link is satisfied if the model meets the rules that `xlms_notes.md` gives
    in the section "How the group evaluates a cross-link against a structure"; otherwise it is
    violated.
- `link_status.tsv` - tab-separated with the columns `model`, `xl_id`, `status`; one row for every
  combination of a model and a cross-link; `status` is `satisfied`, `violated` or `not_evaluable`
  as defined above.
- `selection.json` - `{{"best_model": "model_NN", "ranking": ["model_NN", ...]}}`. `ranking` lists
  every model exactly once, ordered by the total number of satisfied cross-links (inter-protein
  plus intra-protein), largest number first; models with the same number may be given in any
  order. Its first element is `best_model`.
- `mapping.tsv` - tab-separated with the columns `accession`, `position`, `model_chain`,
  `model_resnum`; one row for every distinct cross-linked position (accession + position) that
  occurs in `crosslinks.csv`. `model_chain` lists the chain identifiers of all chains of
  `models/model_01.pdb` that are copies of that protein, separated by commas; `model_resnum` is the
  residue number of that residue in these chains. Write `NA` in both columns if the residue is not
  present in the model.

Constraints:
- Do not modify anything under `/task/input`; write deliverables only under `/task/output`.
- Every model and every cross-link must be reported.
- Do not ask for confirmation. Execute directly.
"""

JUDGE = """\
**Where the truth comes from.** The hidden native is a pair of heterodimers from one RCSB entry.
Every shipped model keeps heterodimer 1 and moves heterodimer 2 as a rigid body, so interface RMSD,
fnat and DockQ of each model are known from the transform. The status of a cross-link in a model
follows the definition in `xlms_notes.md`: surface distance on a 1 A lattice (blocking radius,
tether and cutoff as stated there), minimum over the allowed assignments of the two positions to the
copies. This is a computed quantity. The generator therefore (i) uses a cross-link only if its
surface distance is at least 3 A away from the cutoff in every shipped model, (ii) rejects a model
frame in which a lattice point is exactly at the blocking radius of an atom or exactly at the tether
length from a cross-linked CA atom, and (iii) the release gate requires that two independent
implementations (reference: Biopython + k-d tree + sparse-graph Dijkstra in scipy; second: own
reader + bucketed lattice search in numpy) reproduce every count and every status. True cross-links
are lysine pairs within range in the native complex, false ones are far apart in the native.

**Hard gates (score 0).** Missing/unparseable file; `restraints.tsv` not listing every model exactly
once with integer counts; `link_status.tsv` not holding exactly one row per model and cross-link
with a valid status; `selection.json` without a `best_model` that is a model, without a `ranking`
that is a permutation of all models, or with `ranking[0] != best_model`; duplicated positions in
`mapping.tsv`.

**Components.**

| weight | component | metric | full credit | zero |
|---|---|---|---|---|
| 0.30 | restraint_counts | fraction of models whose four counts are all exact | all models exact (otherwise 0.6 x fraction) | no model exact |
| 0.15 | link_status | three-class Matthews correlation over all model x cross-link statuses | >= 0.995 | <= 0.95 |
| 0.30 | best_model | hidden interface RMSD of `best_model` | near-native (< 2 A) | any other model |
| 0.10 | ranking | ROC AUC of the submitted order for near-native versus wrong models | >= 0.95 | <= 0.70 |
| 0.15 | residue_mapping | fraction of cross-linked positions with the exact set of chains + residue number in model_01 (or NA) | every position exact (otherwise 0.6 x fraction) | no position exact |

**Why these tolerances.** Counts and statuses are exact because every scored surface distance is
at least 3 A from the cutoff under both implementations (measured margins are in the design notes).
One wrong model costs 0.30 x (1 - 0.6 x (N-1)/N) > 0.12 and fails the task. The status table gives
graded credit for the same information; straight-line distances and the one-chain reading of the
copies both land below its zero point. Any near-native model is accepted as `best_model`; the
ranking is scored only on near-native versus wrong pairs. Every near-native model satisfies at least
3 more cross-links than any wrong model; a wrong model leads under straight-line distances and
under "two positions of one protein lie in one chain" (planted, measured by the generator).

**What is deliberately not scored.** The order among wrong models and among near-native models; any
distance table, plot or report.
"""


def _notes(linker: tuple, accs: list[str], fdr_pct: int) -> str:
    name, longname, spacer, cutoff = linker
    return "\n".join([
        "# XL-MS hand-over notes", "",
        "## Sample and cross-linking",
        f"- Sample: complex of {accs[0]} and {accs[1]} (identities blinded for this hand-over), both "
        "full-length, co-expressed, purified by affinity and size-exclusion chromatography. SEC-MALS and "
        "native MS give a 2:2 assembly: two copies of each protein.",
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
        "- A row with the same protein and the same position on both sides comes from two linked peptides "
        "that cover the same lysine.",
        "- Error control: target-decoy, list cut at 5% FDR on the spectrum-match level, which the "
        f"facility estimates to correspond to about {fdr_pct}% FDR on the residue-pair level. Decoy "
        "hits have been removed from the table; there was no manual validation.",
        "",
        "## How the group evaluates a cross-link against a structure",
        "- Distance measure. The cross-linker cannot pass through protein, so the group does not use the "
        "straight-line distance between the two lysines. It uses the surface distance defined below, "
        "measured between the CA atoms of the two linked lysines on the complete structure (all chains "
        f"present). For {name} a pair of residues is within range if its surface distance is at most "
        f"{cutoff:.1f} A.",
        "- Surface distance between two CA atoms P and Q of a structure:",
        "  1. Lattice. The lattice points are the points whose x, y and z coordinates are all whole numbers "
        "of angstrom, in the coordinate frame of the structure file as it is. The lattice is not bounded.",
        "  2. Blocked points. A lattice point is blocked if its distance to the centre of at least one "
        f"non-hydrogen atom of the model is less than {R_BLOCK:.1f} A (every atom of every chain counts, the "
        "two linked lysines included). All other lattice points are open.",
        "  3. Paths. A path is a sequence of open lattice points in which consecutive points differ by at "
        "most 1 in each of x, y and z (a point has up to 26 neighbours). The length of a path is the sum of "
        "the straight-line distances between consecutive points (1, sqrt(2) or sqrt(3) A per step). A path "
        "may consist of a single point and then has length 0.",
        "  4. Tethers. P is joined to the first point of the path and Q to the last point of the path by "
        f"straight segments. Each of the two segments may be at most {TETHER:.1f} A long and is not tested "
        "against the protein (it stands for the lysine side chain).",
        "  5. The surface distance is the smallest value of (segment at P) + (path length) + (segment at Q) "
        "over all choices of an open first point, an open last point and a path between them. If there is "
        "no such choice, the surface distance is infinite and the pair is out of range.",
        "- Copies. The mass spectrometer cannot tell the two copies of a protein apart. A cross-link "
        "between position i of one protein and position j of a protein (the same or the other one) is "
        "satisfied by a structure if at least one assignment of the two positions to chains of the "
        "structure gives a pair of residues that is within range. The two ends of a cross-link are always "
        "two different residues, so for a row with the same protein and the same position on both sides "
        "the only assignment is that position in the two different copies.",
        "- A cross-link to a position that is absent from a structure cannot be evaluated for that "
        "structure; the group leaves it out of both the satisfied count and the total.",
        "- The group ranks structures by the number of cross-links they satisfy.",
        "",
        "## Models",
        "- The models come from the group's structure-prediction pipeline. The prediction input were the "
        "expression constructs used for crystallisation trials, which are shorter than the full-length "
        "proteins (disordered terminal regions and flexible internal loops had been trimmed). Both copies "
        "of a protein were predicted from the same construct.",
        "- The model files are the unmodified output of the pipeline. They were pooled from several runs "
        "that were set up by different people.",
        "- `ranking.json` holds the pipeline's ranking confidence (0.8 x ipTM + 0.2 x pTM); the model files "
        "are named in the order of that confidence.",
    ]) + "\n"


# ----------------------------------------------------------------------------- generator

def _generate(seed: int, level: int, workdir: Path, attempt: int) -> TaskBuild:
    B = _build(seed, attempt, workdir)
    rng, c, pool, prots, mod_lys, lys = B["rng"], B["cutoff"], B["pool"], B["prots"], B["mod_lys"], B["lys"]
    sel = _select(B, rng)
    M, links = sel["M"], sel["links"]
    U, N = len(mod_lys), len(M)
    kbs = B["kbs"]

    # ---- cross-links to positions that are missing from the models (never evaluable)
    def native_ca(k, pos):
        out = []
        for b in (0, 1):
            r = next((x for x in prots[k]["copies"][b] if x["pos"] == pos), None)
            if r is not None:
                out.append(r["xyz"][r["names"].index("CA")])
        return out

    unres = []
    loop_lys = [(k, pos) for k, pos, st in lys if st == "native"]
    none_lys = [(k, pos) for k, pos, st in lys if st == "none"]
    cands = []
    for k2, pos2 in loop_lys:
        for u, (k1, pos1) in enumerate(mod_lys):
            d = min(float(np.linalg.norm(x - y)) for x in native_ca(k1, pos1) for y in native_ca(k2, pos2))
            if d <= c - 12.0:
                cands.append(((k1, pos1), (k2, pos2)))
    for i in rng.permutation(len(cands))[:3]:
        unres.append((*cands[int(i)], "true_unresolved"))
    if none_lys:
        k2, pos2 = none_lys[int(rng.integers(0, len(none_lys)))]
        unres.append((mod_lys[int(rng.integers(0, U))], (k2, pos2), "false_unresolved"))
    if len({x[1] for x in unres}) < 2:
        raise Retry("fewer than two cross-linked positions are missing from the models")

    # ---- file order = confidence order. model_01 is the straight-line look-alike; the order does not
    # favour the near-native models
    while True:
        order = [int(i) for i in rng.permutation(N)]
        order.remove(3)
        order.insert(0, 3)
        rank = {m: i for i, m in enumerate(order)}
        auc = float(np.mean([rank[a] < rank[b] for a in range(3) for b in range(3, N)]))
        if 0.30 <= auc <= 0.60:
            break
    idx = [M[i] for i in order]                               # pool index of model_01, model_02, ...
    near_flag = [pool[i]["near"] for i in idx]
    names = [f"model_{i + 1:02d}" for i in range(N)]
    conf = np.sort(rng.uniform(0.36, 0.89, size=N))[::-1]
    conf[0] = max(conf[0], conf[1] + 0.02)

    # ---- chain layout per model. chain = (protein, copy); order in the file and letters vary
    orders = [((0, 0), (1, 0), (0, 1), (1, 1)), ((0, 0), (0, 1), (1, 0), (1, 1)),
              ((1, 1), (0, 1), (1, 0), (0, 0)), ((1, 0), (1, 1), (0, 1), (0, 0)),
              ((0, 1), (1, 0), (1, 1), (0, 0)), ((1, 0), (0, 0), (1, 1), (0, 1))]
    letters = [("A", "B", "C", "D"), ("B", "C", "D", "E"), ("A", "C", "B", "D")]
    schemes = [(o, l) for o in orders for l in letters]
    het_links = [(u, v) for u, v, _ in links if mod_lys[u][0] != mod_lys[v][0]]

    def l3_inter(i, scheme):
        """Inter-protein count of the level-3 method: one chain per protein (the last one in the file),
        straight-line CA-CA distance."""
        last = {kb[0]: kb[1] for kb in scheme[0]}
        E = B["E"][idx[i]]
        return sum(1 for u, v in het_links if E[2 * u + last[mod_lys[u][0]], 2 * v + last[mod_lys[v][0]]] <= c)

    l3_lead = None
    for _try in range(400):
        first = schemes[int(rng.integers(0, len(schemes)))]
        per_model = [first]
        for i in range(1, N):
            if near_flag[i] or rng.random() < 0.6:
                per_model.append(schemes[int(rng.integers(0, len(schemes)))])
            else:
                per_model.append(first)
        if len(set(per_model)) < 4 or sum(s == first for s in per_model) < 3 or \
                any(s == first for s, nf in zip(per_model, near_flag) if nf):
            continue
        cnt = [l3_inter(i, s) for i, s in enumerate(per_model)]
        l3_lead = max(x for x, nf in zip(cnt, near_flag) if not nf) - max(x for x, nf in zip(cnt, near_flag) if nf)
        if l3_lead >= 1:
            break
    else:
        raise Retry("no chain layout under which the level-3 method ranks a wrong model first")

    # ---- write the models
    mdir = workdir / "models"
    mdir.mkdir(parents=True, exist_ok=True)
    for f in mdir.glob("*.pdb"):
        f.unlink()
    res_m, nums = B["res_m"], B["nums"]
    pl = {}
    for kb, rs in res_m.items():
        bf = np.array([r["b"] for r in rs])
        pl[kb] = 96.0 - 34.0 * (np.argsort(np.argsort(bf)) / max(len(bf) - 1, 1)) ** 2
    files, labels = {}, []
    for name, pi, (order4, labs) in zip(names, idx, per_model):
        lab = dict(zip(order4, labs))
        labels.append(lab)
        coords = pool[pi]["coords"]
        text = _pdb_text([(lab[kb], res_m[kb], coords[kb], nums[kb],
                           np.clip(pl[kb] + rng.normal(0, 1.2, size=len(nums[kb])), 30, 98.5)) for kb in order4])
        (mdir / f"{name}.pdb").write_text(text)
        files[f"models/{name}.pdb"] = mdir / f"{name}.pdb"
        # the written coordinates are the ones the tables were computed on
        xyz, meta = _read_atoms(text)
        back = np.vstack([xyz[[j for j, m in enumerate(meta) if m[0] == lab[kb]]] for kb in kbs])
        if not np.array_equal(back, np.vstack([coords[kb] for kb in kbs])):
            raise RuntimeError("written coordinates differ from the tabulated ones")

    # ---- cross-link table and the truth per model
    sat, het, Sm = sel["sat"], sel["het"], sel["Sm"]
    rows = []
    for u, v, kind in links:
        rows.append((mod_lys[u], mod_lys[v], kind, (u, v)))
    for a_, b_, kind in unres:
        rows.append((a_, b_, kind, None))
    out_rows = []
    for a_, b_, kind, uv in rows:
        if rng.random() < 0.5:
            a_, b_ = b_, a_
        false = kind.startswith("false")
        n_csm = 1 + int(rng.poisson(0.7 if false else 3.5))
        score = float(np.clip(rng.normal(78, 16) if false else rng.normal(112, 24), 41, 240))
        out_rows.append((prots[a_[0]]["acc"], a_[1], prots[b_[0]]["acc"], b_[1], n_csm, score, kind, uv))
    out_rows = [out_rows[i] for i in rng.permutation(len(out_rows))]
    xl_ids = [f"XL{i + 1:03d}" for i in range(len(out_rows))]
    xl_csv = "xl_id,protein_a,position_a,protein_b,position_b,n_csm,score\n" + "".join(
        f"{x},{r[0]},{r[1]},{r[2]},{r[3]},{r[4]},{r[5]:.1f}\n" for x, r in zip(xl_ids, out_rows))
    fdr = sum(1 for r in out_rows if r[6].startswith("false")) / len(out_rows)
    truth_models, status = {}, {}
    margins_sat, n_beyond = [], 0
    for name, pi, lab in zip(names, idx, labels):
        cnt = {"inter_satisfied": 0, "inter_total": 0, "intra_satisfied": 0, "intra_total": 0}
        st = {}
        for x, r in zip(xl_ids, out_rows):
            if r[7] is None:
                st[x] = "not_evaluable"
                continue
            u, v = r[7]
            cls = "inter" if r[0] != r[2] else "intra"
            ok = bool(sat[pi][u, v])
            cnt[cls + "_total"] += 1
            cnt[cls + "_satisfied"] += int(ok)
            st[x] = "satisfied" if ok else "violated"
            if ok:
                margins_sat.append(c - float(Sm[pi][u, v]))
            else:
                n_beyond += 1
        p = pool[pi]
        status[name] = st
        truth_models[name] = {**cnt, "near_native": bool(p["near"]),
                              **{q: p[q] for q in ("irmsd", "fnat", "lrmsd", "dockq")},
                              "ca_contacts_8A": p["ca_contacts_8A"], "ca_clashes_4A": p["ca_clashes_4A"],
                              "chains": {prots[k]["acc"]: sorted(lab[(k, b)] for b in (0, 1)) for k in (0, 1)}}
    near_names = [n for n in names if truth_models[n]["near_native"]]
    wrong_names = [n for n in names if not truth_models[n]["near_native"]]
    total = {n: truth_models[n]["inter_satisfied"] + truth_models[n]["intra_satisfied"] for n in names}
    gap = min(total[n] for n in near_names) - max(total[n] for n in wrong_names)
    if gap < COUNT_GAP:
        raise RuntimeError("count gap lost between selection and files")

    fasta_order = [0, 1] if rng.random() < 0.5 else [1, 0]
    fasta = "".join(f">{prots[k]['acc']} blinded protein, full length, {len(prots[k]['full'])} aa\n"
                    + "\n".join(prots[k]["full"][i:i + 60] for i in range(0, len(prots[k]["full"]), 60)) + "\n"
                    for k in fasta_order)
    lab01 = labels[0]
    acc_k = {p["acc"]: k for k, p in enumerate(prots)}
    positions = sorted({(acc_k[r[a]], r[a + 1]) for r in out_rows for a in (0, 2)})
    mapping = [{"accession": prots[k]["acc"], "position": pos,
                "chains": sorted(lab01[(k, b)] for b in (0, 1)) if pos in prots[k]["resnum"] else None,
                "resnum": prots[k]["resnum"].get(pos)} for k, pos in positions]
    n_unres_pos = sum(1 for m in mapping if m["chains"] is None)
    measured = {**sel["measured"], "count_gap_on_files": gap, "fdr": round(fdr, 3),
                "name_order_auc": round(auc, 3), "level3_method_lead_of_a_wrong_model": int(l3_lead),
                "unresolved_positions": n_unres_pos, "n_positions": len(mapping),
                "min_margin_satisfied_A": round(min(margins_sat), 3), "n_satisfied_calls": len(margins_sat),
                "n_violated_calls_all_at_least_gray_beyond": n_beyond, "attempt": attempt,
                "pool_index_of_models": idx}
    truth = {"cutoff": c, "block_radius": R_BLOCK, "tether": TETHER, "gray": GRAY, "models": truth_models,
             "status": status, "mapping": mapping,
             "crosslinks": [{"xl_id": x, "kind": r[6]} for x, r in zip(xl_ids, out_rows)],
             "source": {"pdb": B["e4"]["pdb"], "bodies": [list(b) for b in B["e4"]["bodies"]]},
             "measured": measured}
    ranking_json = json.dumps({"description": "ranking confidence = 0.8*ipTM + 0.2*pTM", "order": names,
                               "confidence": {n: round(float(v), 3) for n, v in zip(names, conf)}}, indent=1)

    # ---- truth-derived probes / attacks
    def restraints_tsv(tm, order=None, extra=False):
        out = ["model\tinter_satisfied\tinter_total\tintra_satisfied\tintra_total" + ("\ttotal_satisfied" if extra else "")]
        for n in (order or names):
            m = tm[n]
            line = f"{n}\t{m['inter_satisfied']}\t{m['inter_total']}\t{m['intra_satisfied']}\t{m['intra_total']}"
            out.append(line + (f"\t{m['inter_satisfied'] + m['intra_satisfied']}" if extra else ""))
        return "\n".join(out) + "\n"

    def status_tsv(st, order=None, fmt=str):
        return "model\txl_id\tstatus\n" + "".join(f"{n}\t{x}\t{fmt(st[n][x])}\n" for n in (order or names) for x in xl_ids)

    def selection_json(ranking):
        return json.dumps({"best_model": ranking[0], "ranking": ranking}, indent=1)

    def mapping_tsv(entries, na="NA", sep=","):
        return "accession\tposition\tmodel_chain\tmodel_resnum\n" + "".join(
            f"{m['accession']}\t{m['position']}\t{sep.join(m['chains']) if m['chains'] else na}\t"
            f"{m['resnum'] if m['resnum'] is not None else na}\n" for m in entries)

    by_support = sorted(names, key=lambda n: (-total[n], -truth_models[n]["dockq"]))
    prng = np.random.default_rng(seed + 431)
    near_sorted = sorted(near_names, key=lambda n: total[n])
    wrong_shuf = [wrong_names[i] for i in prng.permutation(len(wrong_names))]
    tol_rank = near_sorted[:-1] + [wrong_shuf[0], near_sorted[-1]] + wrong_shuf[1:]
    shuffled = [names[i] for i in prng.permutation(N)]
    one_off = json.loads(json.dumps(truth_models))
    one_off[wrong_names[len(wrong_names) // 2]]["inter_satisfied"] += 1
    top_wrong = max(wrong_names, key=lambda n: total[n])
    n_rows_inter = sum(1 for r in out_rows if r[0] != r[2])
    constant = {n: {"inter_satisfied": 0, "inter_total": n_rows_inter, "intra_satisfied": 0,
                    "intra_total": len(out_rows) - n_rows_inter} for n in names}
    acc_first = prots[fasta_order[0]]["acc"]
    identity_map = [{"accession": m["accession"], "position": m["position"],
                     "chains": ["A" if m["accession"] == acc_first else "B"], "resnum": m["position"]} for m in mapping]
    all_violated = {n: {x: "violated" for x in xl_ids} for n in names}
    truth_map = mapping_tsv(mapping)
    alt = {"solve_alt.py": A4 / "solve_alt.py", "sasd_np.py": A4 / "sasd_np.py"}
    run = "python3 /alegen_script/solve_alt.py "

    probes = [
        Probe("planted_truth", "submission written directly from the generator's truth tables", 0.999,
              files={"restraints.tsv": restraints_tsv(truth_models), "link_status.tsv": status_tsv(status),
                     "selection.json": selection_json(by_support), "mapping.tsv": truth_map}),
        Probe("within_tolerance", "rows shuffled, extra column, upper-case statuses, lower-case NA, chains "
              "separated by ' / ', the weakest near-native model chosen as best, one wrong model ranked above "
              "one near-native model, wrong models in random order", 0.99,
              files={"restraints.tsv": restraints_tsv(truth_models, shuffled, extra=True),
                     "link_status.tsv": status_tsv(status, shuffled, fmt=str.upper),
                     "selection.json": selection_json(tol_rank),
                     "mapping.tsv": mapping_tsv([mapping[i] for i in prng.permutation(len(mapping))], na="na", sep=" / ")}),
        Probe("independent_solution", "second implementation: no Biopython, no scipy; own PDB reader, mapping by "
              "exact matching blocks (difflib), surface distance by a bucketed lattice search in numpy", 0.99,
              script="SOLUTION=2 bash /solution/solve.sh"),
    ]
    attacks = [
        Attack("prior_only_constant", "no analysis: nothing satisfied, totals = rows of the table, best = "
               "model_01, ranking = file order, positions used as residue numbers on chains A/B", 0.1,
               files={"restraints.tsv": restraints_tsv(constant), "link_status.tsv": status_tsv(all_violated),
                      "selection.json": selection_json(names), "mapping.tsv": mapping_tsv(identity_map)}),
        Attack("one_count_off", "everything correct except one satisfied count of one model (+1)", 0.89,
               files={"restraints.tsv": restraints_tsv(one_off), "link_status.tsv": status_tsv(status),
                      "selection.json": selection_json(by_support), "mapping.tsv": truth_map}),
        Attack("best_is_best_wrong_model", "correct counts, statuses and mapping, but the best-scoring wrong "
               "model is put first", 0.72,
               files={"restraints.tsv": restraints_tsv(truth_models), "link_status.tsv": status_tsv(status),
                      "selection.json": selection_json([top_wrong] + [n for n in by_support if n != top_wrong]),
                      "mapping.tsv": truth_map}),
        Attack("l3_reference", "the level-3 reference method on the level-4 inputs: numbering and chain identity "
               "handled per model, but straight-line CA-CA distance, one chain per protein, ranking by the "
               "inter-protein count", 0.75,
               script=run + "--distance straight --copies last --select inter", script_files=alt),
        Attack("naive_straight_line", "ignores P_surface_distance only: copies handled as defined, straight-line "
               "CA-CA distance with the cutoff of the notes", 0.75,
               script=run + "--distance straight", script_files=alt),
        Attack("naive_same_protein_same_chain", "ignores P_homo_oligomer only: surface distance as defined, but "
               "two positions of one protein are looked up in one chain (a same-position link becomes a "
               "residue with itself)", 0.75, script=run + "--copies within", script_files=alt),
        Attack("naive_one_chain_per_protein", "surface distance as defined, but every protein is represented by "
               "one chain of the model", 0.75, script=run + "--copies last", script_files=alt),
        Attack("naive_pooled_counts", "inter- and intra-protein cross-links not separated: pooled numbers "
               "reported as inter, intra left at zero", 0.75, script=run + "--pooled", script_files=alt),
        Attack("naive_positions_as_resnums", "cross-link positions used directly as model residue numbers",
               0.75, script=run + "--numbering raw", script_files=alt),
        Attack("naive_single_offset", "one constant offset per protein (start of the construct), the trimmed "
               "internal loop ignored", 0.75, script=run + "--numbering offset", script_files=alt),
        Attack("naive_constant_chains", "chain letters of model_01 assumed for every model", 0.75,
               script=run + "--chains first", script_files=alt),
        Attack("naive_top_confidence", "correct counts, statuses and mapping, but selection and ranking follow "
               "the pipeline confidence", 0.75, script=run + "--select confidence", script_files=alt),
        Attack("naive_prune_unsatisfied", "cross-links that no model satisfies are discarded as false "
               "positives before counting", 0.75, script=run + "--prune", script_files=alt),
        Attack("naive_default_pipeline", "positions as residue numbers, chain letters of model_01, one chain "
               "per protein, straight-line distance, selection by confidence", 0.60,
               script=run + "--numbering raw --chains first --copies last --distance straight --select confidence",
               script_files=alt),
    ]

    g = [sorted(p["deleted"]) for p in prots]
    m = sel["measured"]
    traps = [
        Trap("numbering_offset_and_gap",
             "Cross-link positions are full-length positions; the models are numbered 1..n over trimmed "
             f"constructs (N-terminal regions of {prots[0]['res'][0]['pos'] - 1} and "
             f"{prots[1]['res'][0]['pos'] - 1} residues missing, internal loops of {len(g[0])} and {len(g[1])} "
             "residues missing), so the shift changes along each chain.",
             "the residue at a cross-linked position is not a lysine in the model (or does not exist); "
             "chain sequences are shorter than the FASTA sequences; xlms_notes.md says the prediction "
             "input were trimmed constructs", "Understanding / domain knowledge gap"),
        Trap("chain_identity_varies",
             "Which chain letters (and which positions in the file) are which protein differs between models.",
             "chain lengths/sequences per chain letter differ between model files; xlms_notes.md says the "
             "models were pooled from several runs", "Execution / unchecked assumption"),
        Trap("misleading_confidence",
             "The model with the highest pipeline confidence (model_01) is wrong; the confidence order "
             f"separates near-native from wrong models with AUC {auc:.2f}.",
             "the satisfied counts disagree with ranking.json", "Approach / wrong strategy"),
        Trap("false_positive_crosslinks",
             f"{fdr:.0%} of the reported residue pairs are false; the best model violates them.",
             "xlms_notes.md states the residue-pair FDR and that there was no manual validation; no model "
             "satisfies every cross-link", "Understanding / domain knowledge gap"),
        Trap("P_surface_distance",
             "The group's criterion is a surface distance (shortest path through solvent on a lattice), not the "
             "straight line. With straight-line distances and the same cutoff a wrong model satisfies "
             f"{m['straight_line_lead_of_a_wrong_model']} more cross-links than the best near-native model, and "
             "model_01 is such a model.",
             "xlms_notes.md, section 'How the group evaluates a cross-link against a structure', defines the "
             "measure step by step (lattice, blocked points, paths, tethers, cutoff)", "Approach / wrong strategy"),
        Trap("P_homo_oligomer",
             "The complex holds two copies of each protein. A cross-link is satisfied if any assignment of its "
             "two positions to chains is within range; a row with the same position on both sides can only run "
             "between the two copies. Looking two positions of one protein up in one chain puts a wrong model "
             f"{m['one_chain_rule_lead_of_a_wrong_model']} cross-link(s) ahead of the near-native models.",
             "every model has four chains with two sequences; the instruction and xlms_notes.md ('Copies') state "
             "the rule; crosslinks.csv contains rows with identical ends", "Understanding / domain knowledge gap"),
    ]
    inputs = {**files, "crosslinks.csv": xl_csv, "sequences.fasta": fasta, "ranking.json": ranking_json,
              "xlms_notes.md": _notes(B["linker"], [prots[k]["acc"] for k in fasta_order], int(round(fdr * 100)))}
    return TaskBuild(
        task_id=f"xlms-model-selection-s{seed:04d}-l{level}", family="xlms_model_selection", seed=seed, level=level,
        title="Selecting a predicted 2:2 complex model with XL-MS restraints and a surface-distance criterion",
        summary="Evaluate predicted models of a complex with two copies of each of two proteins against "
                "lysine-lysine cross-links under the group's surface-distance rule and select the model the data "
                "support.",
        subdomain="3.1 Biomolecular Structure & Design",
        software=["Python", "Biopython", "numpy/scipy/pandas"],
        instruction_md=INSTRUCTION.format(n_models=N),
        agent_dockerfile=dockerlib.agent_dockerfile(conda=CONDA),
        verifier_dockerfile=dockerlib.verifier_dockerfile(pip=[]),
        input_files=inputs,
        solution_files={"solve.sh": A4 / "solve.sh", "xl4_pipeline.py": A4 / "xl4_pipeline.py",
                        "solve_alt.py": A4 / "solve_alt.py", "sasd_np.py": A4 / "sasd_np.py"},
        tests_files={"score.py": A4 / "score.py", "reference/truth.json": json.dumps(truth, indent=1)},
        traps=traps, attacks=attacks, probes=probes,
        capabilities=["implementing a written distance definition exactly (shortest path on a lattice)",
                      "handling copies of a protein in a complex",
                      "mapping experimental residue numbering onto model numbering",
                      "per-file sanity checks instead of assuming a constant layout",
                      "weighing experimental evidence against a tool's self-reported confidence",
                      "exact output contracts"],
        judge_explanation_md=JUDGE, pass_threshold=0.9,
        oracle_steps=13, oracle_runtime_min=B["e4"]["oracle_min"], cpus=2, memory_mb=8000, storage_mb=10000,
        extra_metadata={"n_models": N, "n_crosslinks": len(out_rows)},
    )
