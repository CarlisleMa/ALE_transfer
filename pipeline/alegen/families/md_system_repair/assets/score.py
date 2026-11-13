"""Verifier for md_system_repair: a defective crystal structure -> simulation-ready Amber system.

Everything is recomputed from the submitted topology and coordinates (ParmEd for the structure,
OpenMM for the physics). report.json is only compared with those recomputed values.
The hidden reference (reference/truth.json) is the complete structure the generator started from
plus the list of defects it planted; no tool output is part of the reference.
"""
import difflib
import math
import os
import sys

sys.path.insert(0, "/tests")
from _scoring import CandidateError, ramp, run  # noqa: E402

import numpy as np  # noqa: E402

AA = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())
ALIAS = {"CYX": "CYS", "CYM": "CYS", "HID": "HIS", "HIE": "HIS", "HIP": "HIS", "HSD": "HIS", "HSE": "HIS",
         "HSP": "HIS", "ASH": "ASP", "GLH": "GLU", "LYN": "LYS"}
MAX_ATOMS = 150_000
WEIGHTS = {"sequence": 0.07, "atoms_and_charge": 0.07, "disulfides": 0.08, "termini_and_breaks": 0.08,
           "conformers": 0.06, "solvent_and_ions": 0.10, "physical_sanity": 0.07, "equilibrated": 0.08,
           "report": 0.09, "production_ready": 0.30}

# tolerances (justified with measured numbers in docs/families/md_system_repair.md)
GATE_CA_RMSD = 1.0          # A, as-built CA atoms vs the hidden structure
ATOM_DEV = 1.0              # A, resolved heavy atom counted as "in place"
CHI1_TOL = 40.0             # degrees
SALT_TOL = 0.02             # mol/L
PAD_LOW, PAD_HIGH = 0.5, 5.5  # A below / above the SOP padding
XTAL_RADIUS = 0.5           # A
XTAL_MAX_FRACTION = 0.2
DENSITY_BAND = (0.97, 1.08)   # g/cm3
TEMP_BAND = (250.0, 350.0)    # K, kinetic temperature of the stored velocities
XH_ALONG_BOND_MAX = 5.0       # A/ps, below this the stored velocities come from constrained dynamics
MD_TEMP_BAND = (280.0, 320.0)  # K, after 0.3 ps of verifier MD started from the restart
EQ_CA_RMSD = 3.0            # A
BOX_REPORT_TOL = 0.02       # A
KB = 0.0083144626           # kJ/mol/K


def canon(name):
    n = name.strip().upper()
    if n in AA:
        return n
    if n in ALIAS:
        return ALIAS[n]
    if len(n) == 4 and n[0] in "NC" and (n[1:] in AA or n[1:] in ALIAS):
        return canon(n[1:])
    return None


def kabsch(mobile, target):
    """Rotation R and translation t with mobile @ R.T + t ~ target; returns (R, t, rmsd)."""
    mc, tc = mobile.mean(axis=0), target.mean(axis=0)
    h = (mobile - mc).T @ (target - tc)
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    t = tc - mc @ rot.T
    diff = mobile @ rot.T + t - target
    return rot, t, float(np.sqrt((diff ** 2).sum(axis=1).mean()))


def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    return math.degrees(math.atan2(np.dot(np.cross(b1n, v), w), np.dot(v, w)))


def load_coords(path, what, natom):
    import parmed
    try:
        rst = parmed.amber.Rst7(str(path))
        xyz = np.asarray(rst.coordinates, dtype=float).reshape(-1, 3)
        box = None if rst.box is None else np.asarray(rst.box, dtype=float).reshape(-1)[:6]
        vel = None
        if getattr(rst, "hasvels", False):
            vel = np.asarray(rst.vels, dtype=float).reshape(-1, 3)
    except CandidateError:
        raise
    except Exception as e:
        raise CandidateError(f"{what} is not a readable Amber coordinate/restart file ({type(e).__name__})")
    if xyz.shape[0] != natom:
        raise CandidateError(f"{what} holds {xyz.shape[0]} atoms but system.prmtop has {natom}")
    if not np.all(np.isfinite(xyz)):
        raise CandidateError(f"{what} contains non-finite coordinates")
    return xyz, box, vel


def main(ctx):
    import parmed

    truth = ctx.ref_json("truth.json")
    sop = truth["sop"]
    tres = truth["residues"]
    checks = {}     # component -> {name: value in [0,1]}

    def put(comp, name, value, note=None):
        checks.setdefault(comp, {})[name] = float(max(0.0, min(1.0, value)))
        if note is not None:
            ctx.info[f"{comp}.{name}"] = note

    # ------------------------------------------------------------------ load
    p_prm = ctx.path("system.prmtop", 80_000_000)
    p_crd = ctx.path("system.inpcrd", 80_000_000)
    p_rst = ctx.path("equilibrated.rst7", 80_000_000)
    report = ctx.json("report.json", 200_000)
    head = p_prm.read_bytes()[:4000]
    ctx.gate(b"%FLAG" in head and b"%VERSION" in head, "system.prmtop is an Amber parm7 topology")
    try:
        parm = parmed.amber.AmberParm(str(p_prm))
    except Exception as e:
        raise CandidateError(f"system.prmtop cannot be parsed as an Amber topology ({type(e).__name__})")
    natom = len(parm.atoms)
    ctx.gate(0 < natom <= MAX_ATOMS, f"topology has between 1 and {MAX_ATOMS} atoms")
    xyz0, box0, _ = load_coords(p_crd, "system.inpcrd", natom)
    ctx.gate(box0 is not None and len(box0) >= 6 and np.all(np.isfinite(box0)) and np.all(box0[:3] > 1.0),
             "system.inpcrd carries periodic box dimensions")

    Z = np.array([a.atomic_number for a in parm.atoms])
    mass = np.array([a.mass for a in parm.atoms])
    if np.all(Z <= 0):                                   # topology without ATOMIC_NUMBER: infer from mass
        Z = np.array([1 if m < 3.5 else 6 if m < 13 else 7 if m < 15 else 8 if m < 17 else 11 if m < 24
                      else 16 if m < 33 else 17 if m < 36 else 99 for m in mass])
    charge = np.array([a.charge for a in parm.atoms])
    nbr = [[] for _ in range(natom)]
    for b in parm.bonds:
        nbr[b.atom1.idx].append(b.atom2.idx)
        nbr[b.atom2.idx].append(b.atom1.idx)
    res_of = np.array([a.residue.idx for a in parm.atoms])

    # ------------------------------------------------------------------ classify residues
    waters, ions, solute = [], [], []
    for r in parm.residues:
        idx = [a.idx for a in r.atoms]
        zs = [Z[i] for i in idx]
        heavy = [z for z in zs if z > 1]
        if len(idx) == 1:
            ions.append(r)
        elif heavy == [8] and zs.count(1) == 2 and len(idx) <= 5 and canon(r.name) is None:
            waters.append(r)
        else:
            solute.append(r)
    aa_res = [r for r in solute if canon(r.name) is not None]
    other_res = [r for r in solute if canon(r.name) is None]
    ctx.gate(len(aa_res) >= 5, "topology contains a polypeptide")

    # ------------------------------------------------------------------ map onto the hidden structure
    tnames = [t["name"] for t in tres]
    snames = [canon(r.name) for r in aa_res]
    sm = difflib.SequenceMatcher(None, tnames, snames, autojunk=False)
    pairs = []                                            # (truth index, residue)
    for blk in sm.get_matching_blocks():
        for k in range(blk.size):
            pairs.append((blk.a + k, aa_res[blk.b + k]))
    ctx.gate(len(pairs) >= 0.7 * len(tres), "at least 70% of the expected residues are present in sequence order")

    def atom(res, name):
        for a in res.atoms:
            if a.name == name:
                return a.idx
        return None

    ca_pairs = [(ti, atom(r, "CA")) for ti, r in pairs]
    ca_pairs = [(ti, i) for ti, i in ca_pairs if i is not None]
    ctx.gate(len(ca_pairs) >= 0.7 * len(tres), "matched residues carry CA atoms")
    t_ca = np.array([tres[ti]["atoms"]["CA"] for ti, _ in ca_pairs])
    rot, trans, rmsd0 = kabsch(xyz0[[i for _, i in ca_pairs]], t_ca)
    ctx.info["as_built_ca_rmsd"] = round(rmsd0, 4)
    ctx.gate(rmsd0 <= GATE_CA_RMSD, "as-built CA atoms superpose on the structure that was handed over")
    fit0 = xyz0 @ rot.T + trans                           # as-built coordinates in the input frame
    s_of = {ti: r for ti, r in pairs}                     # truth index -> submitted residue
    t_of = {r.idx: ti for ti, r in pairs}                 # submitted residue index -> truth index

    # ------------------------------------------------------------------ sequence
    ident = len(pairs) / max(len(tres), len(aa_res))
    put("sequence", "residue_identity", ramp(ident, 0.90, 1.0), round(ident, 4))
    put("sequence", "no_extra_or_missing_residues", 1.0 if ident == 1.0 else 0.0)
    put("sequence", "no_foreign_residues", 0.0 if other_res_foreign(other_res, Z) else 1.0,
        sorted({r.name for r in other_res})[:10])

    # ------------------------------------------------------------------ heavy atoms, net charge
    complete = 0
    bad_atoms = []
    for ti, t in enumerate(tres):
        r = s_of.get(ti)
        if r is None:
            continue
        have = sorted(a.name for a in r.atoms if Z[a.idx] > 1)
        if have == sorted(t["atoms"]):
            complete += 1
        elif len(bad_atoms) < 8:
            bad_atoms.append(t["label"])
    frac_complete = complete / len(tres)
    put("atoms_and_charge", "heavy_atoms_complete", ramp(frac_complete, 0.90, 1.0), round(frac_complete, 4))
    if bad_atoms:
        ctx.info["atoms_and_charge.residues_with_wrong_heavy_atoms"] = bad_atoms
    q_solute = float(sum(charge[a.idx] for r in solute for a in r.atoms))
    put("atoms_and_charge", "protein_net_charge", 1.0 if abs(q_solute - truth["net_charge"]) <= 0.05 else 0.0,
        round(q_solute, 4))

    # ------------------------------------------------------------------ disulfides
    found, false_ss = set(), 0
    for b in parm.bonds:
        i, j = b.atom1.idx, b.atom2.idx
        if Z[i] == 16 and Z[j] == 16 and res_of[i] != res_of[j]:
            ti, tj = t_of.get(res_of[i]), t_of.get(res_of[j])
            has_h = any(Z[k] == 1 for k in nbr[i] + nbr[j])
            if ti is None or tj is None or has_h:
                false_ss += 1
            else:
                found.add(tuple(sorted((ti, tj))))
    want = {tuple(sorted(p)) for p in truth["disulfides"]}
    tp = len(found & want)
    jac = tp / max(1, len(found | want) + false_ss)
    put("disulfides", "ss_bonds", jac, f"{tp}/{len(want)} true, {len(found - want) + false_ss} false")

    # ------------------------------------------------------------------ termini, chain breaks, connectivity
    def outside(i):
        return [k for k in nbr[i] if res_of[k] != res_of[i]]

    frag_of = {}
    for fi, (lo, hi) in enumerate(truth["fragments"]):
        for ti in range(lo, hi + 1):
            frag_of[ti] = fi
    first, last = s_of.get(0), s_of.get(len(tres) - 1)
    ok = 0.0
    if first is not None and atom(first, "N") is not None:
        n = atom(first, "N")
        nh = sum(1 for k in nbr[n] if Z[k] == 1)
        ok = 1.0 if (not outside(n) and nh == (2 if tnames[0] == "PRO" else 3)) else 0.0
    put("termini_and_breaks", "n_terminus_charged", ok)
    ok = 0.0
    if last is not None and atom(last, "C") is not None:
        c = atom(last, "C")
        oxy = [k for k in nbr[c] if Z[k] == 8]
        ok = 1.0 if (not outside(c) and len(oxy) == 2 and not any(Z[m] == 1 for k in oxy for m in nbr[k])) else 0.0
    put("termini_and_breaks", "c_terminus_charged", ok)
    n_pep, n_pep_ok, unexpected = 0, 0, 0
    for ti in range(len(tres) - 1):
        if frag_of[ti] != frag_of[ti + 1]:
            continue
        n_pep += 1
        ra, rb = s_of.get(ti), s_of.get(ti + 1)
        if ra is None or rb is None:
            continue
        c, n = atom(ra, "C"), atom(rb, "N")
        if c is not None and n is not None and n in nbr[c]:
            n_pep_ok += 1
    put("termini_and_breaks", "peptide_bonds", ramp(n_pep_ok / max(1, n_pep), 0.9, 1.0),
        f"{n_pep_ok}/{n_pep}")
    for b in parm.bonds:
        i, j = b.atom1.idx, b.atom2.idx
        ti, tj = t_of.get(res_of[i]), t_of.get(res_of[j])
        if res_of[i] == res_of[j] or ti is None or tj is None:
            continue
        if Z[i] == 16 and Z[j] == 16:
            continue
        lo, hi = (i, j) if ti < tj else (j, i)
        is_pep = (abs(ti - tj) == 1 and frag_of[ti] == frag_of[tj]
                  and parm.atoms[lo].name == "C" and parm.atoms[hi].name == "N")
        if not is_pep:
            unexpected += 1
    put("termini_and_breaks", "no_bond_between_unlinked_residues", 1.0 if unexpected == 0 else 0.0, unexpected)
    aa_idx = {r.idx for r in aa_res}

    def cap_ok(res_idx, want_comp, link_z):
        if res_idx in aa_idx:
            return False
        zs = sorted(Z[a.idx] for a in parm.residues[res_idx].atoms if Z[a.idx] > 1)
        return zs == want_comp and link_z
    for fi in range(len(truth["fragments"]) - 1):
        end, start = truth["fragments"][fi][1], truth["fragments"][fi + 1][0]
        ok_c = ok_n = 0.0
        ra, rb = s_of.get(end), s_of.get(start)
        if ra is not None and atom(ra, "C") is not None:
            c = atom(ra, "C")
            out = outside(c)
            oxy = [k for k in nbr[c] if Z[k] == 8]
            if len(out) == 1 and Z[out[0]] == 7 and len(oxy) == 1 and \
                    cap_ok(res_of[out[0]], [6, 7], sum(1 for k in nbr[out[0]] if Z[k] == 1) == 1):
                ok_c = 1.0
        if rb is not None and atom(rb, "N") is not None:
            n = atom(rb, "N")
            out = outside(n)
            if len(out) == 1 and Z[out[0]] == 6 and \
                    cap_ok(res_of[out[0]], [6, 6, 8], any(Z[k] == 8 for k in nbr[out[0]])):
                ok_n = 1.0
        put("termini_and_breaks", f"break{fi + 1}_c_side_n_methylamide", ok_c)
        put("termini_and_breaks", f"break{fi + 1}_n_side_acetyl", ok_n)

    # ------------------------------------------------------------------ resolved atoms / conformers
    n_res_atoms, n_in_place = 0, 0
    for ti, t in enumerate(tres):
        r = s_of.get(ti)
        for name in t["resolved"]:
            n_res_atoms += 1
            i = atom(r, name) if r is not None else None
            if i is not None and np.linalg.norm(fit0[i] - np.array(t["atoms"][name])) <= ATOM_DEV:
                n_in_place += 1
    frac_place = n_in_place / max(1, n_res_atoms)
    put("conformers", "resolved_atoms_in_place", ramp(frac_place, 0.90, 0.99), round(frac_place, 4))
    for w in truth["chi1_watch"]:
        r = s_of.get(w["idx"])
        ok = 0.0
        if r is not None:
            ids = [atom(r, x) for x in ("N", "CA", "CB", w["gamma"])]
            if all(i is not None for i in ids):
                chi = dihedral(*[xyz0[i] for i in ids])
                d = abs((chi - w["chi1"] + 180.0) % 360.0 - 180.0)
                ok = 1.0 if d <= CHI1_TOL else 0.0
                ctx.info[f"conformers.chi1_delta_{tres[w['idx']]['label']}"] = round(d, 1)
        put("conformers", f"major_conformer_{tres[w['idx']]['label']}", ok)

    # ------------------------------------------------------------------ solvent and ions
    q_tot = float(charge.sum())
    put("solvent_and_ions", "neutral", 1.0 if abs(q_tot) <= 0.05 else 0.0, round(q_tot, 4))
    n_na = sum(1 for r in ions if Z[r.atoms[0].idx] == 11 and abs(charge[r.atoms[0].idx] - 1.0) < 0.01)
    n_cl = sum(1 for r in ions if Z[r.atoms[0].idx] == 17 and abs(charge[r.atoms[0].idx] + 1.0) < 0.01)
    n_wat = len(waters)
    put("solvent_and_ions", "only_na_cl_ions", 1.0 if (n_na + n_cl == len(ions) and n_wat > 0) else 0.0,
        f"na={n_na} cl={n_cl} other_ions={len(ions) - n_na - n_cl}")
    conc = (min(n_na, n_cl) / n_wat * 55.5) if n_wat else 0.0
    excess_ok = abs(n_cl - n_na) == abs(int(round(truth["net_charge"])))
    put("solvent_and_ions", "salt_concentration",
        1.0 if (abs(conc - sop["salt_molar"]) <= SALT_TOL and excess_ok) else 0.0, round(conc, 4))
    tip3p = n_wat > 0
    for r in waters[:: max(1, n_wat // 500)]:                # sample up to ~500 waters
        idx = [a.idx for a in r.atoms]
        if len(idx) != 3:
            tip3p = False
            break
        for i in idx:
            want_q = -0.834 if Z[i] == 8 else 0.417
            if abs(charge[i] - want_q) > 0.002:
                tip3p = False
    put("solvent_and_ions", "water_model", 1.0 if tip3p else 0.0)
    sol_idx = [a.idx for r in solute for a in r.atoms]
    ext = xyz0[sol_idx].max(axis=0) - xyz0[sol_idx].min(axis=0)
    pads = (box0[:3] - ext) / 2.0
    ortho = bool(np.all(np.abs(box0[3:6] - 90.0) < 0.01))
    pad_ok = ortho and pads.min() >= sop["padding"] - PAD_LOW and pads.max() <= sop["padding"] + PAD_HIGH
    put("solvent_and_ions", "box_padding", 1.0 if pad_ok else 0.0, [round(float(p), 2) for p in pads])
    sites = np.array(truth["xtal_sites"]) if truth["xtal_sites"] else np.zeros((0, 3))
    kept_frac = 0.0
    if len(sites):
        probe = [a.idx for r in waters for a in r.atoms if Z[a.idx] == 8] + [r.atoms[0].idx for r in ions]
        pts = fit0[probe]
        hit = 0
        for s in sites:
            if len(pts) and np.min(np.abs(pts - s).max(axis=1)) <= XTAL_RADIUS and \
                    np.min(np.linalg.norm(pts - s, axis=1)) <= XTAL_RADIUS:
                hit += 1
        kept_frac = hit / len(sites)
    put("solvent_and_ions", "no_crystal_solvent_or_ions", 1.0 if kept_frac <= XTAL_MAX_FRACTION else 0.0,
        round(kept_frac, 4))

    # ------------------------------------------------------------------ physics (OpenMM, CPU)
    import openmm as mm
    import openmm.app as app
    import openmm.unit as u

    heavy_bonds = [(b.atom1.idx, b.atom2.idx) for b in parm.bonds
                   if Z[b.atom1.idx] > 1 and Z[b.atom2.idx] > 1 and res_of[b.atom1.idx] in {r.idx for r in solute}]
    max_bond = max((float(np.linalg.norm(xyz0[i] - xyz0[j])) for i, j in heavy_bonds), default=0.0)
    put("physical_sanity", "as_built_bond_lengths", 1.0 if 0.0 < max_bond <= 2.6 else 0.0, round(max_bond, 3))

    threads = str(min(2, os.cpu_count() or 1))
    sim_ok = True
    e0 = e_min = float("nan")
    system = None
    try:
        prm = app.AmberPrmtopFile(str(p_prm))
        system = prm.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=0.8 * u.nanometer,
                                  constraints=app.HBonds, rigidWater=True)
    except Exception as e:
        sim_ok = False
        ctx.info["physical_sanity.openmm_error"] = f"{type(e).__name__}: {str(e)[:200]}"
    put("physical_sanity", "openmm_builds_periodic_system", 1.0 if sim_ok else 0.0)

    def context(xyz, box, seed):
        integ = mm.LangevinMiddleIntegrator(300 * u.kelvin, 1.0 / u.picosecond, 0.002 * u.picoseconds)
        integ.setRandomNumberSeed(seed)
        plat = mm.Platform.getPlatformByName("CPU")
        ctxt = mm.Context(system, integ, plat, {"Threads": threads})
        ctxt.setPeriodicBoxVectors(mm.Vec3(box[0] / 10, 0, 0), mm.Vec3(0, box[1] / 10, 0),
                                   mm.Vec3(0, 0, box[2] / 10))
        ctxt.setPositions((xyz / 10.0).tolist())
        return ctxt, integ

    def pot(ctxt):
        return ctxt.getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)

    if sim_ok:
        try:
            ctxt, _ = context(xyz0, box0, 1)
            e0 = pot(ctxt)
            mm.LocalEnergyMinimizer.minimize(ctxt, 10.0, 60)
            e_min = pot(ctxt)
            del ctxt
        except Exception as e:
            ctx.info["physical_sanity.minimise_error"] = f"{type(e).__name__}: {str(e)[:200]}"
    ctx.info["physical_sanity.energy_as_built_kj"] = None if not math.isfinite(e0) else round(e0, 0)
    ctx.info["physical_sanity.energy_minimised_kj"] = None if not math.isfinite(e_min) else round(e_min, 0)
    put("physical_sanity", "minimises_to_negative_energy",
        1.0 if (math.isfinite(e_min) and e_min < 0.0) else 0.0)

    # ------------------------------------------------------------------ equilibrated restart
    eq = {"readable_with_box": 0.0, "density": 0.0, "velocities_thermal": 0.0, "protein_intact": 0.0,
          "md_continues_stably": 0.0}
    try:
        xyz1, box1, vel1 = load_coords(p_rst, "equilibrated.rst7", natom)
    except CandidateError as e:
        ctx.info["equilibrated.error"] = str(e)
        xyz1 = box1 = vel1 = None
    if xyz1 is not None and box1 is not None and len(box1) >= 6 and np.all(np.isfinite(box1)) and \
            np.all(box1[:3] > 1.0) and np.all(np.abs(box1[3:6] - 90.0) < 0.01):
        eq["readable_with_box"] = 1.0
        dens = float(mass.sum() / (box1[0] * box1[1] * box1[2]) * 1.66054)
        ctx.info["equilibrated.density_g_cm3"] = round(dens, 4)
        eq["density"] = 1.0 if DENSITY_BAND[0] <= dens <= DENSITY_BAND[1] else 0.0
        if vel1 is not None and np.all(np.isfinite(vel1)) and system is not None:
            # ParmEd returns Amber-restart velocities in A/ps; kinetic energy in kJ/mol
            ke = 0.5 * float((mass[:, None] * (vel1 * 0.1) ** 2).sum())
            # Were bonds to hydrogen constrained in the run that wrote these velocities? Then the relative
            # velocity of the two atoms along the bond is ~0 (measured 0.8-1.2 A/ps) instead of thermal (~16 A/ps).
            xh = np.array([(b.atom1.idx, b.atom2.idx) for b in parm.bonds
                           if (Z[b.atom1.idx] == 1) != (Z[b.atom2.idx] == 1)], dtype=int).reshape(-1, 2)
            along = 0.0
            if len(xh):
                rvec = xyz1[xh[:, 0]] - xyz1[xh[:, 1]]
                dv = vel1[xh[:, 0]] - vel1[xh[:, 1]]
                along = float(np.sqrt(np.mean((np.einsum("ij,ij->i", dv, rvec) / np.linalg.norm(rvec, axis=1)) ** 2)))
            constrained = along < XH_ALONG_BOND_MAX
            dof = 3 * natom - 3 - (system.getNumConstraints() if constrained else 0)
            t_kin = 2 * ke / (dof * KB)
            ctx.info["equilibrated.kinetic_temperature_K"] = round(t_kin, 1)
            ctx.info["equilibrated.xh_velocity_along_bond_A_ps"] = round(along, 3)
            if TEMP_BAND[0] <= t_kin <= TEMP_BAND[1]:
                eq["velocities_thermal"] = 1.0
        _, _, rmsd1 = kabsch(xyz1[[i for _, i in ca_pairs]], t_ca)
        ctx.info["equilibrated.ca_rmsd"] = round(rmsd1, 3)
        eq["protein_intact"] = 1.0 if rmsd1 <= EQ_CA_RMSD else 0.0
        if sim_ok:
            try:
                ctxt, integ = context(xyz1, box1, 2)
                ctxt.setVelocitiesToTemperature(300 * u.kelvin, 3)
                e_start = pot(ctxt)
                integ.step(150)
                st = ctxt.getState(getEnergy=True, getPositions=True)
                e_end = st.getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)
                ke_end = st.getKineticEnergy().value_in_unit(u.kilojoule_per_mole)
                t_end = 2 * ke_end / ((3 * natom - system.getNumConstraints() - 3) * KB)
                moved = np.array(st.getPositions(asNumpy=True).value_in_unit(u.angstrom)) - xyz1
                max_move = float(np.abs(moved).max())
                ctx.info["equilibrated.md"] = {"e_start_kj": round(e_start, 0), "e_end_kj": round(e_end, 0),
                                               "t_end_K": round(t_end, 1), "max_move_A": round(max_move, 2)}
                if math.isfinite(e_start) and math.isfinite(e_end) and e_start < 0 and e_end < 0 and \
                        MD_TEMP_BAND[0] <= t_end <= MD_TEMP_BAND[1] and max_move < 8.0:
                    eq["md_continues_stably"] = 1.0
            except Exception as e:
                ctx.info["equilibrated.md_error"] = f"{type(e).__name__}: {str(e)[:200]}"
    for k, v in eq.items():
        put("equilibrated", k, v)

    # ------------------------------------------------------------------ report.json vs recomputed values
    rep = report if isinstance(report, dict) else {}
    rp = {"disulfides": 0.0, "protein_net_charge": 0.0, "n_na": 0.0, "n_cl": 0.0, "n_water": 0.0,
          "box_angstrom": 0.0}

    def norm_label(x):
        s = str(x).strip().upper()
        if ":" in s:
            s = s.split(":")[-1]
        return s
    ds = rep.get("disulfides")
    if isinstance(ds, list) and len(ds) <= 100:
        try:
            got = {frozenset(norm_label(x) for x in p) for p in ds if isinstance(p, (list, tuple)) and len(p) == 2}
            wantl = {frozenset(norm_label(x) for x in p) for p in truth["disulfide_labels"]}
            rp["disulfides"] = 1.0 if (got == wantl and len(got) == len(ds)) else 0.0
        except Exception:
            pass

    def num(x):
        if isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x):
            return None
        return float(x)
    v = num(rep.get("protein_net_charge"))
    rp["protein_net_charge"] = 1.0 if (v is not None and abs(v - q_solute) <= 0.05) else 0.0
    for key, val in (("n_na", n_na), ("n_cl", n_cl), ("n_water", n_wat)):
        v = num(rep.get(key))
        rp[key] = 1.0 if (v is not None and v == val) else 0.0
    bx = rep.get("box_angstrom")
    if isinstance(bx, list) and len(bx) == 3 and all(num(x) is not None for x in bx) and box1 is not None:
        rp["box_angstrom"] = 1.0 if all(abs(float(bx[k]) - box1[k]) <= BOX_REPORT_TOL for k in range(3)) else 0.0
    for k, v in rp.items():
        put("report", k, v)

    # ------------------------------------------------------------------ components
    ready = True
    for comp in [c for c in WEIGHTS if c != "production_ready"]:
        vals = checks[comp]
        value = sum(vals.values()) / len(vals)
        failed = [k for k, v in vals.items() if v < 1.0]
        ready = ready and not failed
        ctx.component(comp, value, WEIGHTS[comp], "failed: " + ", ".join(failed) if failed else "all checks passed")
    ctx.component("production_ready", 1.0 if ready else 0.0, WEIGHTS["production_ready"],
                  "1 only if every check of every other component passed")
    ctx.info["counts"] = {"n_na": n_na, "n_cl": n_cl, "n_water": n_wat, "natom": natom,
                          "protein_charge": round(q_solute, 3)}


def other_res_foreign(other_res, Z):
    """True if a solute residue is neither an amino acid nor a small blocking group (<= 3 heavy atoms of C/N/O)."""
    for r in other_res:
        zs = sorted(Z[a.idx] for a in r.atoms if Z[a.idx] > 1)
        if zs not in ([6, 7], [6, 6, 8]):
            return True
    return False


run(main, pass_threshold=0.9)
