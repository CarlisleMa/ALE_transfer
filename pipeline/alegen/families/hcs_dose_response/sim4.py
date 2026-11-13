"""Level-4 additions to the plate simulator (numpy only). Levels 1-3 never import this module.

What changes against sim.py (all of it planted before a pixel exists):
  * curve shapes   one count curve is a deep partial response (plateau 0.28-0.34, so the 0.5 crossing
                   lies far to the right of the inflection point); the marker per cell of the pure
                   cytotoxic compound rises in proportion to the fraction of cells lost
  * nuclear size   two compounds enlarge the nuclei (up to 1.3-1.5x linear) as their count effect grows
  * colonies       nuclei overlap (centre distance down to PARAMS["dmin"] px for nuclei of 12-15 x 10-12 px), the
                   intensity profile is flat-topped, overlapping nuclei are composited by maximum:
                   there is no intensity dip between them, only the outline shows that there are two
  * chromatin      stronger texture and up to two dark nucleoli per nucleus
  * fragments      small bright nuclear fragments in every well (a fraction of the seeded cells) and
                   many more where cells were lost; the SOP does not count them
  * field border   nuclei are cut by the border; the SOP's counting frame counts those on the top and
                   left edge and not those on the bottom and right edge
The planted cell number of a well is the number of rendered nuclei that the SOP counts.
"""
from __future__ import annotations

import numpy as np

from . import sim
from .sim import GAIN, M0, N0, NCOL, NUC_I, OFFSET, READ_SD, SIZE, TOP_DOSE, Compound, Well, active, gauss_blur, is_edge

PAD = 26                         # canvas margin so that nuclei can extend beyond the field
PX_UM = 1.3                      # micrometres per pixel (stated in the assay notes)
PARAMS = dict(
    dmin=9.5,                    # minimal centre distance of two nuclei (x size factor), px
    p_colony=0.90, p_close=0.55, close_span=2.5,
    dome=0.22, edge=0.10,        # flat-topped profile: (1 - dome*rho^4), linear roll-off over `edge`
    tex_fine=0.20, tex_coarse=0.18, n_nucleoli=(0, 2), nucl_depth=(0.3, 0.6), nucl_sigma=(0.8, 1.3),
    inten_sd=0.25, size_sd=0.05,
    deb_base=0.12, deb_lost=0.40,          # fragments per seeded cell / per lost cell
    deb_r=(1.3, 2.5), deb_i=(1.3, 2.8),    # fragment radius (px) and peak level (x NUC_I)
    deb_margin=7.0, deb_sep=8.0, deb_gap=3.5,   # px to the image edge / between fragments / to a nuclear outline
    gap_in=3.0, cut_min=2.5, vis_min=5.0,  # a nucleus stays >= gap_in px inside an edge or extends >= cut_min px
                                           # beyond it; a cut nucleus shows >= vis_min px of its extent
)


# ----------------------------------------------------------------------------------- plate design
def reshape_compounds(rng: np.random.Generator, cs: list[Compound]) -> dict:
    """Level-4 curve shapes and phenotypes on top of the level-3 compound set."""
    u = rng.uniform
    mp = next(c for c in cs if c.cls == "marker_potent")
    mp.count = active(u(1.2, 3.0), u(1.1, 1.4), u(0.28, 0.34))      # deep partial count response
    cy = next(c for c in cs if c.cls == "cytotoxic")
    pc = next(c for c in cs if c.cls == "partial_count")
    return {"marker_rise": {cy.name: float(u(0.5, 0.8))},
            "size_max": {mp.name: float(u(1.35, 1.5)), pc.name: float(u(1.3, 1.45))}}


def plant_wells4(rng: np.random.Generator, wells: list[Well], cs: list[Compound], shapes: dict,
                 q: dict = PARAMS) -> dict:
    """Planted per-well quantities: sim.plant_wells at level 3 plus marker rise, nuclear size factor
    and number of fragments."""
    by_name = {c.name: c for c in cs}
    plate = {"stain_end": float(rng.uniform(0.36, 0.44)), "stain_dir": int(rng.choice([-1, 1])),
             "row_tilt": float(rng.uniform(-0.10, 0.10)), "edge_count": float(rng.uniform(0.50, 0.58)),
             "beta": 0.35}
    g = np.log(plate["stain_end"])
    for w in wells:
        x = w.col / (NCOL - 1.0)
        if plate["stain_dir"] < 0:
            x = 1.0 - x
        w.stain_factor = float(np.exp(g * x + plate["row_tilt"] * (w.row - 3.5) / 7.0))
        w.count_factor = plate["edge_count"] if is_edge(w.row, w.col) else 1.0
        size = 1.0
        if w.role == "vehicle":
            fm, fc, af = 1.0, 1.0, 0.0
        elif w.role == "reference":
            fm, fc, af = 0.05, 1.0, 0.0
        else:
            cp = by_name[w.compound]
            fm, fc = float(cp.marker(w.dose)), float(cp.count(w.dose))
            af = cp.autofluor * w.dose / TOP_DOSE
            if cp.name in shapes["marker_rise"]:          # survivors carry more marker
                fm = fm * (1.0 + shapes["marker_rise"][cp.name] * (1.0 - fc))
            if cp.name in shapes["size_max"]:             # nuclei grow with the count effect
                eff = (1.0 - fc) / (1.0 - cp.count.bottom)
                size = 1.0 + (shapes["size_max"][cp.name] - 1.0) * min(1.0, max(0.0, eff))
        seeded = N0 * w.count_factor * rng.lognormal(0.0, 0.03)
        w.n_cells = int(round(seeded * fc))
        w.extra["n_seeded"] = int(round(seeded))
        w.extra["marker_level"] = float(M0 * w.stain_factor * fm * rng.lognormal(0.0, 0.03))
        w.extra["size"] = float(size)
        lost = max(0, w.extra["n_seeded"] - w.n_cells)
        w.extra["n_debris"] = int(rng.poisson(q["deb_base"] * w.extra["n_seeded"] + q["deb_lost"] * lost))
        bg = plate["beta"] * M0 * w.stain_factor * rng.lognormal(0.0, 0.15)
        w.background = float(bg + af * M0)
    return plate


def choose_annotated(rng: np.random.Generator, wells: list[Well], cs: list[Compound]) -> list[str]:
    """Fields the lab annotated by hand: untreated dense and sparse fields and the treated phenotypes
    (fragment-rich, enlarged nuclei). All in focus."""
    ok = [w for w in wells if w.defocus == 0]
    cls = {c.cls: c for c in cs}

    def series(c):
        return sorted((w for w in ok if w.compound == c.name), key=lambda w: w.dose)

    def pick(pool, k):
        idx = rng.permutation(len(pool))[:k]
        return [pool[int(i)].well for i in idx]

    out = []
    out += pick([w for w in ok if w.role == "vehicle" and not is_edge(w.row, w.col)], 2)
    out += pick([w for w in ok if w.role == "vehicle" and is_edge(w.row, w.col)], 1)
    cy = series(cls["cytotoxic"])
    out += [cy[-1].well, min(cy, key=lambda w: abs(np.log(w.dose / cls["cytotoxic"].count.ic50))).well]
    out += [series(cls["marker_potent"])[-1].well, series(cls["partial_count"])[-1].well]
    out += [series(cls["dual"])[-2].well]
    assert len(set(out)) == 8, out
    return sorted(out)


# ----------------------------------------------------------------------------------- rendering
def _extents(a, b, th):
    c, s = np.cos(th), np.sin(th)
    return np.sqrt(a * a * c * c + b * b * s * s), np.sqrt(a * a * s * s + b * b * c * c)


def _edge_state(c, e, q):
    """(ok, cut_low, cut_high) of one axis. A nucleus either lies clearly inside an edge or clearly
    extends beyond it, and a cut nucleus shows at least `vis_min` px of its extent."""
    lo = c - e + 0.5                       # > 0: inside the low edge (top / left)
    hi = (SIZE - 0.5) - (c + e)            # > 0: inside the high edge (bottom / right)
    ok = (lo >= q["gap_in"] or (lo <= -q["cut_min"] and c + e + 0.5 >= q["vis_min"])) and \
         (hi >= q["gap_in"] or (hi <= -q["cut_min"] and (SIZE - 0.5) - (c - e) >= q["vis_min"]))
    return ok, lo < 0, hi < 0


def place_cells4(rng: np.random.Generator, n_counted: int, size: float, q: dict) -> np.ndarray:
    """Rows: cx, cy, a, b, theta, counted. Nuclei cut by the bottom or right edge are placed but not
    counted."""
    dmin = q["dmin"] * size
    cells = np.zeros((4 * max(n_counted, 1) + 50, 6))
    k = counted = tries = 0
    while counted < n_counted and tries < 600000 and k < len(cells):
        tries += 1
        if k and rng.random() < q["p_colony"]:
            base = cells[int(rng.integers(0, k)), :2]
            d = (rng.uniform(dmin, dmin + q["close_span"] * size) if rng.random() < q["p_close"]
                 else rng.uniform(14.0 * size, 20.0 * size))
            ang = rng.uniform(0, 2 * np.pi)
            cand = base + d * np.array([np.cos(ang), np.sin(ang)])
        else:
            cand = rng.uniform(-4.0, SIZE + 3.0, size=2)
        if not (-4.0 <= cand[0] <= SIZE + 3.0 and -4.0 <= cand[1] <= SIZE + 3.0):
            continue
        if k and np.min(np.hypot(cells[:k, 0] - cand[0], cells[:k, 1] - cand[1])) < dmin:
            continue
        s = size * float(np.clip(rng.lognormal(0.0, q["size_sd"]), 0.9, 1.12))   # shortest axis stays >= 8.8 px (11.5 um)
        a, b, th = rng.uniform(6.0, 7.6) * s, rng.uniform(4.9, 6.1) * s, rng.uniform(0, np.pi)
        ex, ey = _extents(a, b, th)
        okx, _, cut_right = _edge_state(cand[0], ex, q)
        oky, _, cut_bottom = _edge_state(cand[1], ey, q)
        if not (okx and oky):
            continue
        is_counted = not (cut_right or cut_bottom)
        cells[k] = (cand[0], cand[1], a, b, th, float(is_counted))
        k += 1
        counted += is_counted
    return cells[:k]


def render_well4(rng: np.random.Generator, w: Well, illum1: np.ndarray, illum2: np.ndarray, q: dict = PARAMS):
    """Returns (ch1, ch2, points). points = centres of the nuclei the SOP counts (the planted truth
    behind the manual annotation). Sets w.n_cells, w.marker_mean and w.extra['n_debris'] to what was
    rendered."""
    size = w.extra["size"]
    cells = place_cells4(rng, max(w.extra["n_seeded"], w.n_cells), size, q)
    cidx = np.flatnonzero(cells[:, 5] > 0)
    uidx = np.flatnonzero(cells[:, 5] == 0)
    keep_c = rng.permutation(cidx)[: w.n_cells]
    frac = w.n_cells / max(len(cidx), 1)
    keep_u = uidx[rng.random(len(uidx)) < frac]
    cells = cells[np.sort(np.concatenate([keep_c, keep_u]))]
    counted = cells[:, 5] > 0
    w.n_cells = int(counted.sum())                    # planted count = rendered nuclei the SOP counts
    w.extra["n_cut_uncounted"] = int((~counted).sum())

    C = SIZE + 2 * PAD
    yy, xx = np.mgrid[0:C, 0:C]
    nuc = np.zeros((C, C))
    cyt = np.zeros((C, C))
    zbuf = np.full((C, C), np.inf)
    texf = gauss_blur(rng.normal(0, 1, (C, C)), 1.0)
    texf /= texf.std()
    texc = gauss_blur(rng.normal(0, 1, (C, C)), 2.4)
    texc /= texc.std()
    tex2 = gauss_blur(rng.normal(0, 1, (C, C)), 2.0)
    tex2 /= tex2.std()
    levels = w.extra["marker_level"] * rng.lognormal(-0.5 * 0.25 ** 2, 0.25, size=len(cells))
    w.marker_mean = float(levels[counted].mean()) if counted.any() else 0.0
    for i, (cx, cy, a, b, th, _) in enumerate(cells):
        cx, cy = cx + PAD, cy + PAD
        R = int(np.ceil(a + 11.0))
        x0, x1 = int(max(0, cx - R)), int(min(C, cx + R + 1))
        y0, y1 = int(max(0, cy - R)), int(min(C, cy + R + 1))
        dx, dy = xx[y0:y1, x0:x1] - cx, yy[y0:y1, x0:x1] - cy
        pu, pv = dx * np.cos(th) + dy * np.sin(th), -dx * np.sin(th) + dy * np.cos(th)   # px along the axes
        rho = np.sqrt((pu / a) ** 2 + (pv / b) ** 2)
        inten = NUC_I * rng.lognormal(0, q["inten_sd"])
        prof = (1.0 - q["dome"] * np.clip(rho, 0, 1) ** 4) * np.clip((1.03 - rho) / q["edge"], 0, 1)
        tex = 1.0 + q["tex_fine"] * texf[y0:y1, x0:x1] + q["tex_coarse"] * texc[y0:y1, x0:x1]
        for _k in range(int(rng.integers(q["n_nucleoli"][0], q["n_nucleoli"][1] + 1))):
            rr, ang = 0.55 * np.sqrt(rng.random()), rng.uniform(0, 2 * np.pi)
            nu, nv = rr * a * np.cos(ang), rr * b * np.sin(ang)
            sg = rng.uniform(*q["nucl_sigma"]) * size
            tex = tex * (1.0 - rng.uniform(*q["nucl_depth"]) * np.exp(-((pu - nu) ** 2 + (pv - nv) ** 2) / (2 * sg * sg)))
        stamp = inten * prof * np.clip(tex, 0.05, None)
        nuc[y0:y1, x0:x1] = np.maximum(nuc[y0:y1, x0:x1], stamp)
        # cytoplasm: as in sim.render_well, but wider and flatter, so that the zone next to the nucleus has
        # the same level in a packed colony (where neighbours fill it) and around a lone cell
        ca, cb = a + rng.uniform(5.5, 7.5), b + rng.uniform(5.5, 7.5)
        ox, oy = rng.normal(0, 0.8, size=2)
        u2 = ((dx - ox) * np.cos(th) + (dy - oy) * np.sin(th)) / ca
        v2 = (-(dx - ox) * np.sin(th) + (dy - oy) * np.cos(th)) / cb
        rho2 = np.sqrt(u2 * u2 + v2 * v2)
        body = np.clip((1.0 - rho2) / 0.12, 0, 1) * (0.92 + 0.08 * (1 - np.clip(rho2, 0, 1)))
        body = body * np.where(rho < 0.92, 0.35, 1.0)
        zb = zbuf[y0:y1, x0:x1]
        win = (rho2 < zb) & (rho2 < 1.0)
        sub = cyt[y0:y1, x0:x1]
        sub[win] = (levels[i] * body * (1.0 + 0.10 * tex2[y0:y1, x0:x1]))[win]
        zb[win] = rho2[win]

    # nuclear fragments: small, bright, outside the nuclei, away from the field border and from each other
    deb = np.zeros((0, 2))
    tries = 0
    lo, hi = q["deb_margin"], SIZE - 1 - q["deb_margin"]
    while len(deb) < w.extra["n_debris"] and tries < 20000:
        tries += 1
        p = rng.uniform(lo, hi, size=2)
        if len(deb) and np.min(np.hypot(deb[:, 0] - p[0], deb[:, 1] - p[1])) < q["deb_sep"]:
            continue
        if len(cells):
            ddx, ddy = p[0] - cells[:, 0], p[1] - cells[:, 1]
            cu = (ddx * np.cos(cells[:, 4]) + ddy * np.sin(cells[:, 4])) / (cells[:, 2] + q["deb_gap"])
            cv = (-ddx * np.sin(cells[:, 4]) + ddy * np.cos(cells[:, 4])) / (cells[:, 3] + q["deb_gap"])
            if np.min(cu * cu + cv * cv) < 1.0:
                continue
        deb = np.vstack([deb, p])
        r = rng.uniform(*q["deb_r"])
        ecc, th = rng.uniform(0.8, 1.0), rng.uniform(0, np.pi)
        px, py = p[0] + PAD, p[1] + PAD
        x0, x1, y0, y1 = int(px - 6), int(px + 7), int(py - 6), int(py + 7)
        dx, dy = xx[y0:y1, x0:x1] - px, yy[y0:y1, x0:x1] - py
        rho = np.sqrt(((dx * np.cos(th) + dy * np.sin(th)) / r) ** 2
                      + ((-dx * np.sin(th) + dy * np.cos(th)) / (r * ecc)) ** 2)
        stamp = NUC_I * rng.uniform(*q["deb_i"]) * np.clip((1.05 - rho) / 0.25, 0, 1)
        nuc[y0:y1, x0:x1] = np.maximum(nuc[y0:y1, x0:x1], stamp)
    w.extra["n_debris"] = int(len(deb))

    sig1 = gauss_blur(nuc, 1.0) + 25.0
    sig2 = gauss_blur(cyt, 1.15) + w.background
    if w.defocus > 0:
        sig1 = gauss_blur(sig1, w.defocus)
        sig2 = gauss_blur(sig2, w.defocus)
    out = []
    for sig, il in ((sig1, illum1), (sig2, illum2)):
        e = np.clip(sig[PAD:PAD + SIZE, PAD:PAD + SIZE] * il, 0, None) * GAIN
        img = rng.poisson(e) / GAIN + OFFSET + rng.normal(0, READ_SD, e.shape)
        out.append(np.clip(np.round(img), 0, 65535).astype(np.uint16))
    return out[0], out[1], cells[counted, :2].copy()


def annotation_tsv(rng: np.random.Generator, points: np.ndarray) -> str:
    """Manual point annotation of one field: one click per counted nucleus, within about a pixel of
    its centre and always inside the image."""
    p = np.clip(points + rng.normal(0, 0.8, size=points.shape), 1.0, SIZE - 2.0)
    order = np.lexsort((p[:, 0], np.round(p[:, 1] / 16.0)))      # roughly row by row, as a person clicks
    return "x_px\ty_px\n" + "".join(f"{p[i, 0]:.1f}\t{p[i, 1]:.1f}\n" for i in order)


# ----------------------------------------------------------------------------------- whole plate
def build_plate(seed: int, workdir, q: dict = PARAMS) -> dict:
    """Plant and render a level-4 plate into `workdir` (images/, flatfield/, annotations/). Returns
    the planted objects and the map of input files."""
    from pathlib import Path

    import tifffile

    workdir = Path(workdir)
    rng = np.random.default_rng(884000 + seed)          # level 4 does not replay the level-3 plate of the seed
    cs = sim.draw_compounds(rng, 3)
    sim.assign_rows(rng, cs, 3)
    shapes = reshape_compounds(rng, cs)
    wells, ascending = sim.layout(rng, cs)
    plate = plant_wells4(rng, wells, cs, shapes, q)
    blurred = sim.plant_defocus(rng, wells, cs, 3)
    annotated = choose_annotated(rng, wells, cs)

    dirs = {k: workdir / k for k in ("images", "flatfield", "annotations")}
    for d in dirs.values():
        if d.exists():
            for f in list(d.glob("*.tif")) + list(d.glob("*.tsv")):
                f.unlink()
        d.mkdir(parents=True, exist_ok=True)
    il1, il2 = sim.illumination(0.30), sim.illumination(0.26)
    arng = np.random.default_rng(884500 + seed)
    files: dict = {}
    points = {}
    for w in wells:
        c1, c2, pts = render_well4(rng, w, il1, il2, q)
        points[w.well] = pts
        for ch, im in ((1, c1), (2, c2)):
            p = dirs["images"] / f"{w.well}_ch{ch}.tif"
            tifffile.imwrite(p, im, compression="zlib", compressionargs={"level": 9}, predictor=True)
            files[f"images/{p.name}"] = p
        if w.well in annotated:
            p = dirs["annotations"] / f"{w.well}.tsv"
            p.write_text(annotation_tsv(arng, pts))
            files[f"annotations/{p.name}"] = p
    f1, f2, dk = sim.reference_images(rng, il1, il2)
    for name, im in (("flat_ch1", f1), ("flat_ch2", f2), ("dark", dk)):
        p = dirs["flatfield"] / f"{name}.tif"
        tifffile.imwrite(p, im, compression="zlib", predictor=True)
        files[f"flatfield/{p.name}"] = p
    return {"cs": cs, "wells": wells, "ascending": ascending, "plate": plate, "shapes": shapes,
            "blurred": blurred, "annotated": annotated, "files": files, "points": points}
