"""Forward simulator of a 96-well, 2-channel high-content screening plate (numpy only).

Everything is planted: the potency of every compound on both readouts, the number of cells of
every well, the marker level of every cell, the plate-position factors, the background of every
well and which wells were acquired out of focus. Images are rendered from those values.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

ROWS = "ABCDEFGH"
NCOL = 12
SIZE = 320                       # pixels per side
DOSES = [30.0 / 3 ** i for i in range(7, -1, -1)]     # 0.0137 ... 30 uM, ascending
TOP_DOSE = DOSES[-1]
REF_DOSE = 10.0
OFFSET = 100.0                   # camera offset (grey levels)
GAIN = 2.0                       # electrons per grey level
READ_SD = 2.5
NUC_I = 1400.0                   # nuclear stain level
M0 = 600.0                       # cytoplasmic marker level of a vehicle cell at plate factor 1
N0 = 118                         # cells per field in an unaffected vehicle well


def hill(d, ec50, h, bottom):
    return bottom + (1.0 - bottom) / (1.0 + (np.asarray(d, dtype=float) / ec50) ** h)


def ec50_from_abs(ic50, h, bottom):
    """EC50 (inflection) of a top=1 logistic whose 50%-of-vehicle crossing is `ic50`."""
    return ic50 / ((1.0 - bottom) / (0.5 - bottom) - 1.0) ** (1.0 / h)


@dataclass
class Curve:
    """Planted response of one readout: fraction of the vehicle level as a function of dose."""
    kind: str                    # "active" | "none" | "partial" | "beyond_range"
    ic50: float | None = None    # concentration with response 0.5 (None unless active)
    ec50: float = 1.0
    h: float = 1.0
    bottom: float = 1.0

    def __call__(self, d):
        if self.kind == "none":
            return np.ones_like(np.asarray(d, dtype=float))
        return hill(d, self.ec50, self.h, self.bottom)


def active(ic50, h, bottom):
    return Curve("active", float(ic50), float(ec50_from_abs(ic50, h, bottom)), float(h), float(bottom))


@dataclass
class Compound:
    name: str
    cls: str
    marker: Curve
    count: Curve
    autofluor: float = 0.0       # background added at the top dose, in units of M0
    row: int = -1


@dataclass
class Well:
    well: str
    row: int
    col: int                     # 0-based
    role: str                    # vehicle | reference | compound
    compound: str
    dose: float
    n_cells: int = 0
    marker_mean: float = 0.0     # mean planted cytoplasmic marker level of the cells (grey levels)
    stain_factor: float = 1.0
    count_factor: float = 1.0
    background: float = 0.0      # diffuse marker-channel background incl. compound fluorescence
    defocus: float = 0.0         # extra blur sigma (px); 0 = in focus
    extra: dict = field(default_factory=dict)


# ----------------------------------------------------------------------------------- plate design
def draw_compounds(rng: np.random.Generator, level: int) -> list[Compound]:
    u = rng.uniform
    none = Curve("none")
    ids = sorted(int(x) for x in rng.permutation(np.arange(1000, 9999))[:8])
    names = [f"CPD-{i}" for i in ids]
    rng.shuffle(names)
    cs: list[Compound] = []
    # 1 potent marker inhibitor, cytotoxic only far above its marker potency
    cs.append(Compound(names[0], "marker_potent", active(u(0.045, 0.11), u(0.9, 1.2), u(0.0, 0.05)),
                       active(u(2.5, 5.5), u(1.0, 1.3), u(0.12, 0.2))))
    # 2 weak marker inhibitor
    cs.append(Compound(names[1], "marker_weak", active(u(2.5, 5.5), u(0.9, 1.2), u(0.0, 0.05)), none))
    if level >= 2:
        # 3 cytotoxic, marker per cell unchanged
        cs.append(Compound(names[2], "cytotoxic", none, active(u(0.25, 1.2), u(0.9, 1.2), u(0.1, 0.2))))
        # 4 cytotoxic at low dose, marker inhibited only at >= 12x higher dose
        gi = u(0.08, 0.3)
        cs.append(Compound(names[3], "dual", active(min(gi * u(12, 25), 5.5), u(1.0, 1.3), u(0.0, 0.05)),
                           active(gi, u(1.0, 1.3), u(0.20, 0.30))))
    else:
        gi = u(0.6, 1.5)
        cs.append(Compound(names[2], "dual_marker_first", active(gi / u(3, 6), u(0.9, 1.2), u(0.0, 0.05)),
                           active(gi, u(0.9, 1.2), u(0.1, 0.2))))
        gi = u(2.5, 5.5)
        cs.append(Compound(names[3], "dual_marker_first", active(gi / u(10, 25), u(1.0, 1.3), u(0.0, 0.05)),
                           active(gi, u(1.0, 1.3), u(0.20, 0.30))))
    # 5 marker inhibitor that fluoresces in the marker channel (level 3)
    cs.append(Compound(names[4], "marker_fluorescent" if level >= 3 else "marker_mid",
                       active(u(2.5, 5.5), u(0.9, 1.2), u(0.0, 0.05)), none,
                       autofluor=float(u(1.3, 2.0)) if level >= 3 else 0.0))
    # 6 no effect
    cs.append(Compound(names[5], "inactive", none, none))
    # 7 partial marker inhibitor: plateau above 50%
    hb, bb = u(1.0, 1.5), u(0.62, 0.72)
    cs.append(Compound(names[6], "partial_marker", Curve("partial", None, float(u(0.1, 1.0)), float(hb), float(bb)),
                       none))
    # 8 cytostatic plateau above 50%; marker potency beyond the tested range
    hb, bb = u(1.0, 1.5), u(0.60, 0.65)
    cs.append(Compound(names[7], "partial_count",
                       Curve("beyond_range", None, float(u(150, 400)), 1.0, 0.0),
                       Curve("partial", None, float(u(0.1, 1.0)), float(hb), float(bb))))
    return cs


def assign_rows(rng: np.random.Generator, cs: list[Compound], level: int) -> None:
    """One compound per row. From level 2 on, the pure cytotoxic compound and the cytostatic
    compound occupy the two outer rows (seed decides which)."""
    idx = list(range(8))
    if level >= 2:
        cyto = next(i for i, c in enumerate(cs) if c.cls == "cytotoxic")
        part = next(i for i, c in enumerate(cs) if c.cls == "partial_count")
        rest = [i for i in idx if i not in (cyto, part)]
        rest = [rest[i] for i in rng.permutation(len(rest))]
        outer = [cyto, part] if rng.random() < 0.5 else [part, cyto]
        order = [outer[0]] + rest + [outer[1]]
    else:
        order = [idx[i] for i in rng.permutation(8)]
    for r, i in enumerate(order):
        cs[i].row = r


def layout(rng: np.random.Generator, cs: list[Compound]) -> tuple[list[Well], bool]:
    """Rows = compounds. Each row holds 3 vehicle wells (every column ends up with exactly 2),
    one reference-inhibitor well and the 8-point series in column order."""
    offs = [0, 1, 2, 3, 0, 1, 2, 3]
    offs = [offs[i] for i in rng.permutation(8)]
    ascending = bool(rng.random() < 0.5)
    by_row = {c.row: c for c in cs}
    wells: list[Well] = []
    for r in range(8):
        veh = {(offs[r] + 4 * j) % NCOL for j in range(3)}
        ref = (offs[r] + 2 + 4 * int(rng.integers(0, 3))) % NCOL
        free = [c for c in range(NCOL) if c not in veh and c != ref]
        doses = DOSES if ascending else DOSES[::-1]
        dose_of = dict(zip(free, doses))
        for c in range(NCOL):
            name = f"{ROWS[r]}{c + 1:02d}"
            if c in veh:
                wells.append(Well(name, r, c, "vehicle", "DMSO", 0.0))
            elif c == ref:
                wells.append(Well(name, r, c, "reference", "REF-INH", REF_DOSE))
            else:
                wells.append(Well(name, r, c, "compound", by_row[r].name, float(dose_of[c])))
    return wells, ascending


def is_edge(r: int, c: int) -> bool:
    return r in (0, 7) or c in (0, NCOL - 1)


def plant_wells(rng: np.random.Generator, wells: list[Well], cs: list[Compound], level: int) -> dict:
    """Fill in the planted per-well quantities. Returns the plate-level parameters."""
    by_name = {c.name: c for c in cs}
    plate = {"stain_end": 1.0, "stain_dir": 0, "row_tilt": 0.0, "edge_count": 1.0, "beta": 0.35}
    if level >= 3:
        plate["stain_end"] = float(rng.uniform(0.36, 0.44))        # factor at the far end of the gradient
        plate["stain_dir"] = int(rng.choice([-1, 1]))              # +1: decreasing with column number
        plate["row_tilt"] = float(rng.uniform(-0.10, 0.10))        # log factor across the 8 rows
        plate["edge_count"] = float(rng.uniform(0.50, 0.58))       # cell-number factor of the outer wells
    g = np.log(plate["stain_end"])
    for w in wells:
        x = w.col / (NCOL - 1.0)
        if plate["stain_dir"] < 0:
            x = 1.0 - x
        w.stain_factor = float(np.exp(g * x + plate["row_tilt"] * (w.row - 3.5) / 7.0)) if level >= 3 else 1.0
        w.count_factor = plate["edge_count"] if is_edge(w.row, w.col) else 1.0
        if w.role == "vehicle":
            fm, fc, af = 1.0, 1.0, 0.0
        elif w.role == "reference":
            fm, fc, af = 0.05, 1.0, 0.0
        else:
            cp = by_name[w.compound]
            fm, fc = float(cp.marker(w.dose)), float(cp.count(w.dose))
            af = cp.autofluor * w.dose / TOP_DOSE
        seeded = N0 * w.count_factor * rng.lognormal(0.0, 0.03)
        w.n_cells = int(round(seeded * fc))
        w.extra["n_seeded"] = int(round(seeded))
        w.extra["marker_level"] = float(M0 * w.stain_factor * fm * rng.lognormal(0.0, 0.03))
        bg = plate["beta"] * M0 * w.stain_factor * (rng.lognormal(0.0, 0.15) if level >= 3 else 1.0)
        w.background = float(bg + af * M0)
    return plate


def plant_defocus(rng: np.random.Generator, wells: list[Well], cs: list[Compound], level: int) -> list[str]:
    """Out-of-focus acquisitions (level >= 2): two runs of two neighbouring series wells (the two
    doses below the marker potency of the weak marker inhibitor / below the count potency of the
    potent marker inhibitor) and one vehicle well."""
    if level < 2:
        return []
    out = []
    for cls, readout in (("marker_weak", "marker"), ("marker_potent", "count")):
        cp = next(c for c in cs if c.cls == cls)
        curve = cp.marker if readout == "marker" else cp.count
        series = sorted((w for w in wells if w.compound == cp.name), key=lambda w: w.dose)
        below = [w for w in series if w.dose < curve.ic50 and float(curve(w.dose)) >= 0.72]
        for w in below[-2:]:
            w.defocus = float(rng.uniform(4.5, 6.0))
            out.append(w.well)
    used_rows = {w.row for w in wells if w.defocus > 0}
    veh = [w for w in wells if w.role == "vehicle" and w.row not in used_rows and not is_edge(w.row, w.col)]
    w = veh[int(rng.integers(0, len(veh)))]
    w.defocus = float(rng.uniform(4.5, 6.0))
    out.append(w.well)
    return sorted(out)


# ----------------------------------------------------------------------------------- rendering
def gauss_blur(img: np.ndarray, sigma: float) -> np.ndarray:
    if sigma <= 0:
        return img
    r = int(np.ceil(4 * sigma))
    k = np.exp(-0.5 * (np.arange(-r, r + 1) / sigma) ** 2)
    k /= k.sum()
    p = np.pad(img, r, mode="reflect")
    p = np.apply_along_axis(lambda v: np.convolve(v, k, mode="valid"), 0, p)
    p = np.apply_along_axis(lambda v: np.convolve(v, k, mode="valid"), 1, p)
    return p


def illumination(strength: float) -> np.ndarray:
    y, x = np.mgrid[0:SIZE, 0:SIZE]
    c = (SIZE - 1) / 2.0
    r2 = ((x - c) ** 2 + (y - c) ** 2) / (2 * c * c)
    return 1.0 - strength * r2


def place_cells(rng: np.random.Generator, n: int, touching: bool) -> np.ndarray:
    """Cell centres. `touching`: cells grow as colonies in which part of the nuclei abut their
    neighbours; otherwise every nucleus is separated from its neighbours by a clear gap."""
    margin = 10.0
    lo, hi = margin, SIZE - margin
    dmin = 10.5 if touching else 19.0
    pts = np.zeros((max(n, 1), 2))
    k = 0
    tries = 0
    while k < n and tries < 400000:
        tries += 1
        if touching and k and rng.random() < 0.88:
            base = pts[int(rng.integers(0, k))]
            d = rng.uniform(dmin, dmin + 2.0) if rng.random() < 0.45 else rng.uniform(14.0, 20.0)
            ang = rng.uniform(0, 2 * np.pi)
            cand = base + d * np.array([np.cos(ang), np.sin(ang)])
        else:
            cand = rng.uniform(lo, hi, size=2)
        if not (lo <= cand[0] <= hi and lo <= cand[1] <= hi):
            continue
        if k and np.min(np.hypot(pts[:k, 0] - cand[0], pts[:k, 1] - cand[1])) < dmin:
            continue
        pts[k] = cand
        k += 1
    return pts[:k]


def render_well(rng: np.random.Generator, w: Well, touching: bool, illum1: np.ndarray, illum2: np.ndarray):
    """Returns (ch1, ch2) uint16 and fills w.marker_mean with the mean planted marker of the cells."""
    n_seed = max(w.extra["n_seeded"], w.n_cells)
    pts = place_cells(rng, n_seed, touching)
    keep = rng.permutation(len(pts))[: w.n_cells]
    pts = pts[np.sort(keep)]
    w.n_cells = int(len(pts))                 # the planted count is the number of rendered cells
    yy, xx = np.mgrid[0:SIZE, 0:SIZE]
    nuc = np.zeros((SIZE, SIZE))
    cyt = np.zeros((SIZE, SIZE))
    zbuf = np.full((SIZE, SIZE), np.inf)
    tex1 = gauss_blur(rng.normal(0, 1, (SIZE, SIZE)), 1.3)
    tex1 /= tex1.std()
    tex2 = gauss_blur(rng.normal(0, 1, (SIZE, SIZE)), 2.0)
    tex2 /= tex2.std()
    levels = w.extra["marker_level"] * rng.lognormal(-0.5 * 0.25 ** 2, 0.25, size=len(pts))
    w.marker_mean = float(levels.mean()) if len(pts) else 0.0
    for i, (cx, cy) in enumerate(pts):
        a, b = rng.uniform(6.0, 7.6), rng.uniform(4.9, 6.1)
        th = rng.uniform(0, np.pi)
        x0, x1 = int(max(0, cx - 16)), int(min(SIZE, cx + 17))
        y0, y1 = int(max(0, cy - 16)), int(min(SIZE, cy + 17))
        dx, dy = xx[y0:y1, x0:x1] - cx, yy[y0:y1, x0:x1] - cy
        u = (dx * np.cos(th) + dy * np.sin(th)) / a
        v = (-dx * np.sin(th) + dy * np.cos(th)) / b
        rho = np.sqrt(u * u + v * v)
        inten = NUC_I * rng.lognormal(0, 0.18)
        prof = np.clip(1.0 - rho ** 4, 0, None) ** 0.5 * np.clip((1.02 - rho) / 0.12, 0, 1)
        stamp = inten * prof * (1.0 + 0.16 * tex1[y0:y1, x0:x1])
        nuc[y0:y1, x0:x1] = np.maximum(nuc[y0:y1, x0:x1], stamp)
        # cytoplasm: larger ellipse around the nucleus, nearest cell wins where territories meet
        ca, cb = a + rng.uniform(4.0, 6.0), b + rng.uniform(4.0, 6.0)
        ox, oy = rng.normal(0, 0.8, size=2)
        u2 = ((dx - ox) * np.cos(th) + (dy - oy) * np.sin(th)) / ca
        v2 = (-(dx - ox) * np.sin(th) + (dy - oy) * np.cos(th)) / cb
        rho2 = np.sqrt(u2 * u2 + v2 * v2)
        body = np.clip((1.0 - rho2) / 0.18, 0, 1) * (0.75 + 0.25 * (1 - np.clip(rho2, 0, 1)))
        body = body * np.where(rho < 0.92, 0.35, 1.0)
        zb = zbuf[y0:y1, x0:x1]
        win = (rho2 < zb) & (rho2 < 1.0)
        sub = cyt[y0:y1, x0:x1]
        sub[win] = (levels[i] * body * (1.0 + 0.10 * tex2[y0:y1, x0:x1]))[win]
        zb[win] = rho2[win]
    psf1, psf2 = 1.0, 1.15
    sig1 = gauss_blur(nuc, psf1) + 25.0
    sig2 = gauss_blur(cyt, psf2) + w.background
    if w.defocus > 0:
        sig1 = gauss_blur(sig1, w.defocus)
        sig2 = gauss_blur(sig2, w.defocus)
    out = []
    for sig, il in ((sig1, illum1), (sig2, illum2)):
        e = np.clip(sig * il, 0, None) * GAIN
        img = rng.poisson(e) / GAIN + OFFSET + rng.normal(0, READ_SD, e.shape)
        out.append(np.clip(np.round(img), 0, 65535).astype(np.uint16))
    return out[0], out[1]


def reference_images(rng: np.random.Generator, illum1: np.ndarray, illum2: np.ndarray):
    """Flat-field frames (uniform dye slide, average of 50 exposures) and the dark frame."""
    def avg(sig):
        e = sig * GAIN
        acc = np.zeros(sig.shape)
        for _ in range(8):
            acc += rng.poisson(e) / GAIN + rng.normal(0, READ_SD, sig.shape)
        return acc / 8 / np.sqrt(50 / 8) * 0 + sig + (acc / 8 - sig) / np.sqrt(50 / 8)
    f1 = np.round(avg(2400.0 * illum1) + OFFSET)
    f2 = np.round(avg(2100.0 * illum2) + OFFSET)
    dark = np.round(OFFSET + rng.normal(0, READ_SD / np.sqrt(50), illum1.shape))
    return f1.astype(np.uint16), f2.astype(np.uint16), dark.astype(np.uint16)
