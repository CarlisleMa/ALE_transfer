#!/usr/bin/env python3
"""Independent second solution: classical image analysis only (no Cellpose) and no code shared with
hcs_pipeline.py. Every step uses a different method than the oracle:

  nuclei      Laplacian-of-Gaussian seeds + intensity watershed   (oracle: Cellpose nuclei model)
  background  mode of the field histogram                         (oracle: median far from nuclei)
  marker      mean over a ring 1-4 px outside each nucleus        (oracle: 0-4 px ring)
  focus       band ratio of the power spectrum, noise-corrected   (oracle: edge-gradient percentile)
  plate model edge ratio from medians + Theil-Sen slopes          (oracle: least-squares regression)
  curve fit   exhaustive grid search                              (oracle: scipy least_squares)
"""
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi
from skimage.feature import peak_local_max
from skimage.segmentation import expand_labels, watershed

inp, out = Path("/task/input"), Path("/task/output")
out.mkdir(parents=True, exist_ok=True)
lay = pd.read_csv(inp / "plate_layout.tsv", sep="\t")
lay["r"] = lay["well"].str[0].map(lambda s: ord(s) - 65)
lay["c"] = lay["well"].str[1:].astype(int) - 1
lay["edge"] = lay["r"].isin([0, lay["r"].max()]) | lay["c"].isin([0, lay["c"].max()])

dark = tifffile.imread(inp / "flatfield" / "dark.tif").astype(np.float64)
gainmap = {}
for ch in (1, 2):
    f = tifffile.imread(inp / "flatfield" / f"flat_ch{ch}.tif").astype(np.float64) - dark
    f = ndi.uniform_filter(f, 9)
    gainmap[ch] = f / np.median(f)


def read(well, ch):
    return (tifffile.imread(inp / "images" / f"{well}_ch{ch}.tif").astype(np.float64) - dark) / gainmap[ch]


def mad_sd(a):
    return 1.4826 * np.median(np.abs(a - np.median(a)))


def nuclei(img):
    sm = ndi.gaussian_filter(img, 1.2)
    base, sd = np.median(sm), mad_sd(sm)
    top = np.percentile(sm, 99.9)
    if top - base < 15 * sd:
        return np.zeros(img.shape, dtype=int)
    fg = sm > max(base + 0.22 * (top - base), base + 10 * sd)
    fg = ndi.binary_opening(fg, iterations=1)
    log = -ndi.gaussian_laplace(img, 2.6)
    pk = peak_local_max(log, min_distance=5, threshold_abs=0.15 * np.percentile(log[fg], 95) if fg.any() else 0,
                        labels=fg.astype(int), exclude_border=False)
    seeds = np.zeros(img.shape, dtype=int)
    seeds[tuple(pk.T)] = np.arange(1, len(pk) + 1)
    return watershed(-ndi.gaussian_filter(img, 1.8), seeds, mask=fg)


def hist_mode(a):
    lo, hi = np.percentile(a, [0.5, 85])
    width = max((hi - lo) / 120.0, 0.5)
    h, e = np.histogram(a, bins=np.arange(lo, hi + width, width))
    h = ndi.uniform_filter1d(h.astype(float), 7)
    k = int(np.argmax(h))
    return 0.5 * (e[k] + e[k + 1])


yy, xx = np.indices(dark.shape)
fy = np.fft.fftfreq(dark.shape[0])[:, None]
fx = np.fft.fftfreq(dark.shape[1])[None, :]
fr = np.hypot(fy, fx)
LOW, HIGH, NOISE = (fr > 0.015) & (fr < 0.07), (fr > 0.14) & (fr < 0.30), fr > 0.42


def sharpness(img):
    p = np.abs(np.fft.fft2(img - img.mean())) ** 2
    floor = p[NOISE].mean()
    low = p[LOW].mean() - floor
    if low < 30 * floor:
        return np.nan                              # (almost) empty field: cannot be judged, keep
    return max(p[HIGH].mean() - floor, 0.0) / low


rec = []
for w in lay["well"]:
    n1, m2 = read(w, 1), read(w, 2)
    lab = nuclei(n1)
    ids = np.unique(lab)
    ids = ids[ids > 0]
    inner, outer = expand_labels(lab, 1), expand_labels(lab, 4)
    ring = np.where(inner > 0, 0, outer)
    bg = hist_mode(m2)
    sums = ndi.sum(m2, ring, ids) if len(ids) else np.array([])
    area = ndi.sum(np.ones_like(m2), ring, ids) if len(ids) else np.array([])
    has = area > 0                                   # a nucleus enclosed by its neighbours has no ring
    vals = sums[has] / area[has] - bg
    rec.append((w, len(ids), float(np.mean(vals)) if len(vals) else np.nan, sharpness(n1)))
tab = lay.merge(pd.DataFrame(rec, columns=["well", "n_cells", "marker_per_cell", "sharp"]), on="well")
tab["qc_pass"] = np.where(tab["sharp"] < 0.35 * np.nanmedian(tab["sharp"]), 0, 1)
good = tab["qc_pass"] == 1


def theil_sen(x, y):
    i, j = np.triu_indices(len(x), 1)
    dx = x[j] - x[i]
    keep = dx != 0
    return float(np.median((y[j] - y[i])[keep] / dx[keep]))


def vehicle_expectation(col):
    v = tab[(tab["role"] == "vehicle") & good]
    y = np.log(v[col].to_numpy(dtype=float))
    c, r, e = v["c"].to_numpy(float), v["r"].to_numpy(float), v["edge"].to_numpy(bool)
    bc = br = le = 0.0
    for _ in range(6):
        res = y - bc * c - br * r
        le = np.median(res[e]) - np.median(res[~e])
        z = y - le * e
        bc = theil_sen(c, z - br * r)
        br = theil_sen(r, z - bc * c)
    a = np.median(y - le * e - bc * c - br * r)
    return np.exp(a + le * tab["edge"].to_numpy(float) + bc * tab["c"].to_numpy(float) + br * tab["r"].to_numpy(float))


tab["resp_m"] = tab["marker_per_cell"] / vehicle_expectation("marker_per_cell")
tab["resp_n"] = tab["n_cells"] / vehicle_expectation("n_cells")

LE = np.linspace(-3.0, 3.5, 196)[:, None, None, None]
HH = np.linspace(0.5, 3.0, 26)[None, :, None, None]
BB = np.linspace(0.0, 0.9, 37)[None, None, :, None]


def half_point(d, y):
    ok = np.isfinite(y)
    d, y = d[ok], y[ok]
    if len(d) < 4:
        return "inactive"
    pred = BB + (1 - BB) / (1 + (d[None, None, None, :] / 10 ** LE) ** HH)
    sse = ((pred - y[None, None, None, :]) ** 2).sum(axis=3)
    i, j, k = np.unravel_index(int(np.argmin(sse)), sse.shape)
    le, h, b = LE[i, 0, 0, 0], HH[0, j, 0, 0], BB[0, 0, k, 0]
    if b >= 0.5:
        return "inactive"
    x50 = 10 ** le * ((1 - b) / (0.5 - b) - 1) ** (1 / h)
    return f"{x50:.4g}" if x50 <= d.max() else "inactive"


lines = []
for name, s in tab[(tab["role"] == "compound") & good].groupby("compound"):
    d = s["dose_uM"].to_numpy(float)
    lines.append((name, half_point(d, s["resp_m"].to_numpy(float)), half_point(d, s["resp_n"].to_numpy(float))))
pd.DataFrame(lines, columns=["compound", "ic50_marker_uM", "gi50_count_uM"]).to_csv(
    out / "dose_response.tsv", sep="\t", index=False)
wm = tab[["well", "n_cells", "marker_per_cell", "qc_pass"]].copy()
wm["marker_per_cell"] = wm["marker_per_cell"].map(lambda v: f"{v:.2f}" if np.isfinite(v) else "NA")
wm.to_csv(out / "well_measurements.tsv", sep="\t", index=False)
