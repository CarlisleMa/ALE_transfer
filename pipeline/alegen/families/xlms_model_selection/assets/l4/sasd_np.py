"""Surface distance on a lattice, numpy only (bucketed shortest-path search, Dial's algorithm).

Used by the generator (truth) and by the independent second solution. The reference solution
(xl4_pipeline.py) is a separate implementation on scipy (k-d tree + sparse-graph Dijkstra).

Definition implemented here (see xlms_notes.md of a level-4 task):
  * lattice points = points with integer x, y, z (A) in the frame of the model file;
  * a lattice point is blocked if it is closer than `radius` to any heavy atom, otherwise open;
  * steps join open points that differ by at most 1 in every coordinate, step length = their
    straight-line distance (1, sqrt 2, sqrt 3);
  * each end is joined to the lattice by a straight tether from the CA atom to an open lattice
    point at most `tether` away (not checked against the protein);
  * distance = minimum of tether + path + tether; infinite if there is none.
"""
import numpy as np


class Lattice:
    def __init__(self, atoms, radius):
        atoms = np.asarray(atoms, dtype=np.float64)
        pad = int(np.ceil(radius)) + 3            # >= 2 open layers around the excluded volume + 1 border layer
        self.lo = np.floor(atoms.min(0)).astype(np.int64) - pad
        hi = np.ceil(atoms.max(0)).astype(np.int64) + pad
        self.shape = tuple(int(v) for v in (hi - self.lo + 1))
        blocked = np.zeros(self.shape, dtype=bool)
        blocked[0], blocked[-1] = True, True
        blocked[:, 0], blocked[:, -1] = True, True
        blocked[:, :, 0], blocked[:, :, -1] = True, True
        base = np.rint(atoms).astype(np.int64)
        frac = atoms - base
        k = int(np.floor(radius + 0.5)) + 1
        r2 = radius * radius
        rel = base - self.lo
        for dx in range(-k, k + 1):
            for dy in range(-k, k + 1):
                for dz in range(-k, k + 1):
                    if dx * dx + dy * dy + dz * dz > (radius + 0.87) ** 2:
                        continue
                    d2 = (frac[:, 0] - dx) ** 2 + (frac[:, 1] - dy) ** 2 + (frac[:, 2] - dz) ** 2
                    m = d2 < r2
                    if m.any():
                        blocked[rel[m, 0] + dx, rel[m, 1] + dy, rel[m, 2] + dz] = True
        self.open = ~blocked.ravel()
        sy, sz = self.shape[1], self.shape[2]
        offs, wts = [], []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    if dx or dy or dz:
                        offs.append((dx * sy + dy) * sz + dz)
                        wts.append(float(np.sqrt(dx * dx + dy * dy + dz * dz)))
        self.offs, self.wts = offs, wts

    def tether_points(self, anchor, tether):
        """Open lattice points within `tether` of the anchor: (flat indices, straight-line distances)."""
        anchor = np.asarray(anchor, dtype=np.float64)
        t = int(np.ceil(tether)) + 1
        c = np.rint(anchor).astype(np.int64)
        r = np.arange(-t, t + 1)
        gx, gy, gz = np.meshgrid(c[0] + r, c[1] + r, c[2] + r, indexing="ij")
        pts = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], axis=1)
        d = np.sqrt(((pts - anchor) ** 2).sum(1))
        idx3 = pts - self.lo
        ok = (d <= tether) & ((idx3 >= 0) & (idx3 < np.array(self.shape))).all(1)
        flat = np.ravel_multi_index(idx3[ok].T, self.shape)
        d = d[ok]
        op = self.open[flat]
        return flat[op], d[op]

    def field(self, start_idx, start_d, limit):
        """Shortest tether+path length from one anchor to every lattice point (inf beyond `limit`)."""
        dist = np.full(self.open.size, np.inf)
        if not len(start_idx):
            return dist
        dist[start_idx] = start_d
        buckets = {}
        kb = np.floor(start_d).astype(np.int64)
        for b in np.unique(kb):
            buckets.setdefault(int(b), []).append(start_idx[kb == b])
        b = 0
        top = int(np.floor(limit))
        opn = self.open
        while b <= top:
            if b in buckets:
                cand = np.unique(np.concatenate(buckets.pop(b)))
                dc = dist[cand]
                m = (dc >= b) & (dc < b + 1)
                fr, df = cand[m], dc[m]
                if fr.size:
                    for o, w in zip(self.offs, self.wts):
                        nb = fr + o
                        nd = df + w
                        m2 = opn[nb] & (nd < dist[nb]) & (nd <= limit)
                        if m2.any():
                            nbm, ndm = nb[m2], nd[m2]
                            dist[nbm] = ndm
                            hi = ndm >= b + 2
                            if hi.any():
                                buckets.setdefault(b + 2, []).append(nbm[hi])
                                nbm = nbm[~hi]
                            if nbm.size:
                                buckets.setdefault(b + 1, []).append(nbm)
            b += 1
        return dist


def pairwise(lat, anchors, tether, limit):
    """Matrix of surface distances between anchors (inf where > limit or unreachable)."""
    tp = [lat.tether_points(a, tether) for a in anchors]
    n = len(anchors)
    out = np.full((n, n), np.inf)
    for i in range(n):
        if not len(tp[i][0]):
            continue
        dist = lat.field(tp[i][0], tp[i][1], limit)
        for j in range(n):
            if j != i and len(tp[j][0]):
                v = float((dist[tp[j][0]] + tp[j][1]).min())
                if v <= limit:
                    out[i, j] = v
    return out
