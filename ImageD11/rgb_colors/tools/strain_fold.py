"""Does a strained crystal give the same folded direction for equivalent indexings ubi and M.ubi?
Cartesian route: U from the B convention, h = U^T r, fold with orthogonal m-3m matrices.
Lattice route: uvw = UB^T r, fold with the integer m-3m matrices (sector N.uvw >= 0, no metric)."""
from __future__ import print_function, division
import itertools, numpy as np
from ImageD11 import grain

# m-3m as integer matrices (signed permutations); for a cubic lattice they are also the Cartesian ones
OPS = np.array([np.diag(s)[list(p)] for p in itertools.permutations(range(3))
                for s in itertools.product((1, -1), repeat=3)], float)
N = np.array([[1, -1, 0], [-1, 0, 1], [0, 1, 0]], float)   # MTEX m-3m sector normals (hkl)

def fold(v, ops):
    """the image of each row of v inside the sector, as a unit vector of the same basis"""
    imgs = np.einsum("kij,nj->nki", ops, v)
    inside = np.all(np.einsum("mi,nki->nkm", N, imgs) >= -1e-12, axis=2)
    k = np.argmax(inside, axis=1)
    out = imgs[np.arange(len(v)), k]
    return out / np.linalg.norm(out, axis=1)[:, None]

def angle(a, b, g=np.eye(3)):
    """angle between rows of a and b (contravariant) in metric g, deg"""
    ab = np.einsum("ni,ij,nj->n", a, g, b)
    aa = np.einsum("ni,ij,nj->n", a, g, a); bb = np.einsum("ni,ij,nj->n", b, g, b)
    return np.degrees(np.arccos(np.clip(ab / np.sqrt(aa * bb), -1, 1)))

rng = np.random.RandomState(0)
r = rng.normal(size=(20000, 3)); r /= np.linalg.norm(r, axis=1)[:, None]
for eps in (0, 1e-4, 1e-3, 1e-2):
    ubi0 = 4.0 * np.linalg.qr(rng.normal(size=(3, 3)))[0]          # a,b,c rows, random orientation
    ubi = np.dot(ubi0, (np.eye(3) + eps * rng.normal(size=(3, 3))).T)       # strain the lattice vectors
    worst_cart = worst_lat = 0
    h0 = fold(np.dot(r, grain.grain(ubi).U), OPS)                    # Cartesian: h = U^T r
    u0 = fold(np.dot(r, np.linalg.inv(ubi)), OPS)                    # lattice: uvw = UB^T r, as rows r ubi^-1
    for M in OPS[1:]:
        if np.linalg.det(M) < 0:
            M = -M                                              # keep the lattice right handed
        ubiM = np.dot(M, ubi)
        h = fold(np.dot(r, grain.grain(ubiM).U), OPS)
        u = fold(np.dot(r, np.linalg.inv(ubiM)), OPS)
        worst_cart = max(worst_cart, angle(h0, h).max())
        worst_lat = max(worst_lat, np.abs(u0 - u).max())
    print("strain %.0e: Cartesian fold differs by up to %.2e deg; lattice fold uvw differs by up to %.1e"
          % (eps, worst_cart, worst_lat))

print("angle from the folded direction to the sector vertex [111], over the equivalent indexings:")
def gsym(g):
    """the metric averaged over the group: invariant, so every operation is an isometry of it"""
    return np.mean([M.dot(g).dot(M.T) for M in OPS], axis=0)
for eps in (1e-4, 1e-3, 1e-2):
    ubi = 4.0 * np.dot(np.linalg.qr(rng.normal(size=(3, 3)))[0], (np.eye(3) + eps * rng.normal(size=(3, 3))).T)
    v = np.array([[1.0, 1, 1]])
    spread_g, spread_s = [], []
    for M in OPS:
        if np.linalg.det(M) < 0:
            M = -M
        ubiM = np.dot(M, ubi)
        u = fold(np.dot(r[:2000], np.linalg.inv(ubiM)), OPS)
        g = np.dot(ubiM, ubiM.T)                   # g_ij = a_i . a_j, for uvw components
        spread_g.append(angle(u, np.repeat(v, len(u), 0), g))
        spread_s.append(angle(u, np.repeat(v, len(u), 0), gsym(g)))
    print("strain %.0e: own metric spread %.2e deg, group-averaged metric spread %.1e deg" % (
        eps, np.ptp(spread_g, axis=0).max(), np.ptp(spread_s, axis=0).max()))
