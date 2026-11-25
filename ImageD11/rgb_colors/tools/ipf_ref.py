"""Reference implementation of ipf_colour_math.md: MTEX 7.2 TSL and HSV direction keys from keys.json.
Every vector is a lattice vector (uvw) or plane normal (hkl); angles use the symmetrised metric g_s.
For the arithmetic the lattice is mapped once to orthonormal coordinates x = L^T uvw, g_s = L L^T;
only dot and cross products enter, so any right handed factor of g_s gives the same colours."""
import json, numpy as np
from ImageD11 import sym_u

def cell_metric(cell):
    """g_ij = a_i . a_j of a cell a b c alpha beta gamma (deg)"""
    a, b, c = cell[:3]; al, be, ga = np.radians(cell[3:])
    return np.array([[a * a, a * b * np.cos(ga), a * c * np.cos(be)],
                     [a * b * np.cos(ga), b * b, b * c * np.cos(al)],
                     [a * c * np.cos(be), b * c * np.cos(al), c * c]])

def group_uvw(gens):
    """the group as matrices acting on uvw columns (ITA triplets act on uvw as written)"""
    return np.array([o.T for o in sym_u.generate_group(*gens).group])

def symmetrise(g, M):
    """g_s = <M^T g M> for M acting on uvw (equivalently <O g O^T> for O = M^-T acting on ubi)"""
    return np.mean([m.T.dot(g).dot(m) for m in M], axis=0)

class Frame:
    def __init__(self, gs):
        self.L = np.linalg.cholesky(gs)              # gs = L L^T
    def vec(self, t):
        """tagged lattice vector {"uvw": ..} or {"hkl": ..} -> unit orthonormal coordinates"""
        if "uvw" in t:
            x = self.L.T.dot(np.asarray(t["uvw"], float))
        else:
            x = np.linalg.solve(self.L, np.asarray(t["hkl"], float))
        return x / np.linalg.norm(x)
    def uvw(self, u):
        x = np.asarray(u, float).dot(self.L)
        return x / np.linalg.norm(x, axis=-1)[..., None]
    def ops(self, M):
        Li = np.linalg.inv(self.L.T)
        return np.array([self.L.T.dot(m).dot(Li) for m in M])

unit = lambda x: x / np.linalg.norm(x, axis=-1)[..., None]
ang = lambda a, b: np.arccos(np.clip(np.sum(a * b, -1), -1, 1))

def fold_nearest(x, R, c):
    """MTEX project2FundamentalRegion: the image nearest the centre c"""
    imgs = np.einsum("kij,nj->nki", R, x)
    k = np.argmax(np.round(np.dot(imgs, c), 12), axis=1)
    return imgs[np.arange(len(x)), k]

def fold(x, R, N):
    """the image R x inside the sector N.x >= 0 (first one; on a boundary any image is equivalent)"""
    imgs = np.einsum("kij,nj->nki", R, x)
    s = np.einsum("mi,nki->nkm", N, imgs).min(axis=2) if len(N) else np.zeros(imgs.shape[:2])
    k = np.argmax(s >= -1e-9, axis=1)
    return imgs[np.arange(len(x)), k]

def vertices(N):
    """corners of the region N.x >= 0: the pairwise plane intersections that lie inside"""
    V = []
    for i in range(len(N)):
        for j in range(i + 1, len(N)):
            c = np.cross(N[i], N[j])
            if np.linalg.norm(c) < 1e-12:
                continue
            for v in (unit(c), -unit(c)):
                if np.all(np.dot(N, v) >= -1e-9) and not any(np.allclose(v, w, atol=1e-9) for w in V):
                    V.append(v)
    return np.array(V).reshape(-1, 3)

def area_centroid(N, V):
    """unit vector of int_S x dA: half the sum over edges of arc length times inward edge normal"""
    if len(V) < 3:
        raise ValueError
    s = np.zeros(3)
    for n in N:
        on = [v for v in V if abs(np.dot(n, v)) < 1e-9]
        if len(on) == 2:
            s += 0.5 * ang(on[0], on[1]) * n
    return unit(s)

def center(N, V, grid_like_mtex=False):
    """MTEX sphericalRegion.center"""
    if len(N) < 2:
        return unit(N[0])
    if len(V) < 3:   # the two normals furthest apart, averaged
        w = np.arccos(np.clip(np.dot(N, N.T), -1, 1)); i, j = np.unravel_index(np.argmax(w), w.shape)
        return unit(N[i] + N[j])
    if len(V) == 3:
        return unit(V.sum(0))
    return area_centroid(N, V)    # MTEX: mean of a 1 deg grid over the sector

def calc_angle(c, rx, v):
    rx = unit(rx - np.dot(rx, c) * c); ry = unit(np.cross(c, rx))
    dv = v - c
    rho = np.mod(np.arctan2(np.dot(dv, ry), np.dot(dv, rx)), 2 * np.pi)
    return np.where(np.isfinite(rho), rho, 0.0)

def correct_angle(N, V, c, rx, n=1000):
    """MTEX polarCoordinates/correctAngle: hue proportional to the integrated distance to the boundary;
    with 3 vertices each of the 3 edges gets a third of the hue circle"""
    rho = np.linspace(0, 2 * np.pi, n)
    if len(V) == 3:
        rho_v = calc_angle(c, rx, V)
        rho = np.unique(np.concatenate([rho, rho_v]))
    rxp = unit(rx - np.dot(rx, c) * c)
    n0 = unit(np.cross(rxp, c))
    # rotate n0 about c by rho (right handed)
    nn = (np.outer(np.cos(rho), n0) + np.outer(np.sin(rho), np.cross(c, n0)) + np.outer(1 - np.cos(rho), c * np.dot(c, n0)))
    om = np.min([ang(unit(np.cross(Ni, nn)), c) for Ni in N], axis=0)
    d = (om[:-1] + om[1:]) / 2 * np.diff(rho)
    if len(V) == 3:
        arc = np.mod(np.cumsum(np.isin(rho[:-1], rho_v)), 3)
        s = np.bincount(arc, d, minlength=3)
        d = d / s[arc] / 3
    return rho, 2 * np.pi * np.concatenate([[0], np.cumsum(d / d.sum())])

def polar(N, V, x, c, ref, oos):
    """MTEX sphericalRegion.polarCoordinates: r = 1 at the centre, 0 on the boundary; rho the hue angle"""
    if np.all(np.abs(np.dot(N, c)) <= 1e-10):
        r = ang(c, x) / np.pi
    else:
        vxc = unit(np.cross(x, c))
        r = np.full(len(x), np.inf)
        for Ni in N:
            bc = unit(np.cross(vxc, Ni))
            with np.errstate(invalid="ignore", divide="ignore"):
                d = ang(-x, bc) / ang(-c, bc)
            r = np.minimum(r, np.where(np.isnan(d), 1.0, d))
    tnull = lambda rx: np.linalg.norm(rx - np.dot(rx, c) * c) < 1e-10
    rx = oos - c
    if tnull(rx):
        rx = ref - c
    if tnull(rx):
        raise ValueError("degenerate: centre, out-of-screen and hue reference coincide (MTEX falls back to orth())")
    rho = calc_angle(c, rx, x)
    atv = [i for i, v in enumerate(V) if ang(c, v) < 1e-10]
    if atv:
        raise NotImplementedError("centre in a vertex")
    if np.any(np.abs(np.dot(N, c)) <= 1e-10):
        rho = np.mod(rho, np.pi) * 2
    elif len(V):
        irho, om = correct_angle(N, V, c, rx)
        rho = np.interp(rho, irho, om)
    return r, rho

def hsl2rgb(h, s, l):
    """MTEX hsl2hsv then MATLAB hsv2rgb"""
    ll = 2 * l; ss = np.where(ll <= 1, s * ll, s * (2 - ll))
    v = (ll + ss) / 2
    with np.errstate(invalid="ignore", divide="ignore"):
        sv = np.nan_to_num(2 * ss / (ll + ss))
    return hsv2rgb(h, sv, v)

def hsv2rgb(h, s, v):
    h = 6 * np.mod(h, 1); k = np.floor(h).astype(int) % 6; f = h - np.floor(h)
    p, q, t = v * (1 - s), v * (1 - s * f), v * (1 - s * (1 - f))
    return np.choose(k, [np.stack(c, -1).T for c in
                         [(v, t, p), (q, v, p), (p, v, t), (p, q, v), (t, p, v), (v, p, q)]]).T

MTEX_F = False   # True: MTEX's cumsum(f) / sum(f), which starts at F(0) = 1.7e-3

def hue_correction(h):
    z = np.linspace(0, 1, 1000)
    f = 0.5 + sum(np.exp(-200 * (np.mod(z - r + 0.5, 1) - 0.5) ** 2) for r in (0, 0.6666, 0.3333))
    if MTEX_F:
        return np.interp(h, z, np.cumsum(f / f.sum()))
    f = np.concatenate([[0], np.cumsum(f[1:] + f[:-1])])
    return np.interp(h, z, f / f[-1])

def ar2rgb(omega, radius, gray, hue_corr):
    L = (radius - 0.5) * gray + 0.5
    with np.errstate(invalid="ignore", divide="ignore"):
        S = np.nan_to_num(gray * (1 - np.abs(2 * radius - 1)) / (1 - np.abs(2 * L - 1)))
    h = hue_correction(omega) if hue_corr else omega
    return hsl2rgb(h, S, L)

def tsl(key, F, R, x):
    N = np.array([F.vec(t) for t in key["normals"]]).reshape(-1, 3)
    V = vertices(N)
    x = fold(x, R, N)
    c = center(N, V)
    r, rho = polar(N, V, x, c, F.vec(key["hue_ref"]), F.vec(key["out_of_screen"]))
    return ar2rgb(np.mod(rho / 2 / np.pi, 1), 0.5 + r / 2, 1.0, False)

def hsv(key, F, R, x, gray_value=(0.2, 0.5), gray_gradient=0.5):
    N0 = np.array([F.vec(t) for t in key["normals"]]).reshape(-1, 3)
    S = np.array([F.vec(t) for t in key["split"]]).reshape(-1, 3)
    x = fold(x, R, N0)
    if any(abs(np.dot(s, n) + 1) < 1e-9 for s in S for n in N0):
        raise ValueError("degenerate: a split plane is a boundary plane turned inside out, the region is empty")
    N = np.concatenate([N0, S]); V = vertices(N)
    wc = center(N, V)
    white = np.ones(len(x), bool); switch = False
    for s in S:   # the split planes: the far side is reflected in and coloured towards black
        ind = np.dot(x, s) < 1e-5
        x[ind] = x[ind] - 2 * np.outer(np.dot(x[ind], s), s)
        white &= ~ind
        if np.dot(wc, s) < 1e-5:
            wc = wc - 2 * np.dot(wc, s) * s; switch = not switch
    white ^= switch
    r, rho = polar(N, V, x, wc, F.vec(key["hue_ref"]), F.vec(key["out_of_screen"]))
    radius = np.where(white, 0.5 + r / 2, (1 - r) / 2)
    th = radius * np.pi
    hi = th > np.pi / 2
    th = (2 * gray_gradient * th / np.pi + (1 - gray_gradient) * (1 - np.cos(th))) / 2
    gray = np.where(hi, 1 - 2 * gray_value[0] * np.abs(th - 0.5), 1 - 2 * gray_value[1] * np.abs(th - 0.5))
    return ar2rgb(np.mod(rho / 2 / np.pi, 1), th, gray, True)

def colours(group, kind, key, cell, uvw, gs=None):
    M = group_uvw(group["generators"])
    if kind == "tsl":
        M = np.concatenate([M, -M])
    g = cell_metric(cell)
    F = Frame(symmetrise(g, M) if gs is None else gs)
    x = F.uvw(uvw)
    return (tsl if kind == "tsl" else hsv)(key, F, F.ops(M), x)
