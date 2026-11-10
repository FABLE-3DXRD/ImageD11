"""Inverse pole figure (IPF) colours of orientation matrices.

The colour of a grain is the colour of the crystal direction along a sample
axis, uvw = inv(ubi).T . axis (ubi rows are a, b, c), in an MTEX 7.2 colour
key (TSL or HSV) of the crystal's point group.  The keys are tabulated as
integer plane normals (hkl) and directions (uvw) in ipf_keys.py, so nothing
depends on a B matrix convention or on an origin of orientation space:

* the direction is folded into the sector h.uvw >= 0 with the integer
  symmetry operations (TSL: the Laue group, HSV: the point group), so all
  equivalent ubi give one colour.  h is the hkl of a boundary plane, or
  g_s.uvw for a plane normal to a direction (monoclinic HSV split planes);
* angles use the symmetrised metric g_s = <M^T g M> over the group
  (g = ubi.ubi^T), which every equivalent ubi shares.  The strain that g_s
  removes is reported by symmetric_metric; colours are warned about when it
  exceeds maxstrain.

For the arithmetic the metric is factored once per grain, g_s = L.L^T, so a
direction is L^T.uvw and a plane normal inv(L).hkl.  Only dot and cross
products enter, so any right handed factor gives the same colour.

MTEX 7.2.0 TSLDirectionKey and HSVDirectionKey (default options), with two
of MTEX's approximations done exactly: the hue integral (MTEX sums 1000 steps,
up to 6e-5 in rgb) and the centre of a sector with 4 corners (TSL m-3 and 23;
MTEX averages a 1 deg grid, up to 6e-3).  The HSV hue correction starts at 0,
where MTEX's starts at 1.7e-3 and jumps (up to 1e-2 in rgb).
See ipf_colour_math.md.
"""

from __future__ import print_function, division

import logging
import numpy as np

from ImageD11 import sym_u
from ImageD11.rgb_colors.ipf_keys import GROUPS


def point_group(name):
    """The sym_u group of a point group name in ipf_keys.GROUPS (e.g. "6mm", "121", "-3mR")"""
    return sym_u.generate_group(*GROUPS[name]["generators"])


def _opset(ops):
    return frozenset(tuple(int(x) for x in np.round(o).ravel()) for o in ops)


_NAMES = {}


def group_name(grp):
    """The ipf_keys.GROUPS name of a sym_u group, found from its operations"""
    return _name_of(grp.group)


def _name_of(ops):
    if not _NAMES:
        for name in GROUPS:
            _NAMES[_opset(point_group(name).group)] = name
    s = _opset(ops)
    if s not in _NAMES:
        raise ValueError("No colour key for this group (%d operations)" % len(ops))
    return _NAMES[s]


def spacegroup_group(spacegroup, cell_choice="standard"):
    """The ipf_keys.GROUPS name of the point group of a space group (1-230).
    cell_choice "rhombohedral" puts R groups on rhombohedral axes (as unitcell)"""
    from xfab import sg
    spg = sg.sg(sgno=spacegroup)
    # xfab rotations act on uvw columns, sym_u holds their transposes
    name = _name_of([np.transpose(r) for r in spg.rot])
    if cell_choice == "rhombohedral" and spg.name.startswith("R"):
        name = {"3": "3R", "-3": "-3R", "321": "32R", "3m1": "3mR", "-3m1": "-3mR"}[name]
    return name


def _uvw_ops(grp):
    """The operations as matrices M acting on uvw columns (sym_u holds M^T, acting on ubi)"""
    return np.array([np.transpose(o) for o in grp.group])


def symmetric_metric(ubis, grp):
    """
    g_s = <M^T g M> over the group, g = ubi.ubi^T, for ubis (n,3,3).
    Returns g_s (n,3,3) and the strain that g_s removes,
    e = (eigenvalues(inv(g_s).g) - 1)/2, (n,3).  sum(e) = 0.
    """
    ubis = np.asarray(ubis, float).reshape(-1, 3, 3)
    M = _uvw_ops(grp)
    g = np.einsum("nij,nkj->nik", ubis, ubis)
    # g_s is linear in g: (g_s)_il = sum_jk P_iljk g_jk
    P = np.einsum("mji,mkl->iljk", M, M).reshape(9, 9) / len(M)
    gs = np.dot(g.reshape(-1, 9), P.T).reshape(-1, 3, 3)
    e = (np.sort(np.linalg.eigvals(np.linalg.solve(gs, g)).real, axis=-1) - 1) / 2
    return gs, e


class _Key(object):
    """One colour key of one group"""

    def __init__(self, name, kind, view):
        grp = GROUPS[name]
        keys = [k for k in grp[kind] if view in k["views"]]
        if not keys:
            raise ValueError("No %s key for %s with view %s, try one of %s" % (
                kind, name, view, sorted(set(v for k in grp[kind] for v in k["views"]))))
        k = keys[0]
        self.kind = kind
        M = _uvw_ops(point_group(name)).round().astype(int)
        self.ops = np.concatenate([M, -M]) if kind == "tsl" else M
        self.normals = k["normals"]
        self.planes = tuple(k["normals"]) + tuple(k.get("split", ()))
        self.oos = k["out_of_screen"]
        self.hue_ref = k["hue_ref"]
        self.hue_table = None

    def fold(self, uvw, gs):
        """the image M.uvw inside the sector, (n,3): h.uvw >= 0 for each
        boundary plane, h = hkl, or g_s.uvw for a plane normal to a direction"""
        if len(self.normals) == 0:
            return uvw
        n = len(uvw)
        if all(b == "hkl" for b, _ in self.normals):
            h = np.array([v for _, v in self.normals], float)
            # h.(M.uvw) = (M^T.h).uvw, every op and plane in one product
            hm = np.einsum("mi,kij->jkm", h, self.ops).reshape(3, -1)
            s = np.dot(uvw, hm).reshape(n, len(self.ops), len(h)).min(axis=2)
            hlen = np.sqrt((h * h).sum(-1)).max()
        else:
            h = np.stack([np.einsum("nij,j->ni", gs, v) if b == "uvw" else
                          np.tile(np.asarray(v, float), (n, 1)) for b, v in self.normals], axis=1)
            s = np.einsum("nmi,kij,nj->nkm", h, self.ops, uvw).min(axis=2)
            hlen = np.sqrt((h * h).sum(-1)).max(1)
        tol = 1e-9 * hlen * np.sqrt((uvw * uvw).sum(-1))
        k = np.argmax(s >= -tol[:, None], axis=1)
        return np.einsum("nij,nj->ni", self.ops[k], uvw)


def _corners(N):
    """corners of the region N.x >= 0 as (i, j, sign), x = sign * N_i x N_j:
    the plane intersections inside it, found on the first grain; and the
    edges (plane, corner, corner): planes with two corners on them"""
    n = N[0]
    found, V = [], []
    for i in range(len(n)):
        for j in range(i + 1, len(n)):
            c = np.cross(n[i], n[j])
            if np.dot(c, c) < 1e-20:
                continue
            c = c / np.sqrt(np.dot(c, c))
            for sign in (1, -1):
                v = sign * c
                if (np.dot(n, v) >= -1e-9).all() and not any(np.allclose(v, w, atol=1e-9) for w in V):
                    V.append(v)
                    found.append((i, j, sign))
    edges = []
    for i in range(len(n)):
        on = [a for a, v in enumerate(V) if abs(np.dot(n[i], v)) < 1e-9]
        if len(on) == 2:
            edges.append((i, on[0], on[1]))
    return found, edges


# ---------------------------------------------------------------------------
# geometry on the unit sphere, one frame per grain: arrays (n, ..., 3)
# ---------------------------------------------------------------------------

def _dot(a, b):
    return (a * b).sum(-1)


def _unit(x):
    with np.errstate(invalid="ignore", divide="ignore"):
        return x / np.sqrt(_dot(x, x))[..., None]


def _ang(a, b):
    """angle between unit vectors, accurate near 0 and pi (arccos is not)"""
    c = np.cross(a, b)
    return np.arctan2(np.sqrt(_dot(c, c)), _dot(a, b))


def _vec(L, Li, t):
    """tagged lattice vector -> unit vectors (n,3) in each grain's frame"""
    basis, v = t
    v = np.asarray(v, float)
    if basis == "uvw":
        return _unit(np.einsum("nji,j->ni", L, v))
    return _unit(np.einsum("nij,j->ni", Li, v))


def _centre(N, V, edges):
    """MTEX sphericalRegion.center"""
    if len(N[0]) == 1:
        return N[:, 0]
    if V is None or V.shape[1] < 3:
        # lune: the two normals furthest apart
        w = _ang(N[:, :, None], N[:, None, :]).reshape(len(N), -1)
        i, j = np.divmod(np.argmax(w, axis=1), N.shape[1])
        n = np.arange(len(N))
        return _unit(N[n, i] + N[n, j])
    if V.shape[1] == 3:
        return _unit(V.sum(axis=1))
    # area centroid: int x dA = sum over edges of arc length x inward normal / 2
    s = 0
    for i, a, b in edges:
        s = s + _ang(V[:, a], V[:, b])[:, None] * N[:, i]
    return _unit(s)


def _azimuth(c, rx, v):
    """angle of v about the centre c from rx, in [0, 2 pi)"""
    rx = _unit(rx - _dot(rx, c)[..., None] * c)
    ry = _unit(np.cross(c, rx))
    dv = v - c
    rho = np.mod(np.arctan2(_dot(dv, ry), _dot(dv, rx)), 2 * np.pi)
    return np.where(np.isfinite(rho), rho, 0.0)


_GX, _GW = np.polynomial.legendre.leggauss(32)
_PX, _PW = np.polynomial.legendre.leggauss(8)


def _exit(N, c, n0, rho):
    """distance from the centre c to the boundary along azimuth rho (n,k):
    on the great circle through c with normal nn = n0 turned about c by rho,
    the angle from c to N x nn is atan2(|N.c|, N.(nn x c)), and
    nn x c = cos(rho) n0 x c + sin(rho) n0"""
    a = np.abs(np.einsum("nmi,ni->nm", N, c))[:, None]
    A = np.einsum("nmi,ni->nm", N, n0)[:, None]
    B = np.einsum("nmi,ni->nm", N, np.cross(n0, c))[:, None]
    return np.arctan2(a, B * np.cos(rho)[..., None] + A * np.sin(rho)[..., None]).min(axis=-1)


def _hue_table(N, V, c, rx, pieces=64):
    """MTEX polarCoordinates/correctAngle: the hue grows with the integral over
    the azimuth of the distance from the centre to the boundary; with 3 corners
    each edge gets a third of the hue circle.  The distance is smooth between
    corners, so the integral is Gauss-Legendre on pieces of the arcs between
    them (MTEX sums 1000 steps).
    The sector is the same figure in every metric of the group (it is made of
    lattice planes and directions whose angles the group fixes; tested), so
    the table comes from the first grain and serves all: knots (k,), hue at
    the knots (k,), dhue/dintegral on each piece (k-1,) and the figure."""
    N, V, c, rx = N[:1], V[:1], c[:1], rx[:1]
    rxp = _unit(rx - _dot(rx, c)[:, None] * c)
    n0 = _unit(np.cross(rxp, c))
    b = np.concatenate([[0], np.sort(_azimuth(c[:, None], rx[:, None], V)[0]), [2 * np.pi]])
    knots = np.concatenate([np.linspace(lo, hi, pieces + 1)[:-1] for lo, hi in zip(b[:-1], b[1:])] +
                           [[2 * np.pi]])
    arc = np.repeat(np.arange(len(b) - 1), pieces)
    half = np.diff(knots) / 2
    x = (knots[:-1] + half)[:, None] + half[:, None] * _GX
    piece = np.dot(_exit(N, c, n0, x.reshape(1, -1)).reshape(x.shape), _GW) * half
    if V.shape[1] == 3:
        # arcs 0 and 3 are one edge, across rho = 0
        label = arc % 3
        w = 1 / (3 * np.bincount(label, piece)[label])
    else:
        w = np.full(len(piece), 1 / piece.sum())
    hue = 2 * np.pi * np.concatenate([[0], np.cumsum(piece * w)])
    return knots, hue, 2 * np.pi * w, (N, c, n0)


def _even_hue(table, rho):
    """the hue at azimuths rho (n,): the table up to the knot below, then
    Gauss on the rest of that piece"""
    knots, hue, w, (N, c, n0) = table
    p = np.clip(np.searchsorted(knots, rho, side="right") - 1, 0, len(w) - 1)
    half = (rho - knots[p]) / 2
    x = (knots[p] + half)[:, None] + half[:, None] * _PX
    f = _exit(N, c, n0, x.reshape(1, -1)).reshape(x.shape)
    return hue[p] + w[p] * np.dot(f, _PW) * half


def _polar(key, N, V, x, c, ref, oos):
    """MTEX sphericalRegion.polarCoordinates: radius 1 at the centre, 0 on the
    boundary, and the hue angle"""
    vxc = _unit(np.cross(x, c))
    r = np.inf
    with np.errstate(invalid="ignore", divide="ignore"):
        for i in range(N.shape[1]):
            bc = _unit(np.cross(vxc, N[:, i]))
            d = _ang(-x, bc) / _ang(-c, bc)
            r = np.minimum(r, np.where(np.isnan(d), 1.0, d))
    rx = oos - c
    flat = np.sqrt(_dot(rx, rx) - _dot(rx, c) ** 2) < 1e-10
    rx[flat] = (ref - c)[flat]
    rho = _azimuth(c, rx, x)
    if V is not None:
        if key.hue_table is None:
            key.hue_table = _hue_table(N, V, c, rx)
        rho = _even_hue(key.hue_table, rho)
    return r, rho


# ---------------------------------------------------------------------------
# angle and radius -> rgb (MTEX ar2rgb, hsl2hsv, MATLAB hsv2rgb)
# ---------------------------------------------------------------------------

def _hsv2rgb(h, s, v):
    h = 6 * np.mod(h, 1)
    k = np.floor(h).astype(int) % 6
    f = h - np.floor(h)
    p, q, t = v * (1 - s), v * (1 - s * f), v * (1 - s * (1 - f))
    r = np.choose(k, [v, q, p, p, t, v])
    g = np.choose(k, [t, v, v, q, p, p])
    b = np.choose(k, [p, p, t, v, v, q])
    return np.stack([r, g, b], axis=-1)


def _hsl2rgb(h, s, l):
    ll = 2 * l
    ss = np.where(ll <= 1, s * ll, s * (2 - ll))
    v = (ll + ss) / 2
    with np.errstate(invalid="ignore", divide="ignore"):
        sv = np.nan_to_num(2 * ss / (ll + ss))
    return _hsv2rgb(h, sv, v)


_Z = np.linspace(0, 1, 1000)
_F = 0.5 + sum(np.exp(-200 * (np.mod(_Z - r + 0.5, 1) - 0.5) ** 2) for r in (0, 0.6666, 0.3333))
_F = np.concatenate([[0], np.cumsum(_F[1:] + _F[:-1])])   # trapezoid: F(0) = 0, F(1) = 1
_F /= _F[-1]


def _ar2rgb(omega, radius, gray, hue_correction):
    L = (radius - 0.5) * gray + 0.5
    with np.errstate(invalid="ignore", divide="ignore"):
        S = np.nan_to_num(gray * (1 - np.abs(2 * radius - 1)) / (1 - np.abs(2 * L - 1)))
    h = np.interp(omega, _Z, _F) if hue_correction else omega
    return _hsl2rgb(h, S, L)


def _region(L, Li, planes):
    """unit plane normals (n,m,3), corners (n,v,3) or None, edges"""
    N = np.stack([_vec(L, Li, t) for t in planes], axis=1)
    corners, edges = _corners(N)
    if not corners:
        return N, None, edges
    V = np.stack([_unit(sign * np.cross(N[:, i], N[:, j])) for i, j, sign in corners], axis=1)
    return N, V, edges


def _tsl(key, L, Li, x):
    N, V, edges = _region(L, Li, key.normals)
    c = _centre(N, V, edges)
    r, rho = _polar(key, N, V, x, c, _vec(L, Li, key.hue_ref), _vec(L, Li, key.oos))
    return _ar2rgb(np.mod(rho / 2 / np.pi, 1), 0.5 + r / 2, 1.0, False)


def _hsv(key, L, Li, x, gray_value=(0.2, 0.5), gray_gradient=0.5):
    N, V, edges = _region(L, Li, key.planes)
    c = _centre(N, V, edges)
    white = np.ones(len(x), bool)
    switch = np.zeros(len(x), bool)
    for s in N[:, len(key.normals):].transpose(1, 0, 2):
        # beyond a split plane: reflected in, coloured towards black
        ind = _dot(x, s) < 1e-5
        x = np.where(ind[:, None], x - 2 * _dot(x, s)[:, None] * s, x)
        white &= ~ind
        flip = _dot(c, s) < 1e-5
        c = np.where(flip[:, None], c - 2 * _dot(c, s)[:, None] * s, c)
        switch ^= flip
    white ^= switch
    r, rho = _polar(key, N, V, x, c, _vec(L, Li, key.hue_ref), _vec(L, Li, key.oos))
    radius = np.where(white, 0.5 + r / 2, (1 - r) / 2)
    th = radius * np.pi
    hi = th > np.pi / 2
    th = (2 * gray_gradient * th / np.pi + (1 - gray_gradient) * (1 - np.cos(th))) / 2
    gray = np.where(hi, 1 - 2 * gray_value[0] * np.abs(th - 0.5), 1 - 2 * gray_value[1] * np.abs(th - 0.5))
    return _ar2rgb(np.mod(rho / 2 / np.pi, 1), th, gray, True)


def ubi_to_rgb(ubis, group, view="default", key="tsl", axis=(0.0, 0.0, 1.0),
               maxstrain=5e-3, chunk=4096):
    """
    IPF colours (n,3) of orientation matrices ubis (n,3,3), rows a, b, c.

    group : a sym_u group, or a name in ipf_keys.GROUPS.  Its operations
            choose the key, so give the point group: the TSL key folds with its
            Laue group, the HSV key with the point group itself.
    view  : MTEX how2plot, "default" or "out-of-screen|east" such as "+c*|+a"
    key   : "tsl" or "hsv"
    axis  : the sample direction, colours show the crystal direction along it
    maxstrain : a warning is logged when a ubi departs from the symmetry of the
            group by more than this strain (see symmetric_metric)
    """
    name = group_name(group) if hasattr(group, "group") else group
    k = _Key(name, key, view)
    ubis = np.asarray(ubis, float).reshape(-1, 3, 3)
    gs, e = symmetric_metric(ubis, point_group(name))
    bad = np.abs(e).max(axis=1) > maxstrain
    if bad.any():
        logging.warning("%d of %d ubi depart from %s symmetry by a strain up to %.2g > %.2g: "
                        "their colours use the symmetrised metric" % (
                            bad.sum(), len(ubis), name, np.abs(e).max(), maxstrain))
    uvw = np.einsum("nji,j->ni", np.linalg.inv(ubis), np.asarray(axis, float))
    return _uvw_to_rgb(k, uvw, gs, chunk)


def _uvw_to_rgb(k, uvw, gs, chunk=4096):
    """colours (n,3) of directions uvw (n,3) with metrics gs (n,3,3)"""
    rgb = np.empty((len(uvw), 3))
    for i in range(0, len(uvw), chunk):
        sl = slice(i, i + chunk)
        L = np.linalg.cholesky(gs[sl])
        Li = np.linalg.inv(L)
        x = _unit(np.einsum("nji,nj->ni", L, k.fold(uvw[sl], gs[sl])))
        rgb[sl] = (_tsl if k.kind == "tsl" else _hsv)(k, L, Li, x)
    return rgb


def _view_axis(token, T=None):
    """"+c*" -> ("hkl", (0, 0, 1)), "-a" -> ("uvw", (-1, 0, 0)).
    T: rows a_H, b_H, c_H of a group on rhombohedral axes, whose views name
    the hexagonal axes: uvw_R = T^T.uvw_H, hkl_R = inv(T).hkl_H"""
    v = np.zeros(3)
    v["abc".index(token[1])] = -1 if token[0] == "-" else 1
    star = token.endswith("*")
    if T is not None:
        T = np.array(T, float)
        v = np.linalg.solve(T, v) if star else np.dot(T.T, v)
    return ("hkl" if star else "uvw", v)


def _uvw_label(u):
    """[uvw] with small integers, or None"""
    u = u / np.abs(u).max()
    for m in range(1, 13):
        w = u * m
        if np.allclose(w, np.round(w), atol=1e-6):
            return "[%d %d %d]" % tuple(np.round(w).astype(int))
    return None


def key_image(ubi, group, view="default", key="tsl", projection="stereographic", n=301):
    """
    The colour key of ubi_to_rgb as images of the sector, looking down the
    out-of-screen axis of the view with its east axis to the right.

    ubi gives the lattice (only its symmetrised metric is used).
    projection : "stereographic" (equal angle) or "equal area"
    Returns a list of (n,n,4) RGBA images, alpha 0 outside the sector: the upper
    hemisphere, and the lower one (seen from above) when the sector reaches it;
    and a list of (panel, x, y, label) for the corners.
    """
    name = group_name(group) if hasattr(group, "group") else group
    k = _Key(name, key, view)
    gs = symmetric_metric(ubi, point_group(name))[0]
    L = np.linalg.cholesky(gs)
    Li = np.linalg.inv(L)
    oos, east = ("+c*", "+a") if view == "default" else view.split("|")
    T = GROUPS[name].get("hexagonal_axes")
    z = _vec(L, Li, _view_axis(oos, T))[0]
    e = _vec(L, Li, _view_axis(east, T))[0]
    e = _unit(e - np.dot(e, z) * z)
    y = np.cross(z, e)
    X, Y = np.meshgrid(np.linspace(-1, 1, n), np.linspace(-1, 1, n))
    R = np.hypot(X, Y)
    disk = R <= 1
    if projection == "stereographic":
        theta = 2 * np.arctan(np.minimum(R, 1))
    else:
        theta = 2 * np.arcsin(np.minimum(R, 1) / np.sqrt(2))
    phi = np.arctan2(Y, X)
    # directions in the frame of L: x = L^T.uvw, so uvw = x.inv(L)
    N = np.array([_vec(L, Li, t) for t in k.normals]).reshape(-1, 1, 3).transpose(1, 0, 2)
    images = []
    for side in (1, -1):
        d = (np.sin(theta) * np.cos(phi))[..., None] * e + (np.sin(theta) * np.sin(phi))[..., None] * y + \
            (side * np.cos(theta))[..., None] * z
        uvw = np.dot(d.reshape(-1, 3), Li[0])
        inside = disk.ravel() & (np.dot(d.reshape(-1, 3), N[0].T) >= -1e-9).all(axis=1)
        if side == -1 and not (inside & (R.ravel() < 0.999)).any():
            break
        img = np.zeros((n * n, 4))
        img[inside, :3] = _uvw_to_rgb(k, uvw[inside], np.repeat(gs, inside.sum(), axis=0))
        img[inside, 3] = 1
        images.append(img.reshape(n, n, 4))
    labels = []
    corners, _ = _corners(N)
    for i, j, sign in corners:
        v = _unit(sign * np.cross(N[0, i], N[0, j]))
        panel = 0 if np.dot(v, z) >= -1e-9 else 1
        if panel >= len(images):
            continue
        ct = min(abs(np.dot(v, z)), 1.0)
        r = np.tan(np.arccos(ct) / 2) if projection == "stereographic" else np.sqrt(2) * np.sin(np.arccos(ct) / 2)
        p = np.arctan2(np.dot(v, y), np.dot(v, e))
        labels.append((panel, r * np.cos(p), r * np.sin(p), _uvw_label(np.dot(v, Li[0]))))
    return images, labels


def plot_key(ubi, group, view="default", key="tsl", projection="stereographic", n=301, axes=None,
             alpha=1.0, zoom=True):
    """Draw key_image with matplotlib (alpha fades it, zoom shows the sector
    rather than the whole projection), on axes with one per panel.
    Returns the axes used."""
    import matplotlib.pyplot as plt
    images, labels = key_image(ubi, group, view, key, projection, n)
    if axes is None:
        axes = plt.subplots(1, len(images), squeeze=False)[1][0]
    if len(axes) < len(images):
        raise ValueError("This key has %d panels (upper and lower hemisphere): give %d axes" % (
            len(images), len(images)))
    axes = axes[:len(images)]
    name = group_name(group) if hasattr(group, "group") else group
    for i, (ax, img) in enumerate(zip(axes, images)):
        rows, cols = np.nonzero(img[..., 3] > 0)
        img[..., 3] *= alpha
        ax.imshow(img, extent=(-1, 1, -1, 1), origin="lower")
        ax.add_patch(plt.Circle((0, 0), 1, fill=False, color="k", lw=0.5))
        for panel, x, y, s in labels:
            if panel == i and s is not None:
                ax.annotate(s, (x, y), fontsize=8, ha="center", va="bottom")
        g = np.linspace(-1, 1, n)
        x, y = (g[cols], g[rows]) if zoom else (g, g)
        pad = 0.1 * max(np.ptp(x), np.ptp(y))
        ax.set_xlim(x.min() - pad, x.max() + pad)
        ax.set_ylim(y.min() - pad, y.max() + pad)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title("%s %s %s%s" % (name, key.upper(), view, " (lower)" if i else ""), fontsize=9)
    return axes


def ipf_points(ubis, group, view="default", key="tsl", axis=(0.0, 0.0, 1.0), projection="stereographic",
               chunk=4096):
    """
    Where each ubi goes in the drawing of key_image: the crystal direction
    along axis, folded into the sector as for ubi_to_rgb, projected about the
    out-of-screen axis of the view with its east axis to the right.
    Returns x, y (n,) and panel (n,): 1 for the lower hemisphere (seen from above).
    """
    name = group_name(group) if hasattr(group, "group") else group
    k = _Key(name, key, view)
    ubis = np.asarray(ubis, float).reshape(-1, 3, 3)
    gs = symmetric_metric(ubis, point_group(name))[0]
    uvw = np.einsum("nji,j->ni", np.linalg.inv(ubis), np.asarray(axis, float))
    oos, east = ("+c*", "+a") if view == "default" else view.split("|")
    T = GROUPS[name].get("hexagonal_axes")
    out = np.empty((len(ubis), 3))
    for i in range(0, len(ubis), chunk):
        sl = slice(i, i + chunk)
        L = np.linalg.cholesky(gs[sl])
        Li = np.linalg.inv(L)
        x = _unit(np.einsum("nji,nj->ni", L, k.fold(uvw[sl], gs[sl])))
        z = _vec(L, Li, _view_axis(oos, T))
        e = _vec(L, Li, _view_axis(east, T))
        e = _unit(e - _dot(e, z)[:, None] * z)
        cz = _dot(x, z)
        theta = np.arctan2(np.sqrt(np.maximum(1 - cz * cz, 0)), np.abs(cz))
        r = np.tan(theta / 2) if projection == "stereographic" else np.sqrt(2) * np.sin(theta / 2)
        phi = np.arctan2(_dot(x, np.cross(z, e)), _dot(x, e))
        out[sl] = np.array([r * np.cos(phi), r * np.sin(phi), cz < 0]).T
    return out[:, 0], out[:, 1], out[:, 2].astype(int)


def plot_ipf(ubis, group, view="default", key="tsl", axis=(0.0, 0.0, 1.0), projection="stereographic",
             axes=None, n=301, alpha=0.3, **scatter_kwargs):
    """
    Inverse pole figure: each ubi as a point in its colour (ubi_to_rgb) on the
    key (plot_key, faded by alpha, drawn in the metric of the first ubi).
    Returns the axes.
    """
    ubis = np.asarray(ubis, float).reshape(-1, 3, 3)
    rgb = ubi_to_rgb(ubis, group, view, key, axis)
    x, y, panel = ipf_points(ubis, group, view, key, axis, projection)
    axes = plot_key(ubis[0], group, view, key, projection, n, axes, alpha)
    kw = {"s": 4, "edgecolors": "k", "linewidths": 0.2}
    kw.update(scatter_kwargs)
    for i, ax in enumerate(axes):
        # plot_key returns one axes per panel: without a lower panel the
        # sector only touches the equator
        m = (panel == i) | ((panel >= len(axes)) & (i == 0))
        ax.scatter(x[m], y[m], c=rgb[m], **kw)
        ax.set_title(ax.get_title() + ", axis %s" % (np.round(np.asarray(axis, float), 3),), fontsize=9)
    return axes
