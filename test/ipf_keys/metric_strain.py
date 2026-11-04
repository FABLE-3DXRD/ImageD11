"""Symmetrised metric and the strain that separates a measured metric from it."""
from __future__ import print_function, division
import numpy as np
from ImageD11 import sym_u, unitcell

def symmetric_metric(ubi, grp):
    """g_s = <O g O^T> over the group (sym_u ops act as O.ubi), g = ubi.ubi^T"""
    g = np.dot(ubi, ubi.T)
    return np.mean([np.dot(np.dot(o, g), o.T) for o in grp.group], axis=0)

def asymmetric_strain(ubi, grp):
    """principal Green-Lagrange strains of the metric relative to the symmetrised one:
    the stretch of any direction uvw is uvw.g.uvw / uvw.g_s.uvw, extremes are eig(g_s^-1 g)"""
    g = np.dot(ubi, ubi.T)
    gs = symmetric_metric(ubi, grp)
    return 0.5 * (np.sort(np.linalg.eigvals(np.linalg.solve(gs, g)).real) - 1)

rng = np.random.RandomState(1)
for name, grp, cell in (("hexagonal", sym_u.hexagonal(), [3.0, 3.0, 9.0, 90, 90, 120]),
                        ("cubic", sym_u.cubic(), [4.0, 4.0, 4.0, 90, 90, 90])):
    B = unitcell.unitcell(cell).B
    U = np.linalg.qr(rng.normal(size=(3, 3)))[0]
    ubi0 = np.linalg.inv(np.dot(U, B))
    gs0 = symmetric_metric(ubi0, grp)
    print(name, "ideal: max|e| %.1e, O g_s O^T = g_s to %.1e" % (np.abs(asymmetric_strain(ubi0, grp)).max(),
          max(np.abs(np.dot(np.dot(o, gs0), o.T) - gs0).max() for o in grp.group)))
    for eps in (1e-3, 5e-3):
        S = rng.normal(size=(3, 3)); S = 0.5 * (S + S.T); S -= np.trace(S) / 3 * np.eye(3)
        S *= eps / np.abs(np.linalg.eigvalsh(S)).max()          # deviatoric, max principal = eps
        F = np.eye(3) + S                                       # small strain applied in the lab
        ubi = np.dot(ubi0, F.T)
        e = asymmetric_strain(ubi, grp)
        print("   applied deviatoric strain %.0e -> measured principal %s" % (eps, np.round(e, 5)))
    # a strain the lattice allows: c/a and volume change keep the symmetry
    ubi = np.dot(np.diag([1.002, 1.002, 0.995]), ubi0) if name == "hexagonal" else 1.01 * ubi0
    print("   allowed change (c/a or volume): max|e| %.1e" % np.abs(asymmetric_strain(ubi, grp)).max())
