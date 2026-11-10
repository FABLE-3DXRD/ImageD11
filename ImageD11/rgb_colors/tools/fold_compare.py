"""How often does MTEX's fold (image nearest the centre) leave the sector, for the keys in keys.json?"""
from __future__ import print_function, division
import json, numpy as np, ipf_ref as R
K = json.load(open("keys.json"))["groups"]
rng = np.random.RandomState(0)
for name, grp in K.items():
    M = R.group_uvw(grp["generators"]); g = R.cell_metric(grp["test_cell"])
    res = []
    for kind in ("tsl", "hsv"):
        MM = np.concatenate([M, -M]) if kind == "tsl" else M
        F = R.Frame(R.symmetrise(g, MM)); Rc = F.ops(MM)
        key = grp[kind][0]                 # the default view
        N = np.array([F.vec(t) for t in key["normals"]]).reshape(-1, 3)
        if not len(N):
            res.append("%s -" % kind); continue
        c = R.center(N, R.vertices(N))
        x = R.unit(rng.normal(size=(20000, 3)))
        out = np.dot(N, R.fold_nearest(x, Rc, c).T).min(0) < -1e-9
        res.append("%s outside %.4f" % (kind, out.mean()))
    print("%-6s %s" % (name, "  ".join(res)))
