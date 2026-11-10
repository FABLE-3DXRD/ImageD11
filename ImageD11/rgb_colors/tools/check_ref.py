"""ipf_ref.py (keys.json + formulas only) against MATLAB MTEX 7.2 colours on all 500 test directions."""
from __future__ import print_function, division
import json, sys, numpy as np
import ipf_ref
ipf_ref.MTEX_F = True
from make_table import mtex_keys
E = json.load(open("mtex72_keys.json")); K = json.load(open("keys.json"))["groups"]
D = {e["group"]: np.array(e["dirs_uvw"], float).reshape(-1, 3) for e in E if "ops" in e}
MK = mtex_keys(E)
only = sys.argv[1:]
for name, grp in K.items():
    if only and name not in only:
        continue
    line = []
    for kind in ("tsl", "hsv"):
        worst = []
        bv = {v: k[0] for k in MK[name][1][kind] for v in k[2]}
        for key in grp[kind]:
            rgb = bv[key["views"][0]]
            try:
                mine = ipf_ref.colours(grp, kind, key, grp["test_cell"], D[name])
                d = np.nanmax(np.abs(mine - rgb), axis=1)
                worst.append((np.nanmax(d), (d > 0.01).sum()))
            except Exception as err:
                worst.append((np.inf, repr(err)[:40]))
        m = max(worst, key=lambda w: w[0])
        line.append("%s %d keys: max %.1e (n>0.01: %s)" % (kind.upper(), len(worst), m[0], m[1]))
    print("%-6s %s" % (name, " | ".join(line)))
