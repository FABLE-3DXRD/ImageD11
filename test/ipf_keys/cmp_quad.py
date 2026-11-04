"""rgb_colors (Gauss on arcs) vs ipf_ref trapezoid with 1000 (MTEX) and 200000 steps, all keys, random directions"""
from __future__ import print_function, division
import json, sys, numpy as np
sys.path.insert(0, "../../test")
import ipf_ref
from test_rgb_colors import cell_ubi, turn_to_z
from ImageD11 import rgb_colors
K = json.load(open("keys.json"))["groups"]
O = json.load(open("../../test/ipf_oracle.json"))
rng = np.random.RandomState(5)
ca = ipf_ref.correct_angle
w_fine, w_1000 = 0, 0
for name, g in K.items():
    u0 = cell_ubi(g["test_cell"])
    uvw = rng.normal(size=(300, 3)).dot(np.linalg.inv(u0).T)   # random lab directions as uvw
    ubis = [np.dot(u0, turn_to_z(np.dot(u0.T, u), 0.0).T) for u in uvw]
    for kind in ("tsl", "hsv"):
        for key in g[kind]:
            mine = rgb_colors.ubi_to_rgb(ubis, name, key["views"][0], kind)
            ipf_ref.correct_angle = lambda *a, **k: ca(*a, n=200000)
            fine = ipf_ref.colours(g, kind, key, g["test_cell"], uvw)
            ipf_ref.correct_angle = ca
            m1000 = ipf_ref.colours(g, kind, key, g["test_cell"], uvw)
            e1, e2 = np.abs(mine - fine).max(), np.abs(mine - m1000).max()
            if e1 > w_fine: w_fine = e1; at1 = (name, kind, key["views"][0])
            if e2 > w_1000: w_1000 = e2; at2 = (name, kind, key["views"][0])
print("vs 200000 steps: %.2e %s" % (w_fine, at1))
print("vs 1000 steps (MTEX): %.2e %s" % (w_1000, at2))
