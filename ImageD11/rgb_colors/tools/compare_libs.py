"""pymtex (TSL, HSV) and orix (TSL) colours on the directions of mtex72_keys.json, with their sectors in
lattice terms.  Each library gets the direction as uvw turned into its own Cartesian frame."""
from __future__ import print_function, division
import json, sys, warnings, numpy as np
warnings.simplefilter("ignore")
from mtex import crystalFrame
from mtex.plotting.directionkeys import HSVDirectionKey, TSLDirectionKey
import orix.quaternion.symmetry as OS
from orix.plot import DirectionColorKeyTSL
from orix.vector import Vector3d
from diffpy.structure import Lattice

ORIX = {"1": "C1", "-1": "Ci", "211": "C2x", "121": "C2y", "112": "C2z", "m11": "Csx", "1m1": "Csy",
        "11m": "Csz", "112/m": "C2h", "222": "D2", "mm2": "C2v", "mmm": "D2h", "4": "C4", "-4": "S4",
        "4/m": "C4h", "422": "D4", "4mm": "C4v", "-42m": "D2d", "4/mmm": "D4h", "3": "C3", "-3": "S6",
        "321": "D3x", "312": "D3y", "3m1": "C3v", "31m": "C3v", "-3m1": "D3d", "-31m": "D3d", "6": "C6",
        "-6": "C3h", "6/m": "C6h", "622": "D6", "6mm": "C6v", "-6m2": "D3h", "-62m": "D3h", "6/mmm": "D6h",
        "23": "T", "m-3": "Th", "432": "O", "-43m": "Td", "m-3m": "Oh"}   # no orix 2/m11, 12/m1, 2mm, m2m


def lat(v, A):
    """uvw and hkl of Cartesian rows v for lattice rows A"""
    v = np.asarray(v, float).reshape(-1, 3)
    return {"uvw": np.dot(v, np.linalg.inv(A)).tolist(), "hkl": np.dot(v, A.T).tolist()}


def desc(N, V, c, A):
    from make_table import vec
    return "N %s | vertices %s | center %s" % (vec(lat(N, A)["hkl"], "hkl"), vec(lat(V, A)["uvw"], "uvw"),
                                                 vec(lat(c, A)["uvw"], "uvw"))


def data(x):
    x = x() if callable(x) else x
    return np.asarray(getattr(x, "data", x), float).reshape(-1, 3)


E = json.load(open(sys.argv[1]))
out = {}
for g in [e for e in E if "ops" in e]:
    name = g["group"]
    e0 = next(e for e in E if e.get("group") == name and "how2plot" in e)
    abc, ang = np.array(e0["abc"], float), np.array(e0["angles"], float)
    uvw = np.array(g["dirs_uvw"], float).reshape(-1, 3)
    res = out[name] = {}
    # pymtex, default plotting convention
    cs = crystalFrame(name, tuple(abc), np.radians(ang))
    A = np.asarray(cs.axes.data, float).reshape(3, 3)
    x = np.dot(uvw, A); x /= np.linalg.norm(x, axis=1)[:, None]
    r = {}
    for key, K, sym in (("tsl", TSLDirectionKey, cs.Laue), ("hsv", HSVDirectionKey, cs)):
        try:
            sym = sym() if callable(sym) else sym
            sR = sym.fundamentalSector()
            r[key] = {"rgb": np.asarray(K(cs).direction2color(x), float).tolist(),
                      "desc": desc(data(sR.N), data(sR.vertices), data(sR.center), A)}
        except Exception as err:
            r[key] = "error: %r" % err
    res["pymtex"] = r
    # orix, TSL only, frame of diffpy's Lattice (x || a*, z || c)
    if name in ORIX:
        sym = getattr(OS, ORIX[name])
        A = Lattice(*(list(abc) + list(ang))).base
        x = np.dot(uvw, A)
        try:
            k = DirectionColorKeyTSL(sym)
            fs = k.symmetry.fundamental_sector if hasattr(k, "symmetry") else sym.laue.fundamental_sector
            res["orix " + ORIX[name]] = {"tsl": {
                "rgb": np.asarray(k.direction2color(Vector3d(x)), float).tolist(),
                "desc": desc(fs.data, fs.vertices.data, fs.center.data, A)}}
        except Exception as err:
            res["orix " + ORIX[name]] = {"tsl": "error: %r" % err}
    else:
        res["orix"] = {"tsl": "no orix point group in this setting"}
    print(name, end=" ", flush=True)
json.dump(out, open("libs_rgb.json", "w"))
