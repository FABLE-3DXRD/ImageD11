"""MTEX 7.2 colour keys in lattice terms: one entry per distinct key over the valid views
(those whose sector is a fundamental region, see check_keys.py)."""
from __future__ import print_function, division
import json, sys, numpy as np
from fractions import Fraction
import check_keys as C

def vec(a, kind):
    a = np.asarray(a, float).reshape(-1, 3)
    out = []
    for v in a:
        r = rat(v)
        o, c = ("[", "]") if kind == "uvw" else ("(", ")")
        out.append(o + " ".join("%d" % x for x in r) + c if r else
                   "~" + o + " ".join("%.4f" % x for x in v / np.abs(v).max()) + c)
    return " ".join(out) if out else "-"

def rat(v, tol=1e-6, dmax=12):
    """smallest integer multiple of v, or None if v is not rational with small indices"""
    v = np.asarray(v, float)
    m = np.abs(v).max()
    if m < 1e-12:
        return None
    v = v / m
    den = 1
    for x in v:
        f = Fraction(x).limit_denominator(dmax)
        if abs(float(f) - x) > tol:
            return None
        den = den * f.denominator // np.gcd(den, f.denominator)
    w = np.round(v * den).astype(int)
    g = np.gcd.reduce(np.abs(w[w != 0]))
    return tuple(int(x) for x in w // g)

ROLE = {"N": "hkl", "refl": "hkl", "vertices": "uvw", "center": "uvw", "whiteCenter": "uvw", "ref": "uvw"}
FIELDS = {"tsl": ["N", "vertices", "center", "ref"], "hsv": ["N", "refl", "vertices", "whiteCenter", "ref"]}

def hue0(e, part):
    """the hue reference MTEX uses (polarCoordinates.m): out-of-screen unless it is parallel to the centre"""
    c = np.asarray(part["whiteCenter" if "whiteCenter" in part else "center"]["uvw"], float).ravel()
    o = np.asarray(e["outOfScreen"]["uvw"], float).ravel()
    if np.linalg.norm(np.cross(c, o)) > 1e-9 * np.linalg.norm(c) * np.linalg.norm(o):
        return "out-of-screen " + vec(e["outOfScreen"]["hkl"], "hkl")
    return "ref"


def mtex_keys(E):
    """{group: (meta, {"tsl": [[rgb, description, views]], "hsv": [...]})}; views with equal colours are one key"""
    G = {e["group"]: e for e in E if "ops" in e}
    out = {}
    for name, g in G.items():
        es = [e for e in E if e.get("group") == name and "how2plot" in e]
        ops = np.array(g["ops"]).reshape(-1, 3, 3); A = np.array(g["A"]).reshape(3, 3)
        laue = C.unique_ops(np.concatenate([ops, -ops]))
        res = {}
        for key, M in (("tsl", laue), ("hsv", ops)):
            keys = res[key] = []
            for e in es:
                if isinstance(e[key], str) or isinstance(e[key]["rgb"], str):
                    continue
                N = np.array(e[key]["N"]["hkl"], float).reshape(-1, 3)
                if not np.all(C.count_inside(N, M, A, 2000) == 1):
                    continue
                rgb = np.array(e[key]["rgb"], float).reshape(-1, 3)
                for k in keys:
                    if np.nanmax(np.abs(k[0] - rgb)) < 1e-6:
                        k[2].append(e["how2plot"])
                        break
                else:
                    d = " | ".join("%s %s" % (f, vec(e[key][f][ROLE[f]], ROLE[f])) for f in FIELDS[key])
                    keys.append([rgb, d + " | hue0 " + hue0(e, e[key]), [e["how2plot"]]])
        out[name] = (dict(laue=g["laue"], nops=len(ops), abc=es[0]["abc"], angles=es[0]["angles"]), res)
    return out


def match(rgb, keys, tol=1e-3):
    """index of the MTEX key closest to rgb, and the largest colour difference to it"""
    d = [np.nanmax(np.abs(k[0] - rgb)) for k in keys]
    i = int(np.argmin(d))
    return i, d[i]


if __name__ == "__main__":
    E = json.load(open(sys.argv[1]))
    libs = json.load(open(sys.argv[2])) if len(sys.argv) > 2 else {}
    print("# MTEX 7.2.0 direction colour keys, lattice terms\n")
    print("Normals (hkl), directions [uvw]; ~ marks a direction with irrational indices for this cell.")
    print("A view is out-of-screen|east. Only views whose sector is a fundamental region are listed.")
    print("hue0 is the direction hue 0 points to from the centre. Colours compared on 500 random directions.\n")
    for name, (meta, res) in mtex_keys(E).items():
        print("## %s\n\nLaue %s, %d operations, test cell %s %s\n" % (name, meta["laue"], meta["nops"], meta["abc"], meta["angles"]))
        for key in ("tsl", "hsv"):
            keys = res[key]
            notes = [[] for _ in keys]
            extra = []
            for lib, d in libs.get(name, {}).items():
                if key not in d:
                    continue
                if isinstance(d[key], str):
                    extra.append("  - %s: %s" % (lib, d[key]))
                    continue
                rgb = np.array(d[key]["rgb"], float).reshape(-1, 3)
                i, dm = match(rgb, keys)
                if dm < 1e-3:
                    notes[i].append("%s (max |dRGB| %.1e)" % (lib, dm))
                else:
                    extra.append("  - **%s %s** %s  \n    nearest MTEX key %d, max |dRGB| %.3f on %d/%d directions > 0.01" % (
                        lib, key.upper(), d[key].get("desc", ""), i + 1, dm,
                        int((np.nanmax(np.abs(keys[i][0] - rgb), 1) > 0.01).sum()), len(rgb)))
            for j, k in enumerate(keys):
                print("%d. **%s** %s  \n   views: %s%s" % (j + 1, key.upper(), k[1], " ".join(k[2]),
                      ("  \n   same colours: " + ", ".join(notes[j])) if notes[j] else ""))
            if extra:
                print("\n   Not an MTEX key:\n" + "\n".join(extra))
        print()
