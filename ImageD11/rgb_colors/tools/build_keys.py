"""keys.json: the MTEX 7.2.0 direction colour keys in lattice terms, from mtex72_keys.json.
Integer vectors tagged with their basis, ITA generators, a few check colours per key."""
from __future__ import print_function, division
import json, sys, warnings, numpy as np
import ipf_ref
warnings.simplefilter("ignore")
from ImageD11 import sym_u
from make_table import rat, mtex_keys

def tag(S, prefer):
    """{"hkl": ints} or {"uvw": ints}, preferring the basis given"""
    for b in (prefer, "uvw" if prefer == "hkl" else "hkl"):
        r = rat(np.asarray(S[b], float).ravel())
        if r:
            return {b: list(r)}
    raise ValueError("no integer indices: %s" % S)

def rows(S, prefer):
    U = np.asarray(S["uvw"], float).reshape(-1, 3); H = np.asarray(S["hkl"], float).reshape(-1, 3)
    return [tag({"uvw": u, "hkl": h}, prefer) for u, h in zip(U, H)]

def generators(M):
    """a few ITA triplets (acting on uvw) that generate the group of uvw matrices M"""
    want = {tuple(np.round(m).astype(int).ravel()) for m in M}
    if len(want) == 1:
        return ["x,y,z"]
    gens = []
    for m in M:
        if np.allclose(m, np.eye(3)):
            continue
        have = {tuple(np.round(o.T).astype(int).ravel()) for o in sym_u.generate_group(*(gens or ["x,y,z"])).group}
        if tuple(np.round(m).astype(int).ravel()) not in have:
            gens.append(sym_u.m_to_string(m.T))
        have = {tuple(np.round(o.T).astype(int).ravel()) for o in sym_u.generate_group(*(gens or ["x,y,z"])).group}
        if have == want:
            break
    assert have == want, (gens, len(have), len(want))
    return gens or ["x,y,z"]

E = json.load(open(sys.argv[1]))
G = {e["group"]: e for e in E if "ops" in e}
K = mtex_keys(E)
NCHECK = 8
out = {"source": "MTEX 7.2.0 (github mtex-toolbox, commit 8c7ebb7), fundamentalSector, TSLDirectionKey, "
                 "HSVDirectionKey with default options, every view out-of-screen|east whose sector is a fundamental region",
       "excluded": "valid MTEX views whose key is degenerate are listed per group under excluded",
       "conventions": "uvw: direction = u a + v b + w c; hkl: plane normal = h a* + k b* + l c*. "
                      "generators: ITA triplets acting on uvw (sym_u.generate_group). "
                      "check: test directions (uvw, test cell) and MTEX colours.",
       "groups": {}}
for name, (meta, res) in K.items():
    g = G[name]
    M = np.array(g["ops"]).reshape(-1, 3, 3)
    uvw_all = np.array(g["dirs_uvw"], float).reshape(-1, 3)
    uvw = uvw_all[:NCHECK]
    es = {e["how2plot"]: e for e in E if e.get("group") == name and "how2plot" in e}
    entry = {"laue": meta["laue"], "generators": generators(M), "order": len(M),
             "test_cell": meta["abc"] + meta["angles"], "check_uvw": np.round(uvw, 9).tolist(), "tsl": [], "hsv": [],
             "excluded": []}
    for key in ("tsl", "hsv"):
        for rgb, desc, views in res[key]:
            e = es[views[0]]; p = e[key]
            k = {"views": views,
                 "out_of_screen": tag(e["outOfScreen"], "hkl"),
                 "normals": rows(p["N"], "hkl"),
                 "vertices": rows(p["vertices"], "uvw"),
                 "hue_ref": tag(p["ref"], "uvw")}
            if key == "hsv":
                k["split"] = rows(p["refl"], "hkl")
            k["check_rgb"] = np.round(rgb[:NCHECK], 6).tolist()
            try:
                d = np.nanmax(np.abs(ipf_ref.colours(entry, key, k, entry["test_cell"], uvw_all) - rgb), axis=1)
            except ValueError as err:
                entry["excluded"].append({"key": key, "views": views, "reason": str(err)})
                continue
            if (d > 0.01).any():
                entry["excluded"].append({"key": key, "views": views, "reason":
                    "MTEX colours %d of %d test directions outside its own sector" % ((d > 0.01).sum(), len(d))})
                continue
            entry[key].append(k)
    if not entry["excluded"]:
        del entry["excluded"]
    out["groups"][name] = entry
s = json.dumps(out, separators=(",", ":"))
open("keys.json", "w").write(s)
print(len(s), "bytes;", sum(len(v[k]) for v in out["groups"].values() for k in ("tsl", "hsv")), "keys")
