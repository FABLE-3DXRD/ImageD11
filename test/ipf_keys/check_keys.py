"""Which MTEX keys are fundamental regions?  Metric free: v (uvw) is in a sector when N.v >= 0
for every normal N (hkl); a fundamental region holds exactly one image M v of a generic v."""
from __future__ import print_function, division
import json, sys, numpy as np
rng = np.random.RandomState(0)

def unique_ops(M):
    out = []
    for m in M:
        if not any(np.array_equal(m, o) for o in out):
            out.append(m)
    return np.array(out)

def count_inside(N, ops, A, n=4000):
    x = rng.normal(size=(n, 3))
    v = np.dot(x, np.linalg.inv(A))                     # random directions as uvw (A rows a,b,c)
    imgs = np.einsum("kij,nj->kni", ops, v)       # uvw' = M uvw
    inside = np.all(np.einsum("mi,kni->kmn", N, imgs) >= -1e-9, axis=1)  # k images, n dirs
    return inside.sum(0)

def main(E):
  G = {e["group"]: e for e in E if "ops" in e}
  res = {}
  for e in E:
      if "ops" in e:
          continue
      g = G[e["group"]]
      ops = np.array(g["ops"]).reshape(-1, 3, 3)
      A = np.array(g["A"]).reshape(3, 3)
      laue = unique_ops(np.concatenate([ops, -ops]))
      for key, M in (("tsl", laue), ("hsv", ops)):
          if isinstance(e[key], str):
              res.setdefault((e["group"], key), []).append((e["how2plot"], "error"))
              continue
          N = np.array(e[key]["N"]["hkl"], float).reshape(-1, 3)
          c = count_inside(N, M, A)
          ok = "ok" if np.all(c == 1) else "bad(%.2f..%d)" % (c.mean(), c.max())
          res.setdefault((e["group"], key), []).append((e["how2plot"], ok))
  for (g, key), r in res.items():
      d = dict(r)
      bad = [v for v, s in r if s != "ok"]
      print("%-6s %s default:%-14s  ok %2d/%d views   not ok: %s" % (g, key, d["default"], len(r) - len(bad), len(r),
            " ".join(bad) if len(bad) < 10 else "%d views" % len(bad)))


if __name__ == "__main__":
    main(json.load(open(sys.argv[1])))
