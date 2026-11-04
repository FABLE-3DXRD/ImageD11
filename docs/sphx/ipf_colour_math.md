# IPF colour keys in crystal space

The TSL and HSV inverse pole figure keys of MTEX 7.2.0, written with lattice
vectors and a metric only. `ImageD11.rgb_colors` implements this file:
`ubi_to_rgb(ubis, group, view, key, axis)` gives the colours and `plot_key` draws
the legend. The keys are in `ImageD11/ipf_keys.py`. The scripts named below are
in `test/ipf_keys/`: `keys.json` holds the keys as dumped from MATLAB, and
`ipf_ref.py`, an independent reference, follows this file line by line.

Check (`check_ref.py`): MATLAB MTEX 7.2.0 colours for 500 random directions per
key, for all 305 keys (45 point group settings, every view that gives a key).
`ipf_ref.py` with `keys.json` (and `MTEX_F = True`, 5.7) reproduces them to at
most 6.3e-11, with two exceptions: m-3 and 23 TSL (6e-3, the centre, section
5.3) and one direction for 432 HSV (4.5e-3, the fold, section 5.1).

`ImageD11.rgb_colors` is checked (`test/test_rgb_colors.py`) against
`test/ipf_oracle.json`: 8 directions per key, for every key and view, each put
along the sample axis by turning a grain. It matches `ipf_ref.py` with 200000
integration steps to 1e-7. TSL matches MATLAB to 1.6e-5, except m-3 and 23:
4.7e-3 (the centre, 5.3). HSV differs from MATLAB by up to 9.5e-3, from the
hue correction (5.7), and by 5.5e-5 otherwise (the hue integral, 5.5).

Source: github.com/mtex-toolbox/mtex at 8c7ebb7. The files are
`fundamentalSector.m`, `TSLDirectionKey.m`, `HSVDirectionKey.m`,
`polarCoordinates.m`, `sphericalRegion.m` (center, vertices) and `ar2rgb.m`.
The method is from G. Nolze and R. Hielscher, J. Appl. Cryst. 49 (2016) 1786,
doi:10.1107/S1600576716012942.

## 1. Vectors

| symbol | meaning |
|---|---|
| `ubi` | rows are the lattice vectors a, b, c (lab frame) |
| `UB = ubi^-1` | columns are a*, b*, c* |
| r | unit sample direction in the lab |
| uvw | r = u a + v b + w c, so uvw = UB^T r (rows: `r . ubi^-1`) |
| hkl | h = r.a, k = r.b, l = r.c, so hkl = ubi r. A plane normal (hkl) is h a* + k b* + l c* |
| g | g = ubi ubi^T, g_ij = a_i . a_j |

Dot products need no frame:

- direction with direction: u^T g u'
- normal with direction: h . u (no metric)
- normal with normal: h^T g^-1 h'

The angle between a reference direction and r is
`cos = ref_uvw . (ubi r) / (sqrt(ref_uvw^T g ref_uvw) |r|)`.

Cross products: a_i x a_j = V eps_ijk a*_k. So the cross product of two
directions [uvw] x [u'v'w'], taken componentwise, gives a plane normal (hkl)
times V, and (hkl) x (h'k'l') gives a direction [uvw] divided by V. This is why
every sector vertex, the intersection of two lattice planes, has integer uvw.

## 2. Symmetry operations

`generators` in `keys.json` are ITA coordinate triplets. They act on uvw as
written: "-y,x-y,z" sends [u v w] to [-v, u-v, w]. Writing M for that matrix:

- directions: uvw' = M uvw
- plane normals: hkl' = M^-T hkl
- lattice: ubi' = M^-T ubi
- metric: g' = M^T g M

`sym_u.generate_group(*generators)` builds the group. Its matrices are
O = M^T, which act on ubi as `O.ubi` and on hkl as `O.hkl`, so a direction
transforms with O^-T.

Folding uses only these integer matrices. They map lattice vectors to lattice
vectors exactly, whatever the metric.

The **TSL key folds with the Laue group**, {M} together with {-M}. The **HSV key
folds with the point group** as given, improper operations (det -1) included.

## 3. Symmetrised metric and the strain warning

The angles in sections 5.3 to 5.6 use the metric averaged over the group:

    g_s = (1/n) sum_M M^T g M        (equivalently <O g O^T> for O acting on ubi)

g_s is invariant: every operation is an isometry of it. It is the same for every
equivalent indexing of the grain, because g of M^-T ubi averages to the same
g_s. Inversion leaves g unchanged, so a point group and its Laue group give the
same g_s.

The deviation of g from g_s is a strain. The stretch of a direction between the
two metrics is u^T g u / u^T g_s u. Its extremes are the eigenvalues of
g_s^-1 g, so the principal strains are

    e = (eig(g_s^-1 g) - 1) / 2,      with trace(g_s^-1 g) = 3 exactly

These are the strains the symmetry forbids. A volume change, or a change of c/a
in a tetragonal or hexagonal cell, stays in g_s and gives e = 0.

**Warn when max|e| > 5e-3.** The metric is then far from the group, so either
the grain is strongly strained or the group is too high for the crystal.

Measured (`metric_strain.py`, random orientation, ideal cell plus a random
deviatoric strain):

| case | max\|e\| |
|---|---|
| ideal hexagonal or cubic | 3e-16 |
| hexagonal c/a change, cubic volume change | <1e-15 |
| cubic, 1e-3 deviatoric strain | 1.0e-3 |
| cubic, 5e-3 deviatoric strain | 4.99e-3 |
| hexagonal, 5e-3 deviatoric strain | 4.4e-3 (the c/a part of the strain is allowed) |

Why g_s, not the grain's own metric (`strain_fold.py`). The test takes a cubic
grain with random strain and indexes it all 24 ways, M^-T ubi, then compares
over 20000 sample directions:

| strain | fold in a Cartesian frame from the B convention | integer fold in uvw | angle to [111], own metric | angle to [111], g_s |
|---|---|---|---|---|
| 1e-4 | 0.008 deg | 0 | 0.025 deg | 5e-13 deg |
| 1e-3 | 0.26 deg | 0 | 0.14 deg | 3e-12 deg |
| 1e-2 | 1.65 deg | 0 | 1.7 deg | 9e-13 deg |

The colours computed this way (`ipf_ref.colours`, strain 3e-3) agree across all
equivalent indexings to at most 7e-11 for m-3m, 6/mmm and 4/mmm, in both keys.

## 4. What a key holds

Each group in `keys.json` has `generators`, `order`, `laue`, `test_cell` and
`check_uvw`, and lists of `tsl` and `hsv` keys. Each key has:

| field | meaning |
|---|---|
| `views` | MTEX plotting conventions (out-of-screen\|east) that give this key; "default" is MTEX's own |
| `out_of_screen` | the out-of-screen axis, used for the hue reference (5.4) |
| `normals` | sector boundaries n: the sector is n . x >= 0 for every n |
| `vertices` | the corners of the sector (they can be recomputed from `normals`, 5.2) |
| `hue_ref` | used when out-of-screen gives no hue reference. TSL: the first vertex in MTEX's order, or a* with no vertices. HSV: east |
| `split` (HSV) | extra planes that split the sector into a white half and a black half |
| `check_rgb` | MTEX colours of `check_uvw` in `test_cell` |

Every vector is tagged with its basis, `{"hkl": [...]}` or `{"uvw": [...]}`,
and has integer indices. A plane can also be given by the direction normal to it
(`uvw`). It then needs g_s to be turned into a normal, as in the 12/m1 HSV split
plane, which is perpendicular to [100].

The centre is not stored. It follows from the normals and vertices (5.3).

MTEX's order of the vertices, and so the TSL `hue_ref`, comes from `unique()` on
Cartesian coordinates. That is an artefact of MTEX's frame, so `hue_ref` is
stored rather than derived.

## 5. Colour of a direction

Inputs: a group, a key, ubi and r. Compute uvw from section 1 and g_s from
section 3. All lengths, angles and cross products below are taken in g_s.

`ipf_ref.py` does the arithmetic in orthonormal coordinates
x = L^T uvw, with g_s = L L^T (Cholesky). Plane normals become L^-1 hkl. This is
the only place a Cartesian frame appears. Any right handed factor of g_s gives
the same colours, because only dot and cross products enter and every operation
is an exact isometry of g_s. That is not true of the B matrix of the strained
cell.

### 5.1 Fold

Replace x by the image M x that lies inside the sector, n . (M x) >= 0 for every
normal n. The test is metric free (hkl . uvw), except for a plane given by the
direction normal to it, whose normal is g_s uvw. On a boundary any image gives
the same point.

MTEX folds differently: it takes the image nearest the centre (largest dot
product with the symmetric copies of the centre). That is not the sector for
-3, 3 (TSL) and -4, -3, 321 (HSV): 6 to 13 % of directions would land outside.
MTEX special-cases those groups (`project2FundamentalRegion.m`, ids 18, 19, 22,
26), so the sector is what it means.

For the 4-vertex sectors MTEX's centre is approximate (5.3), so its fold leaves
a thin sliver outside the sector. One of the 500 test directions for 432 HSV
falls there, 4.4e-4 outside, and gets a colour differing by 4.5e-3.

### 5.2 Vertices

The vertices are the vectors ±(n_i x n_j), over pairs of normals, that satisfy
every n . v >= 0. When both planes are given by hkl they have integer uvw
(section 1). `rgb_colors` finds them in g_s for every key, because of the
planes given by a direction (monoclinic HSV).

### 5.3 Centre c

The region here is the sector, or for HSV the sector with its `split` planes
added. The centre depends on its shape:

| region | centre | used by |
|---|---|---|
| one plane | that normal | TSL 1 and -1; HSV 1, -1, m11, 1m1, 11m |
| lune (2 vertices) | normalise(n_i/\|n_i\| + n_j/\|n_j\|) for the two normals furthest apart | monoclinic, 222, 3, 31m, 3m1, ... |
| triangle | normalise(sum of the unit vertices) | most keys |
| 4 vertices | the area centroid: normalise(sum over edges of arc length x inward unit normal), since the integral of x dA over a spherical polygon is half that sum | TSL m-3 and 23 |

For the 4-vertex case MTEX takes the mean of a 1 deg grid over the sector. That
is 0.094 deg from the area centroid and changes colours by up to 6e-3. With
MTEX's centre substituted, the formulas here reproduce MATLAB to 8e-12.

### 5.4 Polar coordinates about c

Radius:

    r = 1 - angle(c, x) / T

T is the distance from c, along the great circle through x, to where that circle
leaves the region (the region includes the split planes for HSV). For a plane n,
with tau the unit tangent at c towards x, the crossing is at
t = atan2(n.c, -n.tau), taken in (0, pi]. T is the smallest t over the planes.
So r = 1 at the centre and r = 0 on the boundary.

MTEX writes this as min over n of angle(-x, b)/angle(-c, b), with
b = unit((x x c) x n). The two forms agree to 3.5e-11.

Angles are atan2(|a x b|, a . b). arccos(a . b) loses half the digits near 0
and pi: on a boundary plane, where angle(-x, b) = pi, it gave colours differing
by 2e-8 between two equivalent descriptions.

Hue reference:

    rx = out_of_screen - c,  or hue_ref - c when out_of_screen is parallel to c

rx is then projected perpendicular to c, and ry = c x rx.

Azimuth:

    rho = atan2(ry . x, rx . x)  mod 2 pi

### 5.5 Even hue (MTEX correctAngle)

This applies whenever the region has vertices, which is every key except the
one-plane ones (those use rho as it is). Let w(rho) be the distance from c to
the boundary along azimuth rho. That is exactly T of 5.4 for the tangent
cos(rho) rx + sin(rho) ry, to 1.2e-15.

The hue angle is

    omega(rho) = 2 pi W(rho) / W(2 pi),   W(rho) = integral from 0 to rho of w

A triangle has an extra rule. Let rho_1, rho_2, rho_3 be the azimuths of its
three vertices. Each arc between consecutive vertices gets exactly a third of
the hue circle, and within an arc the hue is proportional to W.

MTEX integrates on 1000 equal steps of rho, with the vertex azimuths added, by
the trapezoid rule, then interpolates linearly. Against 200000 steps this
changes colours by at most 5.7e-5. `ipf_ref.py` uses the same 1000 steps.

`rgb_colors` computes the integral itself. With nn the normal of the great
circle through c at azimuth rho, nn = cos(rho) n0 + sin(rho) c x n0 and
n0 = unit(rx x c), the distance to a plane n is

    w_n(rho) = atan2(|n . c|, n . (nn x c)),   nn x c = cos(rho) n0 x c + sin(rho) n0

and w is the smallest w_n. w is smooth between the vertex azimuths, so each arc
between them (and between 0, the vertices, rho and 2 pi) is integrated by
32 point Gauss-Legendre. Measured over 300 random directions for every key:
16, 24, 32 and 64 points differ from 128 by 1.5e-7, 1.1e-10, 2.5e-13 and
2.6e-14; against `ipf_ref.py` with 200000 steps the colours differ by 1.4e-9.

The table omega(rho) is the same for every grain of a key. The sector is made
of lattice planes and directions, and the metrics a group allows only stretch
it along directions it is symmetric about: the dot products of the unit
vectors the colour uses (planes, vertices, c, the hue reference rx actually
used) are the same in every allowed metric (`TestHueTable`, 20 % strain, to
1e-12). A grain's metric still gives its own x, so its r and rho, exactly.
`rgb_colors` builds the table once per call from the first grain: each arc in
64 pieces, 32 point Gauss per piece; a grain's omega is the table at the knot
below rho plus 8 point Gauss on the rest of that piece. Against the
integral per grain the colours differ by 2.3e-13 (300 strained grains, every
key and view).

### 5.6 TSL colour

    hue = omega / 2 pi,   lightness L = 1/2 + r/2,   saturation 1,   no hue correction

So the centre is white and the boundary is fully saturated. Convert HSL to RGB
with 5.8.

### 5.7 HSV colour

Defaults: grayValue (0.2, 0.5), grayGradient gamma = 0.5, colorStretching
s = 1, no maxAngle.

1. Fold x with the point group into the sector, `normals` only (5.1).
2. For each split plane p: if p . x < 1e-5 (the far side), reflect,
   x = x - 2 (p.x) p for unit p, and mark x black. Directions that are never
   reflected are white. If the white centre itself is on the far side of p,
   reflect it and swap black and white. This does not happen in the table.
3. The white centre c is the centre (5.3) of the sector with the split planes
   added. The hue reference falls back to `hue_ref` = east.
4. Compute r and omega as in 5.4 and 5.5, over the region with split planes.
5. t = 1/2 + r/2 for white, t = (1 - r)/2 for black. Both are 1/2 on the
   boundary and on the split planes. White is t = 1 and black t = 0.
6. theta = pi t, then t' = ((2 gamma theta/pi + (1 - gamma)(1 - cos theta)) / 2)^s.
7. gray = 1 - 2 G |t' - 1/2|, with G = 0.2 when theta > pi/2 (white side) and
   G = 0.5 otherwise.
8. L = (t' - 1/2) gray + 1/2,
   S = gray (1 - |2t' - 1|) / (1 - |2L - 1|) (0 when 0/0),
   hue = Hc(omega / 2 pi). Convert with 5.8.

Hc is MTEX's hue correction, which widens red, green and blue:

    f(z) = 1/2 + sum over m in (0, 0.3333, 0.6666) of exp(-200 (((z - m + 1/2) mod 1) - 1/2)^2)
    on z = 0, 1/999, ..., 1
    F = integral of f from 0 (trapezoid), divided by its value at 1
    Hc(h) = linear interpolation of F at h

MTEX takes F = cumsum(f) / sum(f), so F(0) = f(0)/sum(f) = 1.7e-3, not 0,
while F(1) = 1, and its HSV hue jumps by 1.7e-3 where rho passes 2 pi, on the
half great circle from c along rx. `rgb_colors` starts F at 0, so the colour is
continuous there; it differs from MTEX by up to 9.5e-3 in rgb (measured on the
check directions).

### 5.8 HSL to RGB (MTEX hsl2hsv, then MATLAB hsv2rgb)

    l2 = 2L
    s2 = S l2        if l2 <= 1
       = S (2 - l2)  otherwise
    V = (l2 + s2)/2,   S_v = 2 s2 / (l2 + s2) (0 when 0/0)

Then the standard HSV to RGB with hue in [0, 1).

## 6. Views (MTEX how2plot)

A view is the out-of-screen axis and the east axis. They enter in three places:

- The sector hemisphere: out-of-screen.
- The sector half chosen by east for 2/m (2-fold in the screen), mm2 and mmm.
- The hue reference (5.4): out-of-screen, else `hue_ref`.

In the plot of the key, east only turns the picture.

Measured over out-of-screen in {±a*, ±b*, ±c*} and east in {±a, ±b, ±c}
perpendicular to it (`check_keys.py`):

- Only out-of-screen ±c* gives a fundamental region for every group. The
  exceptions are cubic (every view is valid), 1 and -1, and monoclinic, where the
  other in-plane reciprocal axis also works.
- Other views give sectors that are not fundamental regions, or errors in
  `HSVDirectionKey` (48 cases). They are not in `keys.json`.
- For the valid views, most groups have 2 keys per scheme (+c* and -c*). Cubic
  TSL has 6, because the hue starts at whichever axis is out of screen.
  Monoclinic HSV has up to 8 and triclinic HSV up to 24.

Excluded valid views (listed per group in `keys.json` under `excluded`):

- 1 and -1 TSL with ±a* out of screen: the centre, out-of-screen and the hue
  reference a* coincide. MTEX then uses `orth(center)`, built from its Cartesian
  coordinates.
- 222 and -4 HSV with -a* out of screen: MTEX's split plane is a boundary plane
  of the opposite sign, so the region is empty.
- -4 HSV with +a* out of screen: MTEX's special -4 fold leaves 256 of 500 test
  directions outside its own sector (checked in MATLAB with
  `sR.checkInside`).

### Rhombohedral axes

The R groups 3, -3, 32, 3m and -3m on hexagonal axes are the groups 3, -3, 321,
3m1 and -3m1. On rhombohedral axes (`3R`, `-3R`, `32R`, `3mR`, `-3mR` in
`ipf_keys.py`) every key is the hexagonal one in the other basis. With the
obverse setting, ubi_H = T ubi_R and T = [[1,-1,0],[0,1,-1],[1,1,1]]:

- operations O_R = T^-1 O_H T (acting on ubi; integer),
- directions uvw_R = T^T uvw_H,
- plane normals hkl_R = 3 T^-1 hkl_H, divided by the common factor.

The views keep their hexagonal names (`+c*|+a` looks down c*_H = [111]_R).
The same crystal on either axes gets the same colours, to 1e-9 (tested).

## 7. Other implementations

Measured on the same 500 directions, against the MATLAB keys:

| | same as a MATLAB key | different key |
|---|---|---|
| pymtex 1.0.0 (default view) | everything else, to <=4.2e-4; 23 and m-3 TSL to 6e-3 (its own centre) | HSV 1, 121, 12/m1, 112, 112/m, 23; TSL 112, 11m, 112/m |
| orix 0.15.0 (TSL only, frame x \|\| a*, z \|\| c) | 222 to 4/mmm, 31m, -31m, 622 to 6/mmm, cubic, to 0.005-0.009 (not the correctAngle steps; cause not measured) | -1, all monoclinic, 3, -3, 321, 312, 3m1, -3m1, 6, -6, 6/m |

orix has no HSV key, and no 2/m11, 12/m1, 2mm, m2m or -4m2.

The different keys are listed with their sectors in `ipf_keys.md`. Two examples:

- orix 6/m: wedge [210]-[120] where MTEX's is [100]-[110].
- orix 112/m: boundary through a*, whose uvw ~[1 0.145 0] is not integer. That
  key follows from orix's B convention, not from the lattice. pymtex 112/m has
  vertices ±[010], MTEX ±[100].

## 8. Jumps

Measured on pymtex over 200000 pairs of directions 1e-4 rad apart, before this
table was made:

- TSL has jumps of order 1 for -1, 2/m, 4/m, -3, -3m, 6/m and m-3.
- HSV has them only for -1 and -3. Elsewhere its largest difference is 0.011,
  for example across the edges of 432, which are matched by rotations rather
  than mirrors.

Among the lattice holohedries, TSL jumps for -1, 2/m and -3m (hR), and HSV for
-1 only.

## 9. Drawing the key

The colour is a function on the sphere. The projection used to draw the legend
moves points on the page but does not change any colour: MTEX defaults to equal
area, orix to equal angle (stereographic). The legend is centred on
out-of-screen with east to the right.

`key_image` (drawn by `plot_key`) colours the directions of the sector on a grid
over the unit disk. The polar angle from out-of-screen is theta = 2 atan(R)
(stereographic) or 2 asin(R / sqrt 2) (equal area). A second panel shows the
lower hemisphere, seen from above, when the sector reaches it. The corners are
labelled with their [uvw]. Check: in equal area the sector covers
4 pi / |group folded with| to 5.4 % (101 x 101 pixels, every key).

## 10. Not covered

- maxAngle and the other non-default HSV options.
- Orientation-level checks (ubi plus sample axis against MTEX `ipfHSVKey`).
  Here only the direction step is compared. The orientation step is
  uvw = UB^T r (section 1).

## Files

| file | what |
|---|---|
| `dump_keys.m` | MATLAB: sectors, operations, colours for all groups and views -> `mtex72_keys.json` (59 MB, not for the repository) |
| `check_keys.py` | which views give a fundamental region |
| `make_table.py`, `compare_libs.py` | readable table `ipf_keys.md` with pymtex and orix entries (`libs_rgb.json`) |
| `build_keys.py` | `keys.json` (157 KB; 72 KB without the check colours) |
| `ipf_ref.py`, `check_ref.py` | this file as code, and its check against MATLAB |
| `metric_strain.py`, `strain_fold.py` | section 3 |
| `make_shipped.py` | `ImageD11/ipf_keys.py` and `test/ipf_oracle.json`, R groups added |
| `cmp_quad.py` | `rgb_colors` against `ipf_ref.py` with 1000 and 200000 steps (5.5) |
