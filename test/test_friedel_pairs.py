"""
Friedel pair matching on a simulated scanning-3DXRD peaks_2d (cf_2d) table.

Forward model: point grains off the rotation axis, a beam of finite width,
dty steps symmetric about y0 = 0 and omega frames over 360 degrees. Each
reflection is recorded once per dty row where the grain is in the beam (not
merged over y), on the omega frame containing it. The detector position uses
the true grain position; the geometry is then computed as for measured data
(grain assumed on the axis), so pairs carry realistic g and eta residuals.

Truth: a vertical (omega) pair is -hkl at eta -> 180 - eta, a diagonal (eta)
pair is -hkl at eta -> 180 + eta and a horizontal pair is the same hkl at
eta -> -eta (true eta of the ray, from the forward model).
"""
from __future__ import print_function, division
import os
import sys
import shutil
import tempfile
import logging

import numpy as np
import pytest

from ImageD11 import transform, columnfile, parameters
from ImageD11.sinograms import dataset
from ImageD11 import friedel_pairs as fp

YSTEP = 5.0         # micron
NY = 81             # dty = -200 .. 200
OSTEP = 0.25        # degree
NO = 1440           # 360 degrees
HALF_BEAM = 6.0     # micron: a grain is seen on 2 or 3 dty rows

PARS = dict(distance=150000., y_center=1024., z_center=1024., y_size=75.,
            z_size=75., tilt_x=0., tilt_y=0., tilt_z=0., o11=1, o12=0,
            o21=0, o22=-1, wavelength=0.2846, wedge=0., chi=0., t_x=0.,
            t_y=0., t_z=0., omegasign=1.)
DETECTOR = {k: v for k, v in PARS.items() if k not in ('t_x', 't_y', 't_z')}


def _wrap(a):
    return (a + 180.0) % 360.0 - 180.0


def _fcc_hkls(dsmax, a):
    r = np.arange(-6, 7)
    h = np.array(np.meshgrid(r, r, r, indexing='ij')).reshape(3, -1)
    par = h % 2
    fcc = (par[0] == par[1]) & (par[1] == par[2])
    ds = np.sqrt((h * h).sum(axis=0)) / a
    return h[:, fcc & (ds > 0) & (ds < dsmax)]


def simulate_cf2d(seed=42, ngrains=5, a=4.05, dsmax=1.6):
    rng = np.random.RandomState(seed)
    hkl = _fcc_hkls(dsmax, a)
    cols = {k: [] for k in ('sc', 'fc', 'omega', 'dty', 'sum_intensity',
                            'h', 'k', 'l', 'grain', 'eta_true')}
    ygrid = (np.arange(NY) - NY // 2) * YSTEP
    for gid in range(ngrains):
        U = np.linalg.qr(rng.normal(size=(3, 3)))[0]
        g = np.dot(U / a, hkl)
        sx, sy = rng.uniform(-120, 120, 2)
        Ihkl = 1e4 * np.exp(rng.normal(scale=1.0, size=hkl.shape[1]))
        # Friedel's law: I(hkl) == I(-hkl)
        key = [tuple(sorted([tuple(v), tuple(-v)])) for v in hkl.T]
        first = {}
        for i, k in enumerate(key):
            Ihkl[i] = first.setdefault(k, Ihkl[i])
        tth, etas, omegas = transform.uncompute_g_vectors(g, PARS['wavelength'])
        for eta, om in zip(etas, omegas):
            ok = np.isfinite(om) & np.isfinite(eta) & (np.abs(eta) > 1)
            ok &= (np.abs(np.abs(eta) - 180) > 1)
            w = om[ok] % 360
            wf = (np.floor(w / OSTEP) + 0.5) * OSTEP     # frame centre
            sw, cw = np.sin(np.radians(w)), np.cos(np.radians(w))
            lx = sx * cw - sy * sw
            dty_c = -sx * sw - sy * cw                   # y0 = 0
            for yi, dty in enumerate(ygrid):
                ly = dty - dty_c
                inb = np.abs(ly) < HALF_BEAM
                if not inb.any():
                    continue
                n = inb.sum()
                # ray from the grain position; translations must be scalars
                fc, sc = np.array([transform.compute_xyz_from_tth_eta(
                    np.array([t2]), np.array([e]), np.zeros(1),
                    t_x=x, t_y=y, t_z=0., **DETECTOR)
                    for t2, e, x, y in zip(tth[ok][inb], eta[ok][inb],
                                           lx[inb], ly[inb])])[:, :, 0].T
                prof = np.exp(-0.5 * (ly[inb] / (HALF_BEAM / 2)) ** 2)
                cols['sc'].append(sc + rng.normal(scale=0.1, size=n))
                cols['fc'].append(fc + rng.normal(scale=0.1, size=n))
                cols['omega'].append(wf[inb])
                cols['dty'].append(np.full(n, dty))
                cols['sum_intensity'].append(
                    Ihkl[ok][inb] * prof * rng.uniform(0.9, 1.1, n))
                for j, c in enumerate('hkl'):
                    cols[c].append(hkl[j][ok][inb])
                cols['grain'].append(np.full(n, gid))
                cols['eta_true'].append(eta[ok][inb])
    d = {k: np.concatenate(v) for k, v in cols.items()}
    on = (d['sc'] > 0) & (d['sc'] < 2048) & (d['fc'] > 0) & (d['fc'] < 2048)
    cf = columnfile.colfile_from_dict({k: v[on] for k, v in d.items()})
    cf.parameters = parameters.parameters(**PARS)
    cf.updateGeometry()
    return cf


def make_dataset(folder):
    ds = dataset.DataSet(dataroot=folder, analysisroot=folder,
                         sample='sim', dset='friedel')
    ygrid = (np.arange(NY) - NY // 2) * YSTEP
    ds.dty = np.repeat(ygrid[:, None], NO, axis=1)
    ds.omega = np.repeat(((np.arange(NO) + 0.5) * OSTEP)[None, :], NY, axis=0)
    ds.shape = ds.dty.shape
    ds.guessbins()
    return ds


EXPECTED = {   # pair_type -> (same hkl?, partner true eta from eta)
    'omega': (False, lambda e: 180 - e),
    'eta': (False, lambda e: 180 + e),
    'horizontal': (True, lambda e: -e),
}


def pair_quality(cf, pair_type):
    """(n_correct, n_wrong) pairs against the forward-model truth."""
    same, f = EXPECTED[pair_type]
    i1, i2 = fp.get_pairs(cf, pair_type)
    s = 1 if same else -1
    hkl_ok = ((cf.h[i1] == s * cf.h[i2]) & (cf.k[i1] == s * cf.k[i2]) &
              (cf.l[i1] == s * cf.l[i2]) & (cf.grain[i1] == cf.grain[i2]))
    eta_ok = np.abs(_wrap(cf.eta_true[i2] - f(cf.eta_true[i1]))) < 0.5
    good = hkl_ok & eta_ok
    return int(good.sum()), int((~good).sum())


def n_possible(cf, pair_type):
    """Rows whose true partner reflection is also present (upper bound)."""
    same, f = EXPECTED[pair_type]
    s = 1 if same else -1
    have = set(zip(cf.grain, cf.h, cf.k, cf.l, np.round(_wrap(cf.eta_true), 3)))
    want = zip(cf.grain, s * cf.h, s * cf.k, s * cf.l,
               np.round(_wrap(f(cf.eta_true)), 3))
    return sum(w in have for w in want)


@pytest.fixture(scope='module')
def sim():
    return simulate_cf2d()


@pytest.fixture
def folder():
    d = tempfile.mkdtemp()
    yield d
    shutil.rmtree(d)


def _quiet(func, *args, **kwds):
    logging.getLogger('ImageD11').setLevel(logging.WARNING)
    stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')
    try:
        return func(*args, **kwds)
    finally:
        sys.stdout.close()
        sys.stdout = stdout


def _indexer(cf, ds):
    fpi = fp.FriedelPairIndexer(cf, ds, tol_gv=0.02, tol_eta=0.5)
    fpi.logger.setLevel(logging.WARNING)
    return fpi


@pytest.mark.parametrize('pair_type', sorted(EXPECTED))
def test_global(sim, folder, pair_type):
    cf = sim.copy()
    fpi = _indexer(cf, make_dataset(folder))
    _quiet(fpi.match_friedel_pairs, pair_type=pair_type, doplot=False)
    good, wrong = pair_quality(cf, pair_type)
    npos = n_possible(cf, pair_type)
    print(pair_type, 'global', good, wrong, npos // 2)
    assert wrong == 0
    assert good > 0.85 * npos / 2


CHUNKS = [('omega', 'scans'),
          ('omega', 'frames'),
          ('eta', 'eta_bins'),
          ('horizontal', 'minus_eta_bins')]


@pytest.mark.parametrize('pair_type,chunk_type', CHUNKS)
def test_chunks(sim, folder, pair_type, chunk_type):
    cf = sim.copy()
    fpi = _indexer(cf, make_dataset(folder))
    _quiet(fpi.match_friedel_pairs_by_chunks, chunk_type=chunk_type,
           n_eta_bins=360, n_workers=1, doplot=False)
    good, wrong = pair_quality(cf, pair_type)
    npos = n_possible(cf, pair_type)
    print(pair_type, chunk_type, good, wrong, npos // 2)
    assert wrong == 0
    assert good > 0.85 * npos / 2
