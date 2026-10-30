"""harvest_viewer: see what the harvester did, grain by grain (read only: nothing here changes a result file).

    grain map            click a grain
      -> box sinogram    rows = predicted reflections of that grain, columns = dty (scans); click a row
      -> I(dty) profile  obs (box - bg) vs calc = scale x projection of the grain shape at that omega; click a scan
      -> 1D profiles     s, f, omega profiles (observed, learned) and the dty data used for the intensity
      -> boxes           the (reflection, scan) RE-HARVESTED from the raw frames with the harvester's own kernel:
                         box / annulus on every frame, clipped annulus pixels, other grains' boxes, per-frame
                         box / bg / net, and the totals compared with the stored values

Results are a folder written by ImageD11.sinograms.harvest_pipeline (pass2.h5 or pass1.h5 or harvest.h5, and
grains.h5, shapes.h5, learned.h5 when present) or a single harvest file (Harvester.run(outfile=...)); with no
shapes.h5 the grain shapes are the TensorMap labels of the dataset's grains file.

In a notebook:
    %matplotlib widget
    from ImageD11.nbGui import harvest_viewer as hv
    v = hv.Viewer('dataset.h5', 'phase', {'first harvest': 'run/harvest.h5', 'final': 'run'}, gid=0)
    v.show()

Without the GUI:  src = hv.Results('run', hv.Context('dataset.h5', 'phase'));  g = hv.GrainView(src, gid)
                  r = hv.reharvest(src, g.spot[j], iscan)
"""
from __future__ import print_function
import os

import numpy as np
import h5py

import ImageD11.sinograms.dataset
import ImageD11.parameters
import ImageD11.grain
from ImageD11.sinograms import harvester, profile_fit


# ------------------------------------------------------------------------------------------ data
class Context(object):
    """The dataset, its parameters (phase) and the TensorMap labels of its grains file"""

    def __init__(self, dataset, phase, raw_data=None, tensormap=None):
        self.dsfile, self.phase, self.raw_data = dataset, phase, raw_data
        self.tensormap = tensormap or 'TensorMap_' + phase
        cwd = os.getcwd()
        os.chdir(os.path.dirname(os.path.abspath(dataset)))     # file names in the dataset are relative to it
        try:
            self.ds = ImageD11.sinograms.dataset.load(dataset)
            self.pars = ImageD11.parameters.parameters()
            self.pars.loadparameters(self.ds.parfile, phase_name=phase)
        finally:
            os.chdir(cwd)
        self._labels = None

    def labels(self):
        """grain label image of the TensorMap (reconstruction frame), -1 outside grains"""
        if self._labels is None:
            with h5py.File(self.ds.grainsfile, 'r') as f:
                m = f[self.tensormap]['maps']['labels'][0]
            self._labels = np.flip(np.swapaxes(m, 0, 1), 1).astype(int)
        return self._labels


def _first(folder, names):
    for n in names:
        p = os.path.join(folder, n)
        if os.path.exists(p):
            return p
    return None


class Results(object):
    """One set of results: a harvest file (spots, per_scan, merged) and, when present, the per-grain intensities
    (grains.h5), shapes (shapes.h5), learned profiles (learned.h5) and the harvest holding the observed 1D
    profiles. path: a harvest_pipeline output folder or a harvest file."""

    def __init__(self, path, ctx, name=None):
        self.ctx, self.name = ctx, name or os.path.basename(os.path.normpath(path))
        if os.path.isdir(path):
            self.harvest = _first(path, ('pass2.h5', 'pass1.h5', 'harvest.h5'))
            if self.harvest is None:
                raise ValueError('no harvest file in %s' % path)
            self.grains = _first(path, ('grains.h5',)) if self.harvest.endswith('pass2.h5') else None
            self.shapes = _first(path, ('shapes.h5',)) if self.grains else None
            self.learned = _first(path, ('learned.h5',))
            self.profiles = _first(path, ('pass1.h5',))
        else:
            self.harvest, self.grains, self.shapes, self.learned = path, None, None, None
            self.profiles = path
        self.h5 = h5py.File(self.harvest, 'r')
        f = self.h5
        self.gid = f['gid'][()].astype(int)
        self.members = ([[int(x) for x in m.split()] for m in f['members'][()].astype(str)] if 'members' in f
                        else [[g] for g in self.gid])
        self.k_of = dict((int(g), k) for k, g in enumerate(self.gid))
        self.into = dict((m, int(self.gid[k])) for k, ms in enumerate(self.members) for m in ms)
        sp = dict((k, f['spots'][k][()]) for k in ('grain', 's0', 'f0', 'hs', 'hf', 'ho', 'bw', 'omega'))
        if not np.all(np.diff(sp['grain']) >= 0):
            raise ValueError('spots of a grain are not contiguous in %s' % self.harvest)
        sp['start'] = np.searchsorted(sp['grain'], np.arange(len(self.gid) + 1))
        self.spots = sp
        self.T, self.radius = f['translation'][()], f['radius'][()]
        self.omega_step = float(f.attrs['box_omega_step'])
        self._map, self._blocks, self._harv = None, None, None

    def __repr__(self):
        return 'Results(%s: %s%s%s)' % (self.name, self.harvest, ', ' + self.grains if self.grains else '',
                                        ', ' + self.shapes if self.shapes else '')

    def close(self):
        self.h5.close()

    def index(self, gid):
        if int(gid) not in self.k_of:
            if int(gid) in self.into:
                raise KeyError('grain %d was merged into grain %d' % (gid, self.into[int(gid)]))
            raise KeyError('grain %d is not in %s' % (gid, self.name))
        return self.k_of[int(gid)]

    def _row_blocks(self):
        """per_scan rows are written scan by scan, each block sorted by spot: (block starts, rows of grain k in
        block b = [lo[b, k], lo[b, k + 1]))"""
        if self._blocks is None:
            spot = self.h5['per_scan/spot'][()]
            b0 = np.concatenate([[0], np.flatnonzero(np.diff(spot) < 0) + 1, [len(spot)]])
            st = self.spots['start']
            lo = np.array([b0[b] + np.searchsorted(spot[b0[b]:b0[b + 1]], st) for b in range(len(b0) - 1)])
            self._blocks = lo
        return self._blocks

    def grain(self, gid):
        """dict of (scan, reflection) arrays (iscan, box, bg, var, nframes) and per-reflection columns (refl_*,
        I, sigI) of grain gid; from grains.h5 when there is one (with I_final, the model rows 'model', wI, wvar)"""
        k = self.index(gid)
        a, b = self.spots['start'][k], self.spots['start'][k + 1]
        sp = self.h5['spots']
        z = dict(('refl_' + c, sp[c][a:b]) for c in ('omega', 'tth', 'eta', 'edge', 'touch_mask', 'n_overlap',
                                                     'nbox', 's', 'f', 'hs', 'ho'))
        z['refl_hkl'] = np.column_stack([sp[c][a:b] for c in ('h', 'k', 'l')])
        z['I'], z['sigI'] = self.h5['merged/I'][a:b], self.h5['merged/sigI'][a:b]
        if self.grains:
            with h5py.File(self.grains, 'r') as f:
                g = f[str(int(gid))]
                for c in ('iscan', 'box', 'bg', 'var', 'nframes', 'I_final', 'sigI_final', 'I_sum', 'sigI_sum',
                          'profile_used', 'drop', 'chi2', 'wI', 'wvar', 'overlap', 'omega_edge', 'coverage',
                          'transmission'):
                    if c in g:
                        z[c] = g[c][()]
                z['model'] = g['net'][()]
            z['I'], z['sigI'] = z['I_final'], z['sigI_final']
            return z
        lo = self._row_blocks()
        ps = self.h5['per_scan']
        sl = [slice(lo[i, k], lo[i, k + 1]) for i in range(len(lo)) if lo[i, k + 1] > lo[i, k]]
        col = dict((c, np.concatenate([ps[c][s] for s in sl]) if sl else np.zeros(0))
                   for c in ('spot', 'scan', 'box', 'bg', 'var', 'nframes'))
        iscan = np.unique(col['scan']).astype(int)
        si, ri = np.searchsorted(iscan, col['scan']), col['spot'] - a
        z['iscan'] = iscan
        for c in ('box', 'bg', 'var', 'nframes'):
            arr = np.zeros((len(iscan), b - a), np.float32)
            arr[si, ri] = col[c]
            z[c] = arr
        return z

    def recon(self, gid):
        """grain shape (reconstruction frame): shapes.h5, else the TensorMap labels of its members"""
        if self.shapes:
            with h5py.File(self.shapes, 'r') as f:
                return np.nan_to_num(f[str(int(gid))]['recon'][()].astype(float), nan=0.0, posinf=0.0)
        lab = self.ctx.labels()
        return np.isin(lab, self.members[self.index(gid)]).astype(float)

    def grain_map(self):
        """(rgb, idmap) in the reconstruction frame (plot with origin='lower'). With shapes: each shape scaled
        by its 99.5th percentile, idmap = the largest where it is above 0.2, brightness = that value"""
        if self._map is None:
            if self.shapes:
                top, idmap = None, None
                with h5py.File(self.shapes, 'r') as f:
                    for key in f:
                        r = np.clip(np.nan_to_num(f[key]['recon'][()].astype(np.float32), nan=0.0, posinf=0.0), 0,
                                    None)
                        q = np.percentile(r[r > 0], 99.5) if (r > 0).any() else 1.0
                        r = r / max(q, 1e-30)
                        if top is None:
                            top, idmap = np.zeros_like(r), np.full(r.shape, -1)
                        up = r > top
                        top[up], idmap[up] = r[up], int(key)
                idmap = np.where(top > 0.2, idmap, -1)
                bright = np.clip(top, 0, 1)
            else:
                lab = self.ctx.labels()
                lut = np.full(max(lab.max(), max(self.into)) + 1, -1)
                for m, g in self.into.items():
                    lut[m] = g
                idmap = np.where(lab >= 0, lut[np.maximum(lab, 0)], -1)
                bright = (idmap >= 0).astype(float)
            import colorsys
            ids = np.unique(idmap[idmap >= 0])
            lut = np.full((max(ids.max() if len(ids) else 0, 0) + 1, 3), 0.92)
            for g in ids:
                lut[g] = colorsys.hsv_to_rgb((g * 0.6180339887) % 1.0, 0.7, 1.0)
            rgb = lut[np.maximum(idmap, 0)] * bright[..., None]
            rgb[idmap < 0] = 0.92
            self._map = (np.clip(rgb, 0, 1), idmap)
        return self._map

    def profiles_1d(self, gid, j, z=None):
        """1D profiles of reflection j (local index) of grain gid: observed s, f, omega profiles (net, summed over
        scans), learned profiles (normalised), coordinates relative to the prediction, and the dty data.
        Returns (dict, '') or (None, reason)."""
        if self.profiles is None:
            return None, 'no 1D profiles in %s (harvest without profiles)' % self.name
        k = self.index(gid)
        i = int(self.spots['start'][k] + j)
        with h5py.File(self.profiles, 'r') as f:
            if 'profiles' not in f:
                return None, 'no 1D profiles in %s' % self.profiles
            sp = dict((c, f['spots'][c][i]) for c in ('hs', 'hf', 'ho', 's', 'f', 's0', 'f0', 'h', 'k', 'l'))
            off = f['profiles/offset'][i:i + 2]
            obs = f['profiles/data'][off[0]:off[1]].astype(float)
        lrn, nref, isref = None, 0, False
        if self.learned:
            with h5py.File(self.learned, 'r') as f:
                lo = f['offset'][i:i + 2]
                lrn = f['data'][lo[0]:lo[1]].astype(float)
                nref, isref = int(f['n_ref'][i]), bool(f['reference'][i])
        a, b, c = harvester.profile_slices(sp['hs'], sp['hf'], sp['ho'], 0)
        out = dict(spot=i, hkl=(int(sp['h']), int(sp['k']), int(sp['l'])), n_ref=nref, is_reference=isref)
        for d, sl, h, x0 in (('s', a, sp['hs'], sp['s'] - sp['s0']), ('f', b, sp['hf'], sp['f'] - sp['f0']),
                             ('omega', c, sp['ho'], 0.0)):
            out[d] = (np.arange(2 * h + 1) - h - x0, obs[sl], None if lrn is None else lrn[sl])
        z = self.grain(gid) if z is None else z
        seen = z['nframes'][:, j] > 0
        out['dty'] = dict(dty=self.ctx.ds.ybincens[z['iscan']], seen=seen,
                          net=(z['box'][:, j] - z['bg'][:, j]).astype(float), var=z['var'][:, j].astype(float))
        for c in ('wI', 'wvar', 'model'):
            if c in z:
                out['dty'][c] = z[c][:, j].astype(float)
        for c in ('I_final', 'sigI_final', 'I_sum', 'sigI_sum', 'profile_used', 'chi2', 'drop', 'overlap',
                  'omega_edge', 'coverage', 'transmission', 'refl_touch_mask', 'refl_edge'):
            out[c] = z[c][j] if c in z else None
        if out['I_sum'] is None:
            out['I_sum'], out['sigI_sum'] = float(z['I'][j]), float(z['sigI'][j])
        return out, ''

    # ---------------------------------------------------------------- re-harvest
    def harvester(self):
        """A Harvester set up as the one that wrote the harvest file (mask, distortion, monitor scale, reader,
        background estimator), without prediction: the stored spots table is used."""
        if self._harv is None:
            import scipy.ndimage
            f, ds = self.h5, self.ctx.ds
            grains = [ImageD11.grain.grain(u, translation=t) for u, t in zip(f['ubi'][()], f['translation'][()])]
            datapath = self.ctx.raw_data or f.attrs.get('datapath', None)
            h = harvester.Harvester(ds, grains, self.ctx.pars, ds_max=float(f.attrs['ds_max']), radius=self.radius,
                                    datapath=datapath)
            h._good = np.ones(h.detector_shape(), bool) if h.mask is None else ~h.mask
            h._adj = (np.zeros(h._good.size, np.uint8) if h.mask is None else
                      scipy.ndimage.binary_dilation(h.mask, np.ones((3, 3), bool)).ravel().astype(np.uint8))
            bg = f.attrs.get('box_bg_estimator', 'clipped')
            h.bg_estimator = bg.decode() if isinstance(bg, bytes) else str(bg)
            h.p_clip = float(f.attrs.get('box_p_clip', 1e-3))
            self._harv = h
        return self._harv


# ------------------------------------------------------------------------------------------ one grain
class GrainView(object):
    """One grain of a Results: the stored (reflection, scan) rows plus a calculated I(dty) for every reflection:
    calc = a * P(omega) (and a2 * P + b2), P = projection of the grain shape.

    shape='idmap': the grain shape inside its own region of the grain map; 'recon': the whole shape (clipped at
    floor * max). Arrays are (reflection, scan): box, bg, var, nframes, net, valid, calc (= P, unscaled)."""

    def __init__(self, src, gid, shape='idmap', floor=0.02):
        self.src = src
        ds = src.ctx.ds
        self.gid = int(gid)
        self.k = src.index(gid)
        z = src.grain(gid)
        self.z = z
        self.iscan = z['iscan'].astype(int)
        self.dty = ds.ybincens[self.iscan]
        for name in ('box', 'bg', 'var', 'nframes'):
            setattr(self, name, np.ascontiguousarray(z[name].T.astype(float)))
        self.net = self.box - self.bg
        self.valid = self.nframes > 0
        a, b = src.spots['start'][self.k], src.spots['start'][self.k + 1]
        self.spot = np.arange(a, b)
        assert len(self.spot) == self.box.shape[0], 'grain rows and the harvest spot table disagree'
        self.hkl = z['refl_hkl']
        self.omega, self.tth, self.eta = z['refl_omega'], z['refl_tth'], z['refl_eta']
        self.edge, self.touch = z['refl_edge'], z['refl_touch_mask']
        self.n_overlap, self.nbox = z['refl_n_overlap'], z['refl_nbox']
        self.I, self.sigI = z['I'], z['sigI']
        self.translation = src.T[self.k]
        self.radius = src.radius[self.k]
        self.set_shape(shape, floor)

    def set_shape(self, shape='idmap', floor=0.02):
        ds = self.src.ctx.ds
        rec = self.src.recon(self.gid)
        rec = np.where(rec > floor * rec.max(), rec, 0)
        if shape == 'idmap':
            _, idmap = self.src.grain_map()
            rec = rec * (idmap == self.gid)
        self.shape_mode, self.rho = shape, rec
        P = profile_fit.project(rec, self.omega, ds.y0, ds.ymin, ds.ystep, len(ds.ybincens))
        self.calc = P[:, self.iscan]
        self.fit()
        return self

    def fit(self):
        """per reflection, over its harvested scans: net = a P  and  net = a2 P + b2 (least squares)"""
        m = self.valid.astype(float)
        p, y = self.calc * m, self.net * m
        n, sp_, sy, spp, spy = m.sum(1), p.sum(1), y.sum(1), (p * p).sum(1), (p * y).sum(1)
        with np.errstate(divide='ignore', invalid='ignore'):
            self.a = np.where(spp > 0, spy / spp, 0.0)
            det = n * spp - sp_ * sp_
            ok = det > 1e-12 * np.maximum(n * spp, 1e-300)
            self.a2 = np.where(ok, (n * spy - sp_ * sy) / det, self.a)
            self.b2 = np.where(ok, (spp * sy - sp_ * spy) / det, 0.0)
            sig = np.sqrt(np.maximum(self.var, 0))
            self.sig_row = np.array([np.median(s[v]) if v.any() else np.nan for s, v in zip(sig, self.valid)])
            self.offset_sig = self.b2 / self.sig_row
        self.pmax = self.calc.max(1)
        return self

    def row(self, j):
        """everything about reflection j (local index) for printing"""
        return dict(j=int(j), spot=int(self.spot[j]), hkl=tuple(int(x) for x in self.hkl[j]),
                    omega=float(self.omega[j]), tth=float(self.tth[j]), eta=float(self.eta[j]),
                    I=float(self.I[j]), sigI=float(self.sigI[j]), n_overlap=int(self.n_overlap[j]),
                    edge=bool(self.edge[j]), touch_mask=bool(self.touch[j]), a=float(self.a[j]),
                    a2=float(self.a2[j]), b2=float(self.b2[j]), offset_sig=float(self.offset_sig[j]),
                    nscans=int(self.valid[j].sum()))

    def image(self, mode):
        """(reflection, scan) array for the sinogram panel"""
        with np.errstate(divide='ignore', invalid='ignore'):
            sig = np.sqrt(np.maximum(self.var, 1e-9))
            peak = (self.a * self.pmax)[:, None]
            if mode == 'net / fitted peak':
                im = self.net / np.where(peak > 0, peak, np.nan)
            elif mode == 'net / sigma':
                im = self.net / sig
            elif mode == 'calc (shape)':
                im = self.calc / np.where(self.pmax > 0, self.pmax, np.nan)[:, None]
            elif mode == '(net - a P) / sigma':
                im = (self.net - self.a[:, None] * self.calc) / sig
            elif mode == '(net - a2 P - b2) / sigma':
                im = (self.net - self.a2[:, None] * self.calc - self.b2[:, None]) / sig
            elif mode == 'bg per pixel per frame':
                im = self.bg / (self.nbox[:, None] * self.nframes)
            elif mode == 'nframes':
                im = self.nframes.copy()
            else:
                raise ValueError(mode)
        im[~self.valid] = np.nan
        return im


SINO_MODES = {  # mode: (cmap, vmin, vmax)  None = from data
    'net / fitted peak': ('viridis', -0.3, 1.2),
    'net / sigma': ('viridis', -3, 10),
    'calc (shape)': ('viridis', 0, 1),
    '(net - a P) / sigma': ('RdBu_r', -5, 5),
    '(net - a2 P - b2) / sigma': ('RdBu_r', -5, 5),
    'bg per pixel per frame': ('magma', None, None),
    'nframes': ('cividis', None, None),
}


# colour of (reflection, scan) cells that were not harvested: close to the colour of zero, but not equal to it
BAD_COLOUR = {'viridis': 'black', 'cividis': 'black', 'magma': '#0b2a5b', 'RdBu_r': '#c8c8c8'}


# ------------------------------------------------------------------------------------------ re-harvest
SPOT_COLS = ('grain', 'h', 'k', 'l', 's', 'f', 's0', 'f0', 'hs', 'hf', 'ho', 'bw', 'omega', 'tth', 'eta', 'nbox',
             'nann', 'S1s', 'S1f', 'S2s', 'S2f', 'edge', 'touch_mask', 'n_masked', 'n_overlap')


def spot_row(src, ispot):
    g = src.h5['spots']
    return dict((k, g[k][ispot]) for k in SPOT_COLS)


def stored_row(src, ispot, iscan, z=None):
    """the stored (ispot, iscan) row of the harvest"""
    sp = src.spots
    k = int(sp['grain'][ispot])
    z = src.grain(int(src.gid[k])) if z is None else z
    j = ispot - sp['start'][k]
    c = np.flatnonzero(z['iscan'] == iscan)
    if len(c) == 0:
        return None
    out = dict((n, float(z[n][c[0], j])) for n in ('box', 'bg', 'var', 'nframes'))
    out['net'] = out['box'] - out['bg']
    return out


def neighbours(src, ispot, iscan, frames_omega):
    """other predicted boxes (any grain) overlapping the box + annulus window of ispot in s, f, and the omega
    range shown; with whether their grain is in the beam on this scan"""
    sp = src.spots
    ds = src.ctx.ds
    i = ispot
    rs = sp['hs'][i] + sp['bw'][i]
    rf = sp['hf'][i] + sp['bw'][i]
    step = src.omega_step
    lo, hi = np.min(frames_omega), np.max(frames_omega)
    c = np.flatnonzero((np.abs(sp['s0'] - sp['s0'][i]) <= rs + sp['hs']) &
                       (np.abs(sp['f0'] - sp['f0'][i]) <= rf + sp['hf']) &
                       (sp['omega'] + (sp['ho'] + 0.5) * step >= lo) & (sp['omega'] - (sp['ho'] + 0.5) * step <= hi))
    c = c[c != i]
    g = sp['grain'][c]
    T = src.T[g]
    w = np.radians(sp['omega'][c])
    dty_pred = ds.y0 - T[:, 0] * np.sin(w) - T[:, 1] * np.cos(w)
    inbeam = np.abs(ds.ybincens[iscan] - dty_pred) <= src.radius[g]
    return dict(spot=c, gid=src.gid[g], same_grain=g == sp['grain'][i], inbeam=inbeam,
                s0=sp['s0'][c], f0=sp['f0'][c], hs=sp['hs'][c], hf=sp['hf'][c], ho=sp['ho'][c],
                omega=sp['omega'][c])


def reharvest(src, ispot, iscan, extra=3, z=None):
    """Integrate reflection ispot on scan iscan again from the raw frames, frame by frame, with the
    harvester's kernel and the stored box geometry. Also returns the pixel stack (box + annulus window) of
    the box frames plus `extra` frames either side, the clip threshold per frame and the clipped pixels."""
    h = src.harvester()
    ds = src.ctx.ds
    sp = spot_row(src, ispot)
    om_f = np.asarray(ds.omega[iscan], float)
    order = np.argsort(om_f)
    oms = om_f[order]
    step = float(np.median(np.abs(np.diff(oms))))
    half = (sp['ho'] + 0.5) * step
    lo = int(np.searchsorted(oms, sp['omega'] - half))          # as Harvester.plan_scan
    hi = int(np.searchsorted(oms, sp['omega'] + half))
    box_frames = order[lo:hi]
    show = order[max(lo - extra, 0):min(hi + extra, len(order))]
    T = h.grains[int(sp['grain'])].translation
    w = np.radians(sp['omega'])
    dty_pred = ds.y0 - T[0] * np.sin(w) - T[1] * np.cos(w)
    inbeam = abs(float(np.mean(ds.dty[iscan])) - dty_pred) <= h.radius[int(sp['grain'])]
    s0, f0 = int(sp['s0']), int(sp['f0'])
    hs, hf, bw = int(sp['hs']), int(sp['hf']), int(sp['bw'])
    rs, rf = hs + bw, hf + bw
    nfast = h._good.shape[1]
    good = np.zeros((2 * rs + 1, 2 * rf + 1), bool)
    a0, a1, b0, b1 = max(s0 - rs, 0), min(s0 + rs + 1, h._good.shape[0]), max(f0 - rf, 0), min(f0 + rf + 1, nfast)
    good[a0 - s0 + rs:a1 - s0 + rs, b0 - f0 + rf:b1 - f0 + rf] = h._good[a0:a1, b0:b1]
    inbox = np.zeros_like(good)
    inbox[bw:bw + 2 * hs + 1, bw:bw + 2 * hf + 1] = True
    arr = dict((k, np.array([sp[k]], dtype=np.int64 if k in ('s0', 'f0', 'hs', 'hf', 'bw') else float))
               for k in ('s0', 'f0', 'hs', 'hf', 'bw', 'nbox', 'nann', 'S1s', 'S1f', 'S2s', 'S2f'))
    one = np.zeros(1, np.int64)
    clipped = h.bg_estimator == 'clipped'
    stack = np.zeros((len(show), 2 * rs + 1, 2 * rf + 1), np.float32)
    per = np.zeros((len(show), harvester.NACC))
    thr = np.full(len(show), np.nan)
    used = np.isin(show, box_frames)
    rd = harvester.open_reader(h.scan_file(iscan), ds.limapath, h._good, 'auto')
    try:
        for n, fr in enumerate(show):
            idx, vals = rd(int(fr))
            r, c = idx // nfast, idx % nfast
            k = (np.abs(r - s0) <= rs) & (np.abs(c - f0) <= rf)
            stack[n, r[k] - s0 + rs, c[k] - f0 + rf] = vals[k]              # raw counts (unscaled)
            if sp['nann'] > 0:
                m = stack[n][good & ~inbox].sum() / sp['nann']
                thr[n] = harvester._poisson_upper(m, h.p_clip) if clipped else np.inf
            if used[n]:
                acc = np.zeros((1, harvester.NACC))
                harvester._integrate_frame(idx, vals, nfast, h._adj, one, arr['s0'], arr['f0'], arr['hs'],
                                           arr['hf'], arr['bw'], arr['nbox'], arr['nann'], arr['S1s'],
                                           arr['S1f'], arr['S2s'], arr['S2f'],
                                           np.array([om_f[fr] - sp['omega']]), float(h.scale[iscan, fr]),
                                           clipped, h.p_clip, acc)
                per[n] = acc[0]
    finally:
        rd.close()
    I = harvester._I
    tot = per[used].sum(0)
    tot[I['maxval']] = per[used, I['maxval']].max() if used.any() else 0
    out = dict(ispot=int(ispot), iscan=int(iscan), spot=sp, frames=show, omega=om_f[show], used=used,
               inbeam=bool(inbeam), dty=float(ds.ybincens[iscan]), dty_pred=float(dty_pred),
               scale=np.asarray(h.scale[iscan, show], float), stack=stack, good=good, inbox=inbox, bw=bw,
               threshold=thr, per_frame=dict((k, per[:, I[k]]) for k in harvester.ACC),
               total=dict((k, float(tot[I[k]])) for k in harvester.ACC),
               stored=stored_row(src, ispot, iscan, z), file=h.scan_file(iscan), omega_step=src.omega_step)
    out['total']['net'] = out['total']['box'] - out['total']['bg']
    out['neighbours'] = neighbours(src, ispot, iscan, om_f[show])
    return out


def describe(r):
    """one-paragraph text of a reharvest() result, with the check against the stored row"""
    sp, t, s = r['spot'], r['total'], r['stored']
    lines = ["spot {} grain-index {} hkl ({} {} {})  "
             "tth {:.3f}  eta {:.1f}  omega {:.3f}".format(
                 r['ispot'], int(sp['grain']), int(sp['h']), int(sp['k']), int(sp['l']), sp['tth'],
                 sp['eta'], sp['omega']),
             "scan {} dty {:.2f} (grain centre at this omega {:.2f}; in beam: {})".format(
                 r['iscan'], r['dty'], r['dty_pred'], r['inbeam']),
             "box {}x{} px x {} frames, annulus {} px; "
             "nbox {:.0f}, nann {:.0f}; edge {} touch_mask {} "
             "n_overlap {}".format(
                 2 * int(sp['hs']) + 1, 2 * int(sp['hf']) + 1, int(r['used'].sum()), int(sp['bw']),
                 sp['nbox'], sp['nann'], bool(sp['edge']), bool(sp['touch_mask']), int(sp['n_overlap'])),
             "re-harvest: box {:.1f}  bg {:.1f}  net {:.1f}  sig {:.1f}  "
             "clipped annulus px {:.0f}  nframes {:.0f}".format(
                 t['box'], t['bg'], t['net'], np.sqrt(max(t['var'], 0)), t['nclip'], t['nframes'])]
    nb = r['neighbours']
    inbox = ((np.abs(nb['s0'] - sp['s0']) <= nb['hs'] + sp['hs']) & (np.abs(nb['f0'] - sp['f0']) <= nb['hf'] + sp['hf']) &
             (np.abs(nb['omega'] - sp['omega']) <= (nb['ho'] + sp['ho'] + 1) * r['omega_step']))
    o = inbox & nb['inbeam'] & ~nb['same_grain']
    lines.append("other-grain boxes overlapping the peak box, grain in beam on this scan: {} "
                 "(grains {});  not in beam: {}".format(
                     int(o.sum()), sorted(set(int(x) for x in nb['gid'][o]))[:20],
                     int((inbox & ~nb['inbeam']).sum())))
    if s is None or s['nframes'] == 0:
        lines.append('stored: (not integrated on this scan)')
    else:
        ok = np.isclose(s['box'], t['box'], rtol=1e-4, atol=0.5) and np.isclose(s['bg'], t['bg'], rtol=1e-4, atol=0.5)
        lines.append("stored    : box {:.1f}  bg {:.1f}  net {:.1f}  sig {:.1f}  "
                     "nframes {:.0f}   -> {}".format(
                         s['box'], s['bg'], s['net'], np.sqrt(max(s['var'], 0)), s['nframes'],
                         'MATCH' if ok else 'DIFFERENT'))
    return '\n'.join(lines)


# ------------------------------------------------------------------------------------------ 1D profiles
def draw_profiles(fig, p, why='', scan_dty=None):
    """the 1D profiles of a Results.profiles_1d() result: s, f, omega and dty"""
    fig.clf()
    if p is None:
        ax = fig.add_subplot(111)
        ax.axis('off')
        ax.text(0.5, 0.5, why, ha='center', va='center', fontsize=11)
        return fig
    final = p['I_final'] is not None and np.isfinite(p['I_final'])
    axs = fig.subplots(1, 4, gridspec_kw=dict(width_ratios=[1, 1, 1.2, 2]))
    for ax, d, lab in zip(axs[:3], ('s', 'f', 'omega'), ('s - s_pred (px)', 'f - f_pred (px)', 'omega - omega_pred (frames)')):
        x, obs, lrn = p[d]
        ax.plot(x, obs, 'ko', ms=3, label='observed (net, all scans)')
        if lrn is not None:
            ax.plot(x, lrn * p['I_sum'], 'r-', lw=1.5, label='learned x I_sum')
            if final:
                ax.plot(x, lrn * p['I_final'], 'b--', lw=1, label='learned x I_final')
        ax.axhline(0, color='k', lw=0.5)
        ax.set_xlabel(lab, fontsize=8)
        ax.tick_params(labelsize=7)
    axs[0].legend(fontsize=6, loc='upper left')
    axs[0].set_ylabel('counts', fontsize=8)
    d = p['dty']
    m = d['seen']
    ax = axs[3]
    ax.errorbar(d['dty'][m], d['net'][m], np.sqrt(np.maximum(d['var'][m], 0)), fmt='o', ms=3, color='0.5', lw=0.6,
                label='box sum per scan')
    if 'wI' in d:
        ax.errorbar(d['dty'][m], d['wI'][m], np.sqrt(np.maximum(d['wvar'][m], 0)), fmt='s', ms=3, color='tab:blue',
                    lw=0.6, label='profile-weighted per scan')
    if 'model' in d:
        ax.plot(d['dty'][m], d['model'][m], 'r-', lw=1.5, label='rows used for I (fit or data)')
    ax.axhline(0, color='k', lw=0.5)
    if scan_dty is not None:
        ax.axvline(scan_dty, color='r', lw=0.6)
    ax.set_xlabel('dty', fontsize=8)
    ax.tick_params(labelsize=7)
    ax.legend(fontsize=6)
    t = "spot {} hkl {}: I_sum {:.1f} +- {:.1f}".format(p['spot'], p['hkl'], p['I_sum'], p['sigI_sum'])
    if final:
        used = 'dty profile fit (weak)' if p['profile_used'] else 'box sum'
        t += "  ->  I_final {:.1f} +- {:.1f}  [{}; chi2 {:.2f}".format(p['I_final'], p['sigI_final'], used, p['chi2'])
        if p['transmission'] is not None:
            t += "; transmission {:.3f}; dty coverage {:.3f}".format(p['transmission'], p['coverage'])
        if p['drop']:
            why = [n for n, c in (('contaminated', 'overlap'), ('omega range end', 'omega_edge'),
                                  ('touches mask', 'refl_touch_mask'), ('detector edge', 'refl_edge')) if p[c]]
            t += '; DROPPED (%s)' % (', '.join(why) or 'dty coverage')
        t += ']'
    if p['n_ref']:
        t += "  learned from {} strong refs{}".format(p['n_ref'], '; is a reference' if p['is_reference'] else '')
    fig.suptitle(t, fontsize=8)
    return fig


# ------------------------------------------------------------------------------------------ the GUI
class Viewer(object):
    """ipympl + ipywidgets viewer. Needs %matplotlib widget. Clicks act only when no zoom/pan tool is active.
    sources: {name: results folder or harvest file} (or a list), shown in the 'results' menu."""

    def __init__(self, dataset, phase, sources, gid=None, raw_data=None, map_kind='id', shape='idmap'):
        self.ctx = Context(dataset, phase, raw_data=raw_data)
        if not isinstance(sources, dict):
            sources = dict((os.path.basename(os.path.normpath(p)), p) for p in sources)
        self.paths = sources
        self.cache = {}
        self.source = list(sources)[0]
        self.src = self.results(self.source)
        self.gid = int(self.src.gid[0]) if gid is None else int(gid)
        self.shape = shape
        self.g = None
        self.j = None          # selected reflection (local index)
        self.c = None          # selected scan column
        self.rows = None       # displayed row order (local reflection indices)
        self.r = None

    def results(self, name):
        if name not in self.cache:
            self.cache[name] = Results(self.paths[name], self.ctx, name=name)
        return self.cache[name]

    # ---------------------------------------------------------------- layout
    def show(self):
        import matplotlib.pyplot as plt
        import ipywidgets as W
        from IPython.display import display
        with plt.ioff():
            self.fig = plt.figure(figsize=(13, 9.5))
            self.fig2 = plt.figure(figsize=(13, 8))
            self.fig3 = plt.figure(figsize=(13, 3.2))
        for f in (self.fig, self.fig2, self.fig3):
            f.canvas.header_visible = False
            f.canvas.toolbar_position = 'left'
        gs = self.fig.add_gridspec(2, 2, width_ratios=[1, 1.25], height_ratios=[1.25, 1], hspace=0.3, wspace=0.25)
        self.ax_map = self.fig.add_subplot(gs[0, 0])
        self.ax_sino = self.fig.add_subplot(gs[0, 1])
        self.ax_prof = self.fig.add_subplot(gs[1, 1])
        self.ax_off = self.fig.add_subplot(gs[1, 0])
        self.w_gid = W.IntText(value=self.gid, description='grain', layout=W.Layout(width='150px'))
        self.w_go = W.Button(description='load', layout=W.Layout(width='60px'))
        self.w_src = W.Dropdown(options=list(self.paths), value=self.source, description='results',
                                layout=W.Layout(width='240px'))
        self.w_shape = W.Dropdown(options=['idmap', 'recon'], value=self.shape, description='shape',
                                  layout=W.Layout(width='170px'))
        self.w_mode = W.Dropdown(options=list(SINO_MODES), value='net / fitted peak', description='sinogram',
                                 layout=W.Layout(width='330px'))
        self.w_sort = W.Dropdown(options=['omega', 'tth', 'I/sigI', 'offset/sigma'], value='omega',
                                 description='sort', layout=W.Layout(width='190px'))
        self.w_tth = W.FloatRangeSlider(value=(0, 60), min=0, max=60, step=0.1, description='tth',
                                        continuous_update=False, layout=W.Layout(width='380px'))
        self.w_hide = W.Checkbox(value=False, description='hide edge / touch_mask', indent=False)
        self.w_extra = W.IntSlider(value=3, min=0, max=10, description='extra frames', continuous_update=False)
        self.w_log = W.Checkbox(value=False, description='log frames', indent=False)
        self.w_allnb = W.Checkbox(value=False, description='boxes of grains not in beam', indent=False)
        self.w_prev = W.Button(description='< row')
        self.w_next = W.Button(description='row >')
        self.w_sprev = W.Button(description='< scan')
        self.w_snext = W.Button(description='scan >')
        self.w_status = W.HTML()
        self.w_text = W.Output(layout=W.Layout(border='1px solid #ccc'))
        self.w_go.on_click(lambda _: self.load(self.w_gid.value))
        self.w_src.observe(lambda ch: self.switch(ch['new']), 'value')
        self.w_shape.observe(lambda _: self.load(self.gid, reshape=True), 'value')
        for w in (self.w_mode, self.w_sort, self.w_tth, self.w_hide):
            w.observe(lambda _: self.draw_sino(), 'value')
        for w in (self.w_extra, self.w_log, self.w_allnb):
            w.observe(lambda _: self.select_scan(self.c), 'value')
        self.w_prev.on_click(lambda _: self.step_row(-1))
        self.w_next.on_click(lambda _: self.step_row(1))
        self.w_sprev.on_click(lambda _: self.select_scan(None if self.c is None else self.c - 1))
        self.w_snext.on_click(lambda _: self.select_scan(None if self.c is None else self.c + 1))
        self.fig.canvas.mpl_connect('button_press_event', self.on_click)
        ui = W.VBox([W.HBox([self.w_src, self.w_gid, self.w_go, self.w_shape, self.w_status]),
                     W.HBox([self.w_mode, self.w_sort, self.w_tth, self.w_hide]),
                     self.fig.canvas,
                     W.HBox([self.w_prev, self.w_next, self.w_sprev, self.w_snext, self.w_extra, self.w_log, self.w_allnb]),
                     self.fig3.canvas, self.w_text, self.fig2.canvas])
        display(ui)
        self.draw_map()
        self.load(self.gid)

    def switch(self, name):
        """show another set of results; keeps the grain, or goes to the grain it was merged into"""
        self.status('opening {} ...'.format(name))
        self.source = name
        self.src = self.results(name)
        gid = self.src.into.get(int(self.gid), int(self.gid))
        self.g = None
        self.draw_map()
        self.load(gid)

    def status(self, s):
        self.w_status.value = '<b>{}</b>'.format(s)

    # ---------------------------------------------------------------- events
    def on_click(self, ev):
        tb = getattr(self.fig.canvas, 'toolbar', None)
        if tb is not None and getattr(tb, 'mode', '') not in ('', None):
            return
        if ev.xdata is None:
            return
        if ev.inaxes is self.ax_map:
            _, idmap = self.src.grain_map()
            i, j = int(round(ev.ydata)), int(round(ev.xdata))
            if 0 <= i < idmap.shape[0] and 0 <= j < idmap.shape[1] and idmap[i, j] >= 0:
                self.load(int(idmap[i, j]))
            else:
                self.status('no grain there')
        elif ev.inaxes is self.ax_sino and self.rows is not None and len(self.rows):
            r = int(np.clip(round(ev.ydata), 0, len(self.rows) - 1))
            self.select_row(self.rows[r], col=self.col_of_dty(ev.xdata))
        elif ev.inaxes is self.ax_off and self.g is not None:
            g = self.g
            x, y = g.tth, g.offset_sig
            ok = np.isfinite(y) & self.row_filter()
            if ok.any():
                xr = np.ptp(self.ax_off.get_xlim())
                yr = np.ptp(self.ax_off.get_ylim())
                d = ((x - ev.xdata) / xr) ** 2 + ((y - ev.ydata) / yr) ** 2
                d[~ok] = np.inf
                self.select_row(int(np.argmin(d)))
        elif ev.inaxes is self.ax_prof and self.j is not None:
            self.select_scan(self.col_of_dty(ev.xdata))

    def col_of_dty(self, x):
        return int(np.argmin(np.abs(self.g.dty - x)))

    # ---------------------------------------------------------------- grain
    def load(self, gid, reshape=False):
        self.status('loading grain {} ...'.format(gid))
        try:
            if reshape and self.g is not None:
                self.g.set_shape(self.w_shape.value)
            else:
                self.g = GrainView(self.src, gid, shape=self.w_shape.value)
        except Exception as e:
            self.status('grain {}: {!r}'.format(gid, e))
            return
        self.gid = int(gid)
        self.w_gid.value = self.gid
        self.j = None
        self.c = None
        self.draw_map()
        self.draw_sino()
        self.draw_offsets()
        g = self.g
        strong = np.argsort(np.where(g.valid.any(1), g.I / np.maximum(g.sigI, 1e-9), -np.inf))[::-1]
        self.select_row(int(strong[0]))
        self.status('grain {}: {} reflections x {} scans; shape={}'.format(
            gid, len(g.spot), len(g.iscan), g.shape_mode))

    def draw_map(self):
        ax = self.ax_map
        ax.clear()
        rgb, idmap = self.src.grain_map()
        ax.imshow(rgb, origin='lower', interpolation='nearest')
        if self.g is not None:
            ax.contour(idmap == self.gid, levels=[0.5], colors='k', linewidths=1.2)
            ax.contour(self.g.rho > 0, levels=[0.5], colors='w', linewidths=0.6, linestyles=':')
        ax.set_title('{}: grain map - click a grain; grain {}'.format(self.source, self.gid), fontsize=9, loc='left')
        ax.set_xticks([])
        ax.set_yticks([])
        self.fig.canvas.draw_idle()

    def row_filter(self):
        g = self.g
        lo, hi = self.w_tth.value
        ok = (g.tth >= lo) & (g.tth <= hi) & g.valid.any(1)
        if self.w_hide.value:
            ok &= ~(g.edge | g.touch)
        return ok

    def draw_sino(self):
        if self.g is None:
            return
        import matplotlib
        g = self.g
        ok = np.flatnonzero(self.row_filter())
        key = {'omega': g.omega, 'tth': g.tth, 'I/sigI': g.I / np.maximum(g.sigI, 1e-9),
               'offset/sigma': np.nan_to_num(g.offset_sig, nan=-np.inf)}[self.w_sort.value]
        self.rows = ok[np.argsort(key[ok], kind='stable')]
        mode = self.w_mode.value
        im = g.image(mode)[self.rows]
        cmap, vmin, vmax = SINO_MODES[mode]
        if vmin is None and np.isfinite(im).any():
            vmin, vmax = np.nanpercentile(im, [1, 99])
        ax = self.ax_sino
        ax.clear()
        d = g.dty
        cm = matplotlib.colormaps[cmap].copy()
        cm.set_bad(BAD_COLOUR.get(cmap, 'k'))                  # not harvested: near the zero colour, but distinct
        ax.imshow(im, aspect='auto', origin='lower', interpolation='nearest', cmap=cm, vmin=vmin, vmax=vmax,
                  extent=(d[0] - 0.05, d[-1] + 0.05, -0.5, len(self.rows) - 0.5))
        ax.set_xlabel('dty')
        ax.set_ylabel('reflection (sorted by {})'.format(self.w_sort.value))
        ax.set_title('grain {}: {}   [{:.3g}, {:.3g}]'.format(g.gid, mode, vmin, vmax), fontsize=9, loc='left')

        def fmt(x, y):
            r = int(round(y))
            if not (0 <= r < len(self.rows)):
                return ''
            j = self.rows[r]
            return ('dty {:.2f}  hkl {}  omega {:.2f}  tth {:.2f}  '
                    '{} {:.3g}'.format(x, tuple(g.hkl[j]), g.omega[j], g.tth[j], self.w_sort.value, key[j]))
        ax.format_coord = fmt
        self.sel_h = ax.axhline(0, color='r', lw=0.8)
        self.sel_v = ax.axvline(d[0], color='r', lw=0.8)
        ax.set_xlim(d[0] - 0.05, d[-1] + 0.05)
        ax.set_ylim(-0.5, max(len(self.rows), 1) - 0.5)
        self.mark_selection()
        self.fig.canvas.draw_idle()

    def draw_offsets(self):
        g = self.g
        ax = self.ax_off
        ax.clear()
        ok = np.isfinite(g.offset_sig) & g.valid.any(1)
        snr = np.log10(np.clip(g.I / np.maximum(g.sigI, 1e-9), 0.1, 1e4))
        sc = ax.scatter(g.tth[ok], g.offset_sig[ok], c=snr[ok], s=3, cmap='viridis', vmin=-1, vmax=3)
        ax.axhline(0, color='k', lw=0.5)
        if ok.sum() > 50:                    # running median vs tth
            b = np.quantile(g.tth[ok], np.linspace(0, 1, 21))
            ii = np.digitize(g.tth[ok], b[1:-1])
            med = [np.median(g.offset_sig[ok][ii == k]) for k in range(20)]
            ax.plot(0.5 * (b[1:] + b[:-1]), med, 'r-', lw=1.5, label='median')
            ax.legend(fontsize=8, loc='upper left')
        ax.set_ylim(*np.clip(np.nanpercentile(g.offset_sig[ok], [0.5, 99.5]), -50, 50) if ok.any() else (-5, 5))
        ax.set_xlabel('tth (deg)')
        ax.set_ylabel('offset b2 / sigma')
        ax.set_title('per reflection: constant offset in net = a2 P + b2 (colour log10 I/sigI); click a point',
                     fontsize=8, loc='left')
        self.sel_pt, = ax.plot([], [], 'o', mfc='none', mec='r', ms=9)
        if not hasattr(self, '_cb'):
            self._cb = self.fig.colorbar(sc, ax=ax, pad=0.01)
        else:
            self._cb.update_normal(sc)
        self.fig.canvas.draw_idle()

    # ---------------------------------------------------------------- reflection
    def mark_selection(self):
        if self.j is None or self.rows is None:
            return
        r = np.flatnonzero(self.rows == self.j)
        self.sel_h.set_ydata([r[0] if len(r) else np.nan] * 2)
        if self.c is not None:
            self.sel_v.set_xdata([self.g.dty[self.c]] * 2)
        if hasattr(self, 'sel_pt'):
            self.sel_pt.set_data([self.g.tth[self.j]], [self.g.offset_sig[self.j]])

    def step_row(self, d):
        if self.j is None or self.rows is None or len(self.rows) == 0:
            return
        r = np.flatnonzero(self.rows == self.j)
        r = int(np.clip((r[0] if len(r) else 0) + d, 0, len(self.rows) - 1))
        self.select_row(self.rows[r])

    def select_row(self, j, col=None):
        g = self.g
        self.j = int(j)
        if col is None:
            v = np.where(g.valid[j], g.calc[j], -1)
            col = int(np.argmax(v))
        self.c = col
        ax = self.ax_prof
        ax.clear()
        m = g.valid[j]
        sig = np.sqrt(np.maximum(g.var[j], 0))
        d = g.dty
        ax.plot(d[m], g.box[j][m], '.', color='0.6', ms=3, label='box (sum)')
        ax.plot(d[m], g.bg[j][m], '-', color='tab:blue', lw=1, label='bg (annulus x nbox)')
        ax.errorbar(d[m], g.net[j][m], sig[m], fmt='o', color='k', ms=3, lw=0.8, label='net = box - bg (obs)')
        ax.plot(d, g.a[j] * g.calc[j], '-', color='tab:red', lw=1.5, label='calc  a P  (a={:.3g})'.format(g.a[j]))
        ax.plot(d, g.a2[j] * g.calc[j] + g.b2[j], '--', color='tab:orange', lw=1.2,
                label='a2 P + b2  (b2={:.3g} = {:.2f} sigma)'.format(g.b2[j], g.offset_sig[j]))
        ax.axhline(0, color='k', lw=0.5)
        ax.axvspan(d[col] - 0.05, d[col] + 0.05, color='r', alpha=0.15)
        ax.set_xlabel('dty')
        ax.set_ylabel('counts (monitor scaled)')
        if m.any():
            ax.set_xlim(d[m].min() - 0.3, d[m].max() + 0.3)
        fl = ' EDGE' if g.edge[j] else ''
        fl += ' MASK' if g.touch[j] else ''
        ax.set_title('hkl {}  omega {:.2f}  tth {:.2f}  '
                     'eta {:.1f}  I/sig {:.1f}  n_overlap {}{}'
                     '  - click a scan'.format(
                         tuple(int(x) for x in g.hkl[j]), g.omega[j], g.tth[j], g.eta[j],
                         g.I[j] / max(g.sigI[j], 1e-9), g.n_overlap[j], fl), fontsize=8, loc='left')
        ax.legend(fontsize=7, loc='upper right')
        self.mark_selection()
        self.fig.canvas.draw_idle()
        try:
            p, why = self.src.profiles_1d(self.gid, self.j, g.z)
        except Exception as e:
            p, why = None, 'profiles: {!r}'.format(e)
        self.p = p
        draw_profiles(self.fig3, p, why, scan_dty=g.dty[col])
        self.fig3.canvas.draw_idle()
        self.select_scan(col)

    # ---------------------------------------------------------------- scan / boxes
    def select_scan(self, col):
        if col is None or self.j is None:
            return
        g = self.g
        col = int(np.clip(col, 0, len(g.iscan) - 1))
        self.c = col
        self.mark_selection()
        for l in self.ax_prof.patches[:]:
            l.remove()
        self.ax_prof.axvspan(g.dty[col] - 0.05, g.dty[col] + 0.05, color='r', alpha=0.15)
        self.fig.canvas.draw_idle()
        self.status('re-harvesting spot {} on scan {} from the raw frames ...'.format(g.spot[self.j], g.iscan[col]))
        try:
            self.r = reharvest(self.src, int(g.spot[self.j]), int(g.iscan[col]), extra=self.w_extra.value, z=g.z)
        except Exception as e:
            self.status('re-harvest failed: {!r}'.format(e))
            return
        self.w_text.clear_output()
        with self.w_text:
            print(describe(self.r))
            print("calc at this scan: a P = {:.1f}   "
                  "a2 P + b2 = {:.1f}   "
                  "(obs net {:.1f})".format(
                      g.a[self.j] * g.calc[self.j, col], g.a2[self.j] * g.calc[self.j, col] + g.b2[self.j],
                      g.net[self.j, col]))
        draw_boxes(self.fig2, self.r, log=self.w_log.value, all_neighbours=self.w_allnb.value)
        self.fig2.canvas.draw_idle()
        self.status('grain {}  reflection {}  scan {} (dty {:.2f})'.format(self.gid, self.j, g.iscan[col], g.dty[col]))


def draw_boxes(fig, r, log=False, ncols=None, all_neighbours=False):
    """Frames of a reharvest() result: box (red), annulus (cyan), clipped annulus px (x), other boxes
    (dashed: orange other grain in beam, grey not in beam, green same grain); then per-frame sums."""
    from matplotlib.patches import Rectangle
    from matplotlib.colors import LogNorm, Normalize
    fig.clf()
    st = r['stack']
    n = len(st)
    ncols = ncols or min(n, 12)
    nr = int(np.ceil(n / ncols))
    gs = fig.add_gridspec(nr + 1, ncols, height_ratios=[1] * nr + [1.3], hspace=0.35, wspace=0.08)
    sp = r['spot']
    bw, hs, hf = r['bw'], int(sp['hs']), int(sp['hf'])
    rs, rf = hs + bw, hf + bw
    sbox = st[:, r['inbox']]
    vmax = max(float(np.percentile(sbox, 99.5)) if sbox.size else 1, 1)
    norm = LogNorm(vmin=0.5, vmax=max(vmax, 1)) if log else Normalize(0, vmax)
    nb = r['neighbours']
    step = r['omega_step']
    ann = r['good'] & ~r['inbox']
    for i in range(n):
        ax = fig.add_subplot(gs[i // ncols, i % ncols])
        img = st[i] if not log else np.where(st[i] > 0, st[i], np.nan)
        ax.imshow(img, norm=norm, cmap='viridis', interpolation='nearest', origin='upper')
        bad = ~r['good']
        if bad.any():
            ax.imshow(np.ma.masked_where(~bad, np.ones_like(st[i])), cmap='gray', vmin=0, vmax=2, alpha=0.8)
        if r['used'][i]:
            ax.add_patch(Rectangle((bw - 0.5, bw - 0.5), 2 * hf + 1, 2 * hs + 1, fill=False, ec='r', lw=1.2))
            ax.add_patch(Rectangle((-0.5, -0.5), 2 * rf + 1, 2 * rs + 1, fill=False, ec='c', lw=1.2))
            t = r['threshold'][i]
            if np.isfinite(t):
                yy, xx = np.nonzero(ann & (st[i] > t))
                ax.plot(xx, yy, 'x', color='w', ms=3, mew=0.8)
        ax.plot(sp['f'] - sp['f0'] + rf, sp['s'] - sp['s0'] + rs, '+', color='w', ms=7)
        for k in range(len(nb['spot'])):
            if abs(r['omega'][i] - nb['omega'][k]) > (nb['ho'][k] + 0.5) * step:
                continue
            if not (all_neighbours or nb['inbeam'][k] or nb['same_grain'][k]):
                continue
            col = 'lime' if nb['same_grain'][k] else ('orange' if nb['inbeam'][k] else '0.6')
            x0 = nb['f0'][k] - sp['f0'] + rf - nb['hf'][k] - 0.5
            y0 = nb['s0'][k] - sp['s0'] + rs - nb['hs'][k] - 0.5
            ax.add_patch(Rectangle((x0, y0), 2 * nb['hf'][k] + 1, 2 * nb['hs'][k] + 1, fill=False, ec=col, lw=0.9,
                                   ls='--'))
            ax.text(x0 + 0.5, y0 + 0.5, str(int(nb['gid'][k])), color=col, fontsize=6, va='top', clip_on=True)
        ax.set_xlim(-0.5, 2 * rf + 0.5)
        ax.set_ylim(2 * rs + 0.5, -0.5)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title("{:.2f}".format(r['omega'][i]) + (' *' if r['used'][i] else ''), fontsize=7,
                     color='r' if r['used'][i] else 'k')
    ax = fig.add_subplot(gs[nr, :])
    pf = r['per_frame']
    u = r['used']
    x = r['omega']
    raw_box = np.array([s[r['inbox'] & r['good']].sum() for s in st]) * r['scale']
    raw_ann = np.array([s[ann].sum() for s in st]) / max(float(sp['nann']), 1) * float(sp['nbox']) * r['scale']
    ax.plot(x, raw_box, 'o-', color='0.6', ms=3, lw=0.8, label='box sum (all frames)')
    ax.plot(x, raw_ann, 's:', color='tab:cyan', ms=3, lw=0.8, label='unclipped annulus mean x nbox')
    ax.plot(x[u], pf['bg'][u], 's-', color='tab:blue', ms=3, lw=1, label='bg used (clipped annulus) x nbox')
    ax.plot(x[u], (pf['box'] - pf['bg'])[u], 'o-', color='k', ms=4, lw=1.2, label='net per frame')
    ax.axhline(0, color='k', lw=0.5)
    ax.axvspan(sp['omega'] - (sp['ho'] + 0.5) * step, sp['omega'] + (sp['ho'] + 0.5) * step, color='r', alpha=0.07)
    ax.axvline(sp['omega'], color='r', lw=0.6)
    ax.set_xlabel('omega (deg)  [red band = harvester omega window, * frames used]')
    ax.set_ylabel('counts')
    ax.legend(fontsize=7, ncol=4, loc='upper left')
    t = r['total']
    s = r['stored']
    chk = '' if s is None else "   stored box {:.0f} bg {:.0f} net {:.0f}".format(s['box'], s['bg'], s['net'])
    fig.suptitle("spot {} ({} {} {})  scan {} dty {:.2f}"
                 "  in beam {}   re-harvest box {:.0f} bg {:.0f} net {:.0f} "
                 "sig {:.0f}{}\n[red box / cyan annulus / x clipped / dashed: other boxes "
                 "orange=other grain in beam, grey=not in beam, green=same grain]".format(
                     r['ispot'], int(sp['h']), int(sp['k']), int(sp['l']), r['iscan'], r['dty'], r['inbeam'],
                     t['box'], t['bg'], t['net'], np.sqrt(max(t['var'], 0)), chk), fontsize=8)
    return fig
