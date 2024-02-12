#!/usr/bin/env python
# coding: utf-8

# Try to build the point-by-point mapping code ...
from __future__ import print_function, division

import multiprocessing
method = multiprocessing.set_start_method("forkserver")
from ImageD11 import cImageD11  # NO!!! patches to forkserver

import sys, os
from multiprocessing.managers import SharedMemoryManager
import time, random
import numpy as np
import numba
import pprint
from ImageD11 import sym_u, unitcell, parameters
import ImageD11.sinograms.dataset
import ImageD11.indexing


@numba.njit(boundscheck=True)
def hkluniq(ubi, gx, gy, gz, eta, m, tol, hmax):
    """count uniq peaks - move to cImageD11 in the future"""
    # index from -hmax to hmax
    hcount = np.zeros((2 * hmax + 1, 2 * hmax + 1, 2 * hmax + 1, 2), dtype=np.int32)
    tcount = 0
    for i in range(len(gx)):
        if ~m[i]:
            continue
        H = ubi[0][0] * gx[i] + ubi[0][1] * gy[i] + ubi[0][2] * gz[i]
        K = ubi[1][0] * gx[i] + ubi[1][1] * gy[i] + ubi[1][2] * gz[i]
        L = ubi[2][0] * gx[i] + ubi[2][1] * gy[i] + ubi[2][2] * gz[i]
        ih = np.round(H)
        ik = np.round(K)
        il = np.round(L)
        if eta[i] > 0:
            s = 1
        else:
            s = 0
        dh2 = (H - ih) ** 2 + (K - ik) ** 2 + (L - il) ** 2
        if dh2 < tol * tol:
            # fixme sign(yl)
            if abs(ih) < hmax and abs(ik) < hmax and abs(il) < hmax:
                hcount[int(ih) + hmax, int(ik) + hmax, int(il) + hmax, s] = 1
            tcount += 1
    ucount = 0
    for i in range(hcount.size):
        if hcount.flat[i] > 0:
            ucount += 1
    return tcount, ucount


def dtymask(i, j, cosomega, sinomega, dtyi, step, y0):
    """
    Selects the peaks using integer bins in dty.
    Needs a UI / plot to see what is happening (map the i,j to the sinogram and recon)
    """
    dty_icalc = np.round((step * (i * cosomega + j * sinomega) + y0) / step).astype(int)
    return dtyi == dty_icalc


def idxpoint(
    i,
    j,
    isel,
    cosomega,
    sinomega,
    dtyi,
    gx,
    gy,
    gz,
    eta,
    ystep=2.00,
    y0=0.0,
    minpks=1000,
    hkl_tol=0.1,
    ds_tol=0.005,
    cosine_tol=np.cos(np.radians(90 - 0.1)),
    forgen=None,
    hmax=-1,  # auto
    uniqcut=0.75,
):
    """
    Indexing function called at one point in space.
    Selects peaks from the sinogram and attempts to index.
    More of this could be indexing.py I guess.

    i,j = integer step from the rotation axis center

    Effectively a columnfile:
        isel = column of rings to use
        cosomega = column of cosine omega
        sinomega = column of sin omega
        dtyi = column of dty on integer bins
        gx/gy/gz/eta = peak co-ordinates

    ystep, y0 = compute the dtyi bins to select peaks (dtymask function above)

    forgen = the rings to use for INDEXING / orientation searching (a speedup effect)

    minpks = number of peaks to find with rings
    hkl_tol = ImageD11 tolerance of indexed == | h-int(h) | < hkl_tol
    cosine_tol = ImageD11 tolerance for finding pairs of peaks

    uniqcut = unique peaks cutoff (uniq> ucut * max are returned)

    returns a list of :
        (npks, nuniq, ubi)
    """
    cImageD11.cimaged11_omp_set_num_threads(1)
    # args    isel, cosomega, sinomega, dtyi, gx, gy, gz
    mr = isel
    m = dtymask(i, j, cosomega, sinomega, dtyi, ystep, y0)
    # the arrays to use for indexing (this point in space)
    igx = gx[m & mr]
    igy = gy[m & mr]
    igz = gz[m & mr]
    ind = ImageD11.indexing.indexer(
        unitcell=ucglobal,
        wavelength=parglobal.get("wavelength"),
        gv=np.array((igx, igy, igz)).T,
    )
    ind.minpks = minpks
    ind.hkl_tol = hkl_tol
    ind.cosine_tol = cosine_tol  # degrees
    ind.ds_tol = ds_tol
    ind.assigntorings()
    for ind.ring_1 in forgen:
        for ind.ring_2 in forgen:
            ind.find()
            ind.scorethem()
    if ImageD11.indexing.loglevel <= 3:
        sys.stdout.flush()
    global symglobal
    ind.ubis = [sym_u.find_uniq_u(ubi, symglobal) for ubi in ind.ubis]
    # count the unique peaks per grain
    uniqs = np.empty((len(ind.ubis), 2), int)
    for i, ubi in enumerate(ind.ubis):
        # pass in the dty mask for peak counting with all peaks
        uniqs[i] = hkluniq(ubi, gx, gy, gz, eta, m, hkl_tol, hmax)
    if len(ind.ubis) == 0:
        return [
            (0, 0, np.eye(3)),
        ]
    if len(ind.ubis) == 1:
        return [
            (uniqs[0][0], uniqs[0][1], ind.ubis[0]),
        ]
    order = np.lexsort(uniqs.T)[::-1]  # decreasing ...
    best = uniqs[order[0]][1]  # number of unique spots
    last = int(best * uniqcut)
    grains = []
    for i in order:
        ntotal, nuniq = uniqs[i]
        if nuniq < last:
            break
        grains.append((ntotal, nuniq, ind.ubis[i]))
    return grains


def proxy(args):
    """Wrapper function for multiprocessing"""
    i, j, memcolf, idxopts = args
    sm = {}
    for col in "gx", "gy", "gz", "dtyi", "isel", "cosomega", "sinomega", "eta":
        shape, dtype, mem = memcolf[col]
        sm[col] = np.ndarray(shape, dtype, buffer=mem.buf)
    g = idxpoint(
        i,
        j,
        sm["isel"],
        sm["cosomega"],
        sm["sinomega"],
        sm["dtyi"],
        sm["gx"],
        sm["gy"],
        sm["gz"],
        sm["eta"],
        **idxopts
    )
    return i, j, g


# one of these per process:
ucglobal = None
symglobal = None
parglobal = None


def initializer(parfile, symmetry, loglevel=3):
    global ucglobal, symglobal, parglobal
    parglobal = parameters.read_par_file(parfile)
    ucglobal = unitcell.unitcell_from_parameters(parglobal)
    symglobal = sym_u.getgroup(symmetry)()
    ImageD11.indexing.loglevel = loglevel


class PBP:
    def __init__(
        self,
        parfile,
        dset,
        hkl_tol=0.01,
        fpks=0.7,
        ds_tol=0.005,
        etacut=0.1,
        ifrac=None,
        gmax=5,
        cosine_tol=np.cos(np.radians(90 - 0.1)),
        y0=0,
        symmetry="cubic",
        loglevel=3,
        foridx=None,
        forgen=None,
        uniqcut=0.75,
    ):
        """
        parfile = ImageD11 parameter file (for the unit cell + geometry)
        dsname = name of dset file, or object with "ybincens" array

        """
        self.parfile = parfile
        self.hkl_tol = hkl_tol
        self.fpks = fpks
        self.ds_tol = ds_tol
        self.etacut = etacut
        self.symmetry = symmetry
        self.foridx = foridx
        self.forgen = forgen
        self.gmax = gmax
        self.dset = dset
        self.ybincens = np.array(dset.ybincens)
        if ifrac is None:
            self.ifrac = 1.0 / len(self.ybincens)
        else:
            self.ifrac = ifrac
        self.ystep = abs(self.ybincens[-1] - self.ybincens[0]) / (
            len(self.ybincens) - 1
        )
        self.y0 = y0
        self.uniqcut = uniqcut
        self.cosine_tol = cosine_tol
        self.loglevel = loglevel

    def setpeaks(self, colf):
        """
        This is called on the main parent process
        The workers only need the init and they will pick up the shared
        memory defined in here
        """
        # Load the peaks
        colf.parameters.loadparameters(self.parfile)
        colf.updateGeometry()  # for ds
        #
        uc = unitcell.unitcell_from_parameters(colf.parameters)
        uc.makerings(colf.ds.max())
        # peaks that are on rings
        sel = np.zeros(colf.nrows, bool)
        # rings to use for indexing
        isel = np.zeros(colf.nrows, bool)
        npks = 0
        if self.foridx is None:
            self.foridx = range(len(uc.ringds))
        hmax = 0
        for i in range(len(uc.ringds)):
            ds = uc.ringds[i]
            rm = abs(colf.ds - ds) < self.ds_tol
            hkls = uc.ringhkls[ds]
            if i in self.foridx and rm.sum() > 0:
                icut = np.max(colf.sum_intensity[rm]) * self.ifrac
                rm = rm & (colf.sum_intensity > icut)
                isel |= rm
                npks += len(hkls)
                print(
                    i,
                    "%.4f" % (ds),
                    hkls[-1],
                    len(hkls),
                    rm.sum(),
                    "used, sum_intensity>",
                    icut,
                )
                sel |= rm
            else:
                print(i, "%.4f" % (ds), hkls[-1], len(hkls), rm.sum(), "skip")
            for hkl in hkls:
                hmax = max(np.abs(hkl).max(), hmax)
        if self.forgen is None:
            self.forgen = self.foridx

        isel = isel & ((np.abs(np.sin(np.radians(colf.eta)))) > self.etacut)

        colf.addcolumn(isel, "isel")  # peaks selected for indexing
        colf.addcolumn(np.sin(np.radians(colf.omega)), "sinomega")
        colf.addcolumn(np.cos(np.radians(colf.omega)), "cosomega")
        colf.addcolumn(np.round((colf.dty - self.y0) / self.ystep).astype(int), "dtyi")
        # peaks that are on any rings
        self.colf = colf
        self.icolf = colf.copy()
        self.icolf.filter(sel)
        self.npks = npks
        if self.fpks < 1:
            self.minpks = int(npks * self.fpks)
        else:
            self.minpks = self.fpks
        print(
            "Using for indexing:",
            self.icolf.nrows,
            "npks, minpks, forgen",
            self.npks,
            self.minpks,
            self.forgen,
        )
        self.hmax = hmax

    def iplot(self):
        import pylab as pl

        f, ax = pl.subplots(2, 1)
        ax[0].plot(self.colf.tth, self.colf.sum_intensity, ",")
        ax[0].plot(self.icolf.tth, self.icolf.sum_intensity, ",")
        ax[0].set(yscale="log", ylabel="sum intensity", xlabel="tth")
        histo = self.dset.sinohist(
            omega=self.icolf.omega, dty=self.icolf.dty, weights=self.icolf.sum_intensity
        )
        f.colorbar(
            ax[1].pcolormesh(
                self.dset.obinedges,
                self.dset.ybinedges,
                histo.T,
                norm=pl.matplotlib.colors.LogNorm(),
            ),
            ax=ax[1],
        )
        ax[1].set(ylabel="dty", xlabel="omega")
        return f, ax

    def point_by_point(
        self,
        grains_filename,
        nprocs=None,
        gridstep=1,
        debugpoints=None,
        loglevel=None,
    ):
        if loglevel is not None:
            self.loglevel = loglevel
        start = time.time()
        # FIXME - look at dataset ybincens and roi_iradon output_size ?
        if nprocs is None:
            nprocs = len(os.sched_getaffinity(os.getpid()))

        idxopt = {
            "ystep": self.ystep,
            "y0": self.y0,
            "minpks": self.minpks,
            "hkl_tol": self.hkl_tol,
            "ds_tol": self.ds_tol,
            "cosine_tol": self.cosine_tol,
            "forgen": self.forgen,
            "hmax": self.hmax,
            "uniqcut": self.uniqcut,
        }
        pprint.pprint(idxopt)

        nlo = np.floor((self.ybincens.min() - self.y0) / self.ystep).astype(int)
        nhi = np.ceil((self.ybincens.max() - self.y0) / self.ystep).astype(int)
        rng = range(nlo, nhi + 1, gridstep)

        with SharedMemoryManager() as smm:  # this runs in the parent process
            memcolf = {}
            for col in "gx", "gy", "gz", "dtyi", "isel", "cosomega", "sinomega", "eta":
                ar = getattr(self.icolf, col)
                mem = smm.SharedMemory(size=ar.itemsize * ar.size)
                mar = np.ndarray(ar.shape, dtype=ar.dtype, buffer=mem.buf)
                mar[:] = ar[:]  # copy into shm
                memcolf[col] = (
                    ar.shape,
                    ar.dtype,
                    mem,
                )  # this is a shared memory reference

            if debugpoints is None:
                args = [(i, j, memcolf, idxopt) for i in rng for j in rng]
                random.shuffle(args)
            else:
                args = [(i, j, memcolf, idxopt) for i, j in debugpoints]
            gmap = {}
            t0 = time.time()
            with open(grains_filename, "w") as gout:
                gout.write(
                    "#  i  j  ntotal  nuniq  ubi00  ubi01  ubi02  ubi10  ubi11  ubi12  ubi20  ubi21  ubi22\n"
                )
                done = 0
                ng = 0
                with multiprocessing.Pool(
                    nprocs,
                    initializer=initializer,
                    initargs=(self.parfile, self.symmetry, self.loglevel),
                ) as p:
                    for i, j, g in p.imap_unordered(proxy, args):
                        gmap[i, j] = g
                        done += 1
                        ng += len(g)
                        if done % nprocs == 0:
                            sys.stdout.flush()  # before!
                            dt = time.time() - t0
                            print(
                                "Done %7.3f %%, average grains/point %6.2f, %.3f /s/point, total %.1f /s"
                                % (
                                    100.0 * done / len(args),
                                    float(ng) / done,
                                    dt / done,
                                    dt,
                                ),
                                end="\r",
                            )
                        sys.stdout.flush()
                        for k in range(min(self.gmax, len(g))):
                            # only output gmax grains per point max
                            n, u, ubi = g[k]
                            # print(g[k])
                            gout.write("%d  %d  %d  %d  " % (i, j, n, u))
                            gout.write(("%f " * 9) % tuple(ubi.ravel()) + "\n")
                            gout.flush()

        end = time.time()
        print(end - start, "seconds", (end - start) / len(args), "s per point")
        print(end - t0, "seconds", (end - t0) / len(args), "s per point without setup")
