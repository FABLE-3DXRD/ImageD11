
#include "blobs.h" /* Disjoint sets thing for blob finding */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

DLL_LOCAL
void boundscheck(int jpk, int n2, int ipk, int n1) {
    if ((jpk < 0) || jpk >= n2) {
        printf("Bounds check error, jpk, n2\n");
        exit(0);
    }
    if ((ipk < 0) || ipk >= n1) {
        printf("Bounds check error, jpk, n1\n");
        exit(0);
    }
}

/*
# ImageD11_v1.x Software for beamline ID11
# Copyright (C) 2005-2017  Jon Wright
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/

/* ==== connectedpixels ======================================== */

/* F2PY_WRAPPER_START
    function connectedpixels( data, labels, threshold, &
                              verbose, con8, ns, nf)
        intent(c) connectedpixels
!DOC connectedpixels Determines which pixels in data are above the
!DOC user supplied threshold and assigns them into connected objects
!DOC which are output in labels. Connectivity is 3x3 box (8) by default
!DOC and reduces to a +(4) is con8==0
        intent(c)
        real, intent(in) :: data(ns,nf)
        integer, intent(inout) :: labels(ns,nf)
        integer, intent(hide), depend(data) :: ns=shape(data,0)
        integer, intent(hide), depend(data) :: nf=shape(data,1)
        integer, optional :: con8 = 1
        integer, optional ::  verbose = 0
        real threshold
        ! Returns
        integer :: connectedpixels
    end function connectedpixels
F2PY_WRAPPER_END */

// Exported function is in the C wrapper, not here!
DLL_LOCAL
int connectedpixels(float *data, int32_t *labels, float threshold, int verbose,
                    int eightconnected, int ns, int nf) {

    int i, j, irp, ir, ipx;
    int32_t k, *S, *T, np;

    if (verbose) {
        printf("Welcome to connectedpixels ");
        if (eightconnected)
            printf("Using connectivity 8\n");
        else
            printf("Using connectivity 4\n");
    }

    /* lots of peaks possible */
    S = dset_initialise(16384);

    /* To simplify later we hoist the first row and first pixel
     * out of the loops
     *
     * Algorithm scans image looking at stuff previously seen
     */

    /* First point */
    /*  i = 0;   j = 0; */
    if (data[0] > threshold) {
        S = dset_new(&S, &labels[0]);
    } else {
        labels[0] = 0;
    }
    /* First row */
    for (j = 1; j < nf; j++) {
        labels[j] = 0; /* initialize */
        if (data[j] > threshold) {
            if (labels[j - 1] > 0) {
                labels[j] = labels[j - 1];
            } else {
                S = dset_new(&S, &labels[j]);
            }
        }
    }

    /* === Mainloop ============================================= */
    for (i = 1; i < ns; i++) { /* i-1 prev row always exists, see above */
        ir = i * nf;           /* this row */
        irp = ir - nf;         /* prev row */
        /* First point */
        /* j=0; */
        labels[ir] = 0;
        if (data[ir] > threshold) {
            if (labels[irp] > 0) {
                labels[ir] = labels[irp];
            }
            if (eightconnected && (labels[irp + 1] > 0)) {
                match(labels[ir], labels[irp + 1], S);
            }
            if (labels[ir] == 0) {
                S = dset_new(&S, &labels[ir]);
            }
        }
        /* Run along row to just before end */
        for (j = 1; j < nf - 1; j++) {
            ipx = ir + j;
            irp = ipx - nf;
            labels[ipx] = 0;
            if (data[ipx] > threshold) {
                /* Pixel needs to be assigned */
                if (eightconnected && (labels[irp - 1] > 0)) {
                    match(labels[ipx], labels[irp - 1], S);
                }
                if (labels[irp] > 0) {
                    match(labels[ipx], labels[irp], S);
                }
                if (eightconnected && (labels[irp + 1] > 0)) {
                    match(labels[ipx], labels[irp + 1], S);
                }
                if (labels[ipx - 1] > 0) {
                    match(labels[ipx], labels[ipx - 1], S);
                }
                if (labels[ipx] == 0) { /* Label is new ! */
                    S = dset_new(&S, &labels[ipx]);
                }
            } /* (val > threshold) */
        } /* Mainloop j */
        /* Last pixel on the row */
        ipx = ir + nf - 1;
        irp = ipx - nf;
        labels[ipx] = 0;
        if (data[ipx] > threshold) {
            if (eightconnected && (labels[irp - 1] > 0)) {
                match(labels[ipx], labels[irp - 1], S);
            }
            if (labels[irp] > 0) {
                match(labels[ipx], labels[irp], S);
            }
            if (labels[ipx - 1] > 0) {
                match(labels[ipx], labels[ipx - 1], S);
            }
            if (labels[ipx] == 0) { /* Label is new ! */
                S = dset_new(&S, &labels[ipx]);
            }
        }
    }
    /* Now compress the disjoint set to make single list of
     * unique labels going from 1->n
     */
    T = dset_compress(&S, &np);
    /* Now scan through image re-assigning labels as needed */
#pragma omp parallel for private(j, ipx, k) shared(labels)
    for (i = 0; i < ns; i++) {
        for (j = 0; j < nf; j++) {
            ipx = i * nf + j;
            k = labels[ipx];
            if (k > 0) {
                if (T[k] == 0) {
                    printf("Error in connectedpixels\n");
                }
                if (T[k] != k) {
                    labels[i * nf + j] = T[k];
                }
            }
        }
    }
    free(S);
    free(T);
    return np;
}

/* F2PY_WRAPPER_START
    subroutine blobproperties( data, labels, np, omega, &
                               verbose, ns, nf, results)
        intent(c) blobproperties
!DOC blobproperties fills the array results with properties of each labelled
!DOC object described by data (pixel values) and labels. The omega value
!DOC is the angle for this frame.
!DOC results are FIXME
        intent(c)
        real, intent(in) :: data(ns, nf)
        integer, intent(in) :: labels(ns, nf)
        integer, intent(hide), depend(data) :: ns=shape(data,0)
        integer, intent(hide), depend(data) :: nf=shape(data,1)
        integer, intent(in) :: np
        double precision, intent(out) :: results( np, NPROPERTY )
        real, intent(in), optional :: omega = 0
        integer, optional :: verbose = 0
        threadsafe
    end subroutine blobproperties
F2PY_WRAPPER_END */
void blobproperties(float *data, int32_t *labels, int32_t npk, float omega,
                    int verbose, int ns, int nf, double *res) {
    int i, j, bad, ipx;
    double fval;
    int32_t ipk;
    if (verbose) {
        printf("Computing blob moments, ns %d, nf %d, npk %d\n", ns, nf, npk);
    }
    /* Initialise the results */
    for (i = 0; i < npk; i++) {
        for (j = 0; j < NPROPERTY; j++) {
            res[i * NPROPERTY + j] = 0.;
        }
        /* Set min to max +1 and vice versa */
        res[i * NPROPERTY + bb_mn_f] = nf + 1;
        res[i * NPROPERTY + bb_mn_s] = ns + 1;
        res[i * NPROPERTY + bb_mx_f] = -1;
        res[i * NPROPERTY + bb_mx_s] = -1;
        /* All pixels have the same omega in this frame */
        res[i * NPROPERTY + bb_mx_o] = omega;
        res[i * NPROPERTY + bb_mn_o] = omega;
    }
    if (verbose != 0)
        printf("Scanning image\n");

    bad = 0;
    /* i,j is looping along the indices data array */
    for (i = 0; i < ns; i++) {
        for (j = 0; j < nf; j++) {
            ipx = i * nf + j;
            ipk = labels[ipx];
            if (ipk > 0 && ipk <= npk) {
                fval = (double)data[ipx];
                add_pixel(&res[NPROPERTY * (ipk - 1)], i, j, fval, omega);
            } else {
                if (ipk != 0) {
                    bad++;
                    if (bad < 10) {
                        printf("Found %d in your blob image at i=%d, j=%d\n",
                               ipk, i, j);
                    }
                }
            }
        } /* j */
    } /* i */
    if (verbose) {
        printf("\nFound %d bad pixels in the blob image\n", bad);
    }
}

/* F2PY_WRAPPER_START
    function bloboverlaps( labels1, npk1, results1,    &
                           labels2, npk2, results2,    &
                           verbose, ns, nf)
!DOC bloboverlaps determines the overlaps between labels1 and labels2
!DOC for an image series. Peaks in labels2 may be merged if they were
!DOC joined by a peak on labels1. Results in results1 are accumulated
!DOC into results2 if peaks are overlapped.
        intent(c) bloboverlaps
        intent(c)
        integer :: bloboverlaps
        integer, intent( inout ) :: labels1( ns, nf )
        integer, intent( inout ) :: labels2( ns, nf )
        integer, intent(hide), depend(labels1) :: ns=shape(labels1,0)
        integer, intent(hide), depend(labels1) :: nf=shape(labels1,1)
        integer, intent(in) :: npk1, npk2
        double precision, intent( inout ) :: results1( :, NPROPERTY )
        double precision, intent( inout ) :: results2( :, NPROPERTY )
        integer, intent(in) :: verbose = 0
        threadsafe
    end subroutine bloboverlaps
F2PY_WRAPPER_END */
int bloboverlaps(int32_t *b1, int32_t n1, double *res1, int32_t *b2, int32_t n2,
                 double *res2, int verbose, int ns, int nf) {

    int i, j, safelyneed, ipx;
    int32_t *link, p1, p2, ipk, jpk, npk, *T;

    /* Initialise a disjoint set in link
     * image 2 has peak[i]=i ; i=1->n2
     * image 1 has peak[i]=i+n1+1 ; i=1->n1
     *                          0, 1, 2, 3, n2
     *                          4, 5, 6, 7, 8, 9, n1   need 0, 4 == n1+n2+3
     * link to hold 0->n1-1 ; n1->n2+n2-1
     * */

    /* This is a disjoint merge operation ... */
    if (verbose)
        printf("Enter bloboverlaps\n");
    safelyneed = n1 + n2 + 3;
    link = (int *)malloc(safelyneed * sizeof(int32_t));
    link[0] = safelyneed;
    for (i = 1; i < safelyneed; i++) {
        link[i] = i;
    }
    /* flag the start of image number 2 */
    link[n2 + 1] = -99999; /* ==n2=0 Should never be touched by anyone */
    /* results lists of pairs of number */
    /* link holds a disjoint set, we label directly overlapping pixels (i==j) */
    for (i = 0; i < ns; i++) {
        for (j = 0; j < nf; j++) {
            ipx = i * nf + j;
            if ((p1 = b1[ipx]) == 0)
                continue;
            if ((p2 = b2[ipx]) == 0)
                continue;
            if (link[p2] < 0 || link[p1 + n2 + 1] < 0) {
                printf("Whoops!!\n");
                return 0;
            }
            dset_makeunion(link, p2, p1 + n2 + 1);
            if (verbose > 10)
                printf("link %d %d\n", p2, p1 + n2 + 1);
        }
    }
    if (verbose)
        printf("Scanning images\n");
    /* Now we re-label and merge peaks scanning disjoint set */
    for (i = 1; i < safelyneed; i++) {
        if (link[i] != i && i != n2 + 1) {
            j = dset_find(i, link);
            if ((i > n2 + 1) && (j < n2 + 1)) { /* linking between images */
                jpk = j - 1;
                ipk = i - n2 - 2;
                /* Bounds checking */
                boundscheck(jpk, n2, ipk, n1);
                merge(&res2[NPROPERTY * jpk], &res1[NPROPERTY * ipk]);
                if (verbose > 10)
                    printf("merged res2[%d] res1[%d]\n", jpk, ipk);
                continue;
            }
            if ((i > n2 + 1) &&
                (j > n2 + 1)) { /* linking on the same image (1) */
                jpk = j - n2 - 2;
                ipk = i - n2 - 2;
                boundscheck(jpk, n1, ipk, n1);
                assert((n1 > jpk) && (jpk >= 0) && (n1 > ipk) && (ipk >= 0));
                merge(&res1[NPROPERTY * jpk], &res1[NPROPERTY * ipk]);
                if (verbose > 10)
                    printf("merge res1[%d] res1[%d]\n", jpk, ipk);
                continue;
            }
            if (i < n2 + 1 && j < n2 + 1) { /* linking on the same image (2) */
                jpk = j - 1;
                ipk = i - 1;
                boundscheck(jpk, n2, ipk, n2);
                merge(&res2[NPROPERTY * jpk], &res2[NPROPERTY * ipk]);
                if (verbose > 10)
                    printf("merge res2[%d] res2[%d]\n", jpk, ipk);
                continue;
            }
            assert("I am not here!");
        }
    }

    /* This is the case where two spots on the current image become linked by
     * by a spot overlap on the previous one
     *
     * The labels are now wrong, in fact there is a single peak in 3D and
     * two peaks in 2D
     *
     * Thanks to Stine West from Riso for finding this subtle bug
     */

    /* First, work out the new labels */
    /* Make each T[i] contain the unique ascending integer for the set */
    if (verbose)
        printf("Compress set\n");
    T = (int32_t *)(calloc((n2 + 3), sizeof(int32_t)));
    assert(T != NULL);

    npk = 0;
    for (i = 1; i < n2 + 1; i++) {
        if (link[i] == i) {
            npk = npk + 1;
            T[i] = npk;
        } else {
            j = dset_find(i, link);
            assert(j < i);
            T[i] = T[j];
        }
    }
    if (verbose) {
        printf("n1 = %d n2 = %d ", n1, n2);
        for (i = 0; i < n2 + 3; i++)
            printf("T[%d]=%d ", i, T[i]);
    }
    /* T is now compressed, merge the peaks */
    if (verbose)
        printf("Merge peaks in res2\n");
    for (i = 1; i < n2 + 1; i++) {
        /* dest  = T[i]; */
        /* src   = link[i]; */
        if (link[i] == T[i])
            continue;
        if (T[i] < link[i]) {   /* copy and zero out */
            if (link[i] == i) { /* This is the place accumulating */
                for (j = 0; j < NPROPERTY; j++) {
                    res2[NPROPERTY * (T[i] - 1) + j] =
                        res2[NPROPERTY * (link[i] - 1) + j];
                    res2[NPROPERTY * (link[i] - 1) + j] = 0;
                }
            } else {
                /* assert this is empty */
                assert(res2[NPROPERTY * (i - 1) + s_1] < 0.01);
            }
            if (verbose > 1) {
                printf("np i %d j %d %f \n", T[i], link[i],
                       res2[NPROPERTY * T[i] + 1]);
            }
        } else {
            assert("Bad logic in bloboverlaps");
        }
    }
    if (verbose)
        printf("Relabel image of blobs\n");
    /* Relabel the image where they change */
    for (i = 0; i < ns; i++) {
        for (j = 0; j < nf; j++) {
            ipx = i * nf + j;
            if ((p2 = b2[ipx]) == 0)
                continue;
            ipk = T[p2];
            if (ipk != p2) {
                assert(ipk > 0 && ipk <= n2);
                b2[ipx] = ipk;
            }
        }
    }
    if (verbose)
        printf("Done relabelling\n");
    free(T);
    free(link);
    return npk;
}

/* F2PY_WRAPPER_START
    subroutine blob_moments( results, np )
!DOC blob_moments fills in the reduced moments in results array.
!DOC ... FIXME - this would be clearer in python, fast anyway.
        intent(c) blob_moments
        intent(c)
        double precision, intent( inout ) :: results( np, NPROPERTY )
        integer, intent(hide), depend(results) :: np=shape(results,0)
        threadsafe
    end subroutine blob_moments
F2PY_WRAPPER_END */
void blob_moments(double *res, int np) { compute_moments(res, np); }

/* F2PY_WRAPPER_START
    function clean_mask( msk, ret, ns, nf )
!DOC clean_mask removes pixels which are not 4 connected from msk
!DOC while copying into ret.
        intent(c) clean_mask
        intent(c)
        integer*1, intent(in)  :: msk( ns, nf )
        integer, intent(hide), depend(msk) :: ns=shape(msk,0)
        integer, intent(hide), depend(msk) :: nf=shape(msk,1)
        integer*1, intent(inout), dimension(ns, nf) :: ret
        ! returns an int
        integer :: clean_mask
    end function clean_mask
F2PY_WRAPPER_END */
int clean_mask(const int8_t *restrict msk, int8_t *restrict ret, int ns,
               int nf) {
    /* cleans pixels with no 4 connected neighbors */
    int i, j, q, npx;
    int8_t t;
    npx = 0;
#pragma omp parallel for private(i)
    for (i = 0; i < ns * nf; i++) {
        if (msk[i] > 0) {
            ret[i] = 1;
        } else {
            ret[i] = 0;
        }
    }
    i = 0;
    for (j = 0; j < nf; j++) {
        q = i * nf + j;
        if (ret[q] > 0) {
            t = msk[q + nf];
            if (j > 0)
                t += msk[q - 1];
            if (j < (nf - 1))
                t += msk[q + 1];
            if (t > 0) {
                npx++;
            } else {
                ret[q] = 0;
            }
        }
    }
#pragma omp parallel for private(i, j, q, t) reduction(+ : npx)
    for (i = 1; i < (ns - 1); i++) {
        /* j==0 */
        q = i * nf;
        if (ret[q] > 0) {
            t = msk[q - nf] + msk[q + nf] + msk[q + 1];
            if (t > 0) {
                npx++;
            } else {
                ret[q] = 0;
            }
        }
        for (j = 1; j < (nf - 1); j++) {
            q = i * nf + j;
            if (ret[q] > 0) {
                t = msk[q - nf] + msk[q + nf] + msk[q - 1] + msk[q + 1];
                if (t > 0) {
                    npx++;
                } else {
                    ret[q] = 0;
                }
            }
        }
        q = (i + 1) * nf - 1;
        if (ret[q] > 0) {
            t = msk[q - nf] + msk[q + nf] + msk[q - 1];
            if (t > 0) {
                npx++;
            } else {
                ret[q] = 0;
            }
        }
    }
    i = ns - 1;
    for (j = 0; j < nf; j++) {
        q = i * nf + j;
        if (ret[q] > 0) {
            t = msk[q - nf];
            if (j > 0)
                t += msk[q - 1];
            if (j < (nf - 1))
                t += msk[q + 1];
            if (t > 0) {
                npx++;
            } else {
                ret[q] = 0;
            }
        }
    }
    return npx;
}

/* F2PY_WRAPPER_START

    function make_clean_mask( img, cut, msk, ret, ns, nf )
!DOC make_clean_mask is a lot like clean msk but it generates
!DOC the msk using img and cut.
!DOC Beware: work in progress
        intent(c) make_clean_mask
        intent(c)
        real, intent(in), dimension(ns,nf) :: img
        real, intent(in) :: cut
        integer*1, intent(in)  :: msk( ns, nf )
        integer, intent(hide), depend(msk) :: ns=shape(msk,0)
        integer, intent(hide), depend(msk) :: nf=shape(msk,1)
        integer*1, intent(inout), dimension(ns, nf) :: ret
        ! returns an int
        integer :: make_clean_mask
    end function make_clean_mask
F2PY_WRAPPER_END */
int make_clean_mask(float *restrict img, float cut, int8_t *restrict msk,
                    int8_t *restrict ret, int ns, int nf) {
    /* cleans pixels with no 4 connected neighbors */
    int i;
#pragma omp parallel for
    for (i = 0; i < ns * nf; i++) {
        if (img[i] > cut) {
            msk[i] = 1;
        } else {
            msk[i] = 0;
        }
    }
    return clean_mask(&msk[0], &ret[0], ns, nf);
}
