
#include "cImageD11.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * neighbormax assigns values in "l" to address highest pixel
 * of the 8 neighbors.
 *  l == 0 => this pixel is the max
 *  otherwise o[ l-1 ] give the offset
 */

#define pick(A, B, I, J)                                                       \
    if ((A) > (B)) {                                                           \
        (B) = (A);                                                             \
        (I) = (J);                                                             \
    }

DLL_LOCAL
int neighbormax(const float *restrict im, // input
                int32_t *restrict lout,   // output
                uint8_t *restrict l,      // workspace temporary
                int dim0,                 // Image dimensions
                int dim1) {
    // Fill in the first label as npeaks on this row
    int i, j, p, k0, k1, k2, npks;
    float mx0, mx1, mx2;

    npks = 0;
    /* The image borders */
    for (i = 0; i < dim1; i++) { // first and last row here:
        lout[i] = 0;             // ends of rows are below
        l[i] = 0;
        lout[dim1 * (dim0 - 1) + i] = 0;
        l[dim1 * (dim0 - 1) + i] = 0;
    }
#pragma omp parallel for private(j, p, k0, k1, k2, mx0, mx1, mx2)              \
    reduction(+ : npks)
    for (i = dim1; i < (dim0 - 1) * dim1; i = i + dim1) {
        // skipping 1 pixel border
        lout[i] = 0; // set edges to zero: pixel j=0:
        l[i] = 0;

        p = i + 1;
        mx0 = im[p - 1 - dim1];
        k0 = 1;
        pick(im[p - 1], mx0, k0, 2);
        pick(im[p - 1 + dim1], mx0, k0, 3);
        k1 = 4;
        mx1 = im[p - dim1];
        pick(im[p], mx1, k1, 5);
        pick(im[p + dim1], mx1, k1, 6);
        for (j = 1; j < dim1 - 1; j++) {
            p = i + j;
            mx2 = im[p + 1 - dim1];
            k2 = 7;
            pick(im[p + 1], mx2, k2, 8);
            pick(im[p + 1 + dim1], mx2, k2, 9);
            pick(mx1, mx0, k0, k1);
            pick(mx2, mx0, k0, k2);
            l[p] = k0;
            if (k0 == 5) { /* This peak is a maximum, count in lout[i] */
                lout[i]++;
            }
            mx0 = mx1;
            k0 = k1 - 3;
            mx1 = mx2;
            k1 = k2 - 3;
        }
        // set edges to zero: pixel j=dim1:
        lout[i + dim1 - 1] = 0;
        l[i + dim1 - 1] = 0;
        npks += lout[i];
    } // i
    return npks;
}

/* F2PY_WRAPPER_START
    function localmaxlabel( data, labels, wrk, ns, nf )
        intent(c) localmaxlabel
!DOC localmaxlabel assigns a label for each pixel so they are grouped
!DOC to the local maximum. Equal values choose to assign towards the earlier
!DOC value in memory.
!DOC cpu arg (1)0=C, (1)1=SSE2, (1)2=AVX2; if > 9 prints timing
        intent(c)
        real, intent(in) :: data(ns,nf)
        integer*4, intent(inout), note(hello) :: labels(ns,nf)
        integer*1, intent(inout) :: wrk(ns,nf)
        integer, intent(hide), depend(data) :: ns=shape(data,0)
        integer, intent(hide), depend(data) :: nf=shape(data,1)
        ! Returns
        integer :: localmaxlabel
    end function localmaxlabel
F2PY_WRAPPER_END */
int localmaxlabel(const float *restrict im, // input
                  int32_t *restrict lout,   // output
                  uint8_t *restrict l,      // workspace temporary
                  int dim0,                 // Image dimensions
                  int dim1) {
    // old msvc for python 2.7 requires ALL variables declared up here.
    int i, j, p, k, q, npk, t, nt, tid, lo, hi, npks;
    //   int noisy=0;
#define noisy 0
    double tic, toc;
    int o[10] = {0, // special case
                 -1 - dim1, -1,        -1 + dim1, -dim1,    0,
                 +dim1,     +1 - dim1, +1,        +1 + dim1}; // 7,8,9

    if (noisy)
        printf("Not using intrinsics\n");
    npks = neighbormax(im, lout, l, dim0, dim1); //, o);
    if (noisy) {
        toc = my_get_time();
        printf("    neighbormax %.3f ms %d peaks\n", 1000 * (toc - tic), npks);
        tic = toc;
    }
    // Cumulate lout[i] so that lout[i] holds cumulative sums
    npk = 0;
    for (i = 0; i < dim0 * dim1; i = i + dim1) {
        t = npk;
        npk += lout[i];
        lout[i] = t;
    }
    if (noisy) {
        toc = my_get_time();
        printf("    cumsum %.3f ms %d\n", 1000 * (toc - tic), npk);
        tic = toc;
    }
    // Now pass with row offsets in place
#pragma omp parallel for private(nt, t, j, p) schedule(dynamic)
    for (i = 0; i < (dim0 - 1); i++) {
        t = lout[i * dim1];
        nt = lout[(i + 1) * dim1];
        if (t == nt)
            continue;
        // Break early if there is no peak on this row or you already got it
        for (j = 1; j < (dim1 - 1); j++) {
            p = dim1 * i + j;
            if (l[p] == 5) { // pointing to self : k+1==5
                t++;
                lout[p] = t; // final label
                l[p] = 0;    // done, so tagged as zero
                if (t == nt)
                    break; // no more peaks
            }
        }
    }
    /* overwrite front border as zero again */
    for (i = 0; i < dim0 * dim1; i = i + dim1) {
        lout[i] = 0;
    }
    if (noisy) {
        toc = my_get_time();
        printf("    relabel %.3f ms\n", 1000 * (toc - tic));
        tic = toc;
    }
    //
    // Now make all point to their max
    // If we re-write the paths we cannot run in parallel
    //  ... this was the slowest part, so try openmp
    // 105 ms to always walk to max
    // 40 ms if path relabelled
    //
    // This is the same for all versions (optimised or not)
    //  ... perhaps re-write to be a manual loop and fill in
    //  ... the steps that are thread local
#pragma omp parallel private(q, i, tid, nt, k, lo, hi)
    {
#ifdef _OPENMP
        tid = omp_get_thread_num();
        nt = omp_get_num_threads();
#else
        tid = 0;
        nt = 1;
#endif
        lo = dim0 * dim1 * tid / nt;
        hi = dim0 * dim1 * (tid + 1) / nt;
        for (i = lo; i < hi; i++) {
            if (l[i] == 0)
                continue; // done
            // Now we need to walk to find a label
            k = 0;
            q = i + o[l[i]];
            while (l[q]) {
                q = q + o[l[q]];
                k++;
            } // Now q addresses a max or a correct label
            lout[i] = lout[q]; // take label from max
            if (k >
                0) { // relabel the path taken while we know the top label value
                q = i + o[l[i]];
                while (l[q]) {
                    if ((q >= lo) && (q < hi)) {
                        l[q] = 0;
                        lout[q] = lout[i];
                    }
                    q = q + o[l[q]];
                }
            }
            l[i] = 0; // sharing problems??
        }
    }
    if (noisy) {
        toc = my_get_time();
        printf("    write %.3f ms\n", 1000 * (toc - tic));
        tic = toc;
    }
    return npk;
}
