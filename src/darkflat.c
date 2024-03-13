/* stdint */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <float.h>
#include <math.h>

#include "cImageD11.h"

/* To compare to numpy (numba ought to do this anyway)

__m256i _mm256_cvtepu16_epi32 (__m128i a)
__m256 _mm256_cvtepi32_ps (__m256i a)

SSE: __m128 _mm_cvtpu16_ps (__m64 a)


AVX
__m256 _mm256_sub_ps (__m256 a, __m256 b)
__m256 _mm256_mul_ps (__m256 a, __m256 b)
  Architecture	Latency	Throughput (CPI)
  Skylake	4	0.5
  [divide is much slower]


  _m256_load_ps (float const * mem_addr)
  _m256_store_ps (float * mem_addr, __m256 a)
  aligned load/store (32 byte boundary)

  _mm256_loadu_ps (float const * mem_addr)
  _mm256_storeu_ps (float const * mem_addr)
  unaligned load/store (32 byte boundary)
*/

/* F2PY_WRAPPER_START
    subroutine uint16_to_float_darksub( img, drk, data, npx )
!DOC uint16_to_float_darksub subtracts image drk(float32) from
!DOC raw data in data (uint16) and returns in img.
        intent(c) uint16_to_float_darksub
        intent(c)
        real, intent(inout), dimension(npx) :: img
        real, intent(in), dimension(npx) :: drk
        integer(kind=-2), intent(in), dimension(npx) :: data
        integer, intent(hide), depend( img ) :: npx
    end subroutine uint16_to_float_darksub
F2PY_WRAPPER_END */
void uint16_to_float_darksub(float *restrict img, const float *restrict drk,
                             const uint16_t *restrict data, int npx) {
    int i;
#ifdef GOT_OMP_SIMD
#pragma omp parallel for simd
#else
#pragma omp parallel for
#endif
    for (i = 0; i < npx; i++) {
        img[i] = ((float)data[i]) - drk[i];
    }
}

/* F2PY_WRAPPER_START
    subroutine uint16_to_float_darkflm( img, drk, flm, data, npx )
!DOC uint16_to_float_darkflm subtracts image drk(float32) from
!DOC raw data in data (uint16), multiples by flm(float32) and returns in img.
        intent(c) uint16_to_float_darkflm
        intent(c)
        real, intent(inout), dimension(npx) :: img
        real, intent(in), dimension(npx) :: drk, flm
        integer(kind=-2), intent(in), dimension(npx) :: data
        integer, intent(hide), depend( img ) :: npx
    end subroutine uint16_to_float_darkflm
F2PY_WRAPPER_END */
void uint16_to_float_darkflm(float *restrict img, const float *restrict drk,
                             const float *restrict flm,
                             const uint16_t *restrict data, int npx) {
    int i;
#ifdef GOT_OMP_SIMD
#pragma omp parallel for simd
#else
#pragma omp parallel for
#endif
    for (i = 0; i < npx; i++) {
        img[i] = (((float)data[i]) - drk[i]) * flm[i];
    }
}

/* F2PY_WRAPPER_START
    subroutine frelon_lines(img, ns, nf, cut)
!DOC frelon_lines Subtracts the average value of (pixels < cut) per row
        intent(c) frelon_lines
        intent(c)
        real, intent(inout), dimension(ns,nf) :: img
        integer, intent(hide), depend(img) :: ns = shape(img,0)
        integer, intent(hide), depend(img) :: nf = shape(img,1)
        real, intent(in) :: cut
    end subroutine frelon_lines
F2PY_WRAPPER_END */
void frelon_lines(float *img, int ns, int nf, float cut) {
    int i, j, p, npx;
    float rowsum, avg;
    avg = img[0];
#pragma omp parallel for private(i, j, rowsum, npx, p) firstprivate(avg)
    for (i = 0; i < ns; i++) {
        rowsum = 0.;
        npx = 0;
        p = i * nf;
#ifdef GOT_OMP_SIMD
#pragma omp simd
#endif
        for (j = 0; j < nf; j++) {
            if (img[p + j] < cut) {
                rowsum += img[p + j];
                npx++;
            }
        }
        if (npx > 0)
            avg = rowsum / npx;
#ifdef GOT_OMP_SIMD
#pragma omp simd
#endif
        for (j = 0; j < nf; j++)
            img[p + j] = img[p + j] - avg;
    }
}

/* F2PY_WRAPPER_START
    subroutine frelon_lines_sub(img, drk, ns, nf, cut)
!DOC frelon_lines_sub Subtracts drk from img and then same as frelon_lines
        intent(c) frelon_lines_sub
        intent(c)
        real, intent(inout), dimension(ns,nf) :: img, drk
        integer, intent(hide), depend(img) :: ns = shape(img,0)
        integer, intent(hide), depend(img) :: nf = shape(img,1)
        real, intent(in) :: cut
        threadsafe
    end subroutine frelon_lines
F2PY_WRAPPER_END */
void frelon_lines_sub(float *restrict img, float *restrict drk, int ns, int nf,
                      float cut) {
    int i, j, p, npx;
    float rowsum, avg;
    avg = img[0];
#pragma omp parallel for private(i, j, rowsum, npx, p) firstprivate(avg)
    for (i = 0; i < ns; i++) {
        rowsum = 0.;
        npx = 0;
        p = i * nf;
#ifdef GOT_OMP_SIMD
#pragma omp simd
#endif
        for (j = 0; j < nf; j++) {
            img[p + j] = img[p + j] - drk[p + j];
            if (img[p + j] < cut) {
                rowsum += img[p + j];
                npx++;
            }
        }
        if (npx > 0)
            avg = rowsum / npx;
#ifdef GOT_OMP_SIMD
#pragma omp simd
#endif
        for (j = 0; j < nf; j++)
            img[p + j] = img[p + j] - avg;
    }
}

/* F2PY_WRAPPER_START
    subroutine array_mean_var_cut( img, npx, mean, var, n, cut, verbose )
!DOC array_mean_var_cut computes the mean and variance of an image
!DOC with pixels above the value mean+cut*stddev removed. This is iterated
!DOC n times as the mean and variance change as pixels are removed.
        intent(c) array_mean_var_cut
        real, intent(in,c), dimension(npx) :: img
        integer, intent(hide,c), depend(img) :: npx = shape(img,0)
        real, intent(out) :: mean, var
        integer, intent(in,c), optional :: n = 3
        integer, intent(in,c), optional :: verbose = 0
        real, intent(in,c), optional :: cut = 3.
    end subroutine array_mean_var_cut
F2PY_WRAPPER_END */
void array_mean_var_cut(float *restrict img, int npx, float *mean, float *std,
                        int n, float cut, int verbose) {
    int i, nactive;
    float t, s1, s2, wt, y0;
    y0 = img[0];
    s1 = 0;
    s2 = 0;
    if (verbose)
        printf("Args, img[0] %f npx %d n %d cut %f verbose %d\n", img[0], npx,
               n, cut, verbose);
#ifdef GOT_OMP_SIMD
#pragma omp parallel for simd private(t) reduction(+ : s1, s2)
#else
#pragma omp parallel for private(t) reduction(+ : s1, s2)
#endif
    for (i = 0; i < npx; i++) {
        t = img[i] - y0;
        s1 = s1 + t;
        s2 = s2 + t * t;
    }
    /* mean and std */
    *mean = (float)(s1 / npx + y0);
    *std = sqrtf((float)((s2 - (s1 * s1 / npx)) / npx));
    if (verbose > 0)
        printf("n=%d Mean %f, Std %f\n", n, *mean, *std);
    while (--n > 0) {
        y0 = *mean;
        wt = y0 + cut * (*std);
        s1 = 0;
        s2 = 0;
        nactive = 0;
#ifdef GOT_OMP_SIMD
#pragma omp parallel for simd private(t) reduction(+ : s1, s2, nactive)
#else
#pragma omp parallel for private(t) reduction(+ : s1, s2, nactive)
#endif
        for (i = 0; i < npx; i++) {
            if (img[i] < wt) {
                t = img[i] - y0;
                s1 = s1 + t;
                s2 = s2 + t * t;
                nactive++;
            }
        }
        *mean = (float)(s1 / nactive + *mean);
        *std = sqrtf(((s2 - (s1 * s1 / nactive)) / nactive));
        if (verbose > 0)
            printf("n=%d Mean %f, Std %f\n", n, *mean, *std);
    }
}

/* F2PY_WRAPPER_START
    subroutine array_mean_var_msk( img, msk, npx, mean, var, n, cut, verbose )
!DOC array_mean_var_msk computes the mean and variance of an image
!DOC with pixels above the value mean+cut*stddev removed. This is iterated
!DOC n times as the mean and variance change as pixels are removed.
        intent(c) array_mean_var_msk
        real, intent(in,c), dimension(npx) :: img
        integer(kind=-1) , intent(in,c), dimension(npx) :: msk
        integer, intent(hide,c), depend(img) :: npx = shape(img,0)
        real, intent(out) :: mean, var
        integer, intent(in,c), optional :: n = 3
        integer, intent(in,c), optional :: verbose = 0
        real, intent(in,c), optional :: cut = 3.
    end subroutine array_mean_var_msk
F2PY_WRAPPER_END */
void array_mean_var_msk(float *restrict img, uint8_t *restrict msk, int npx,
                        float *mean, float *std, int n, float cut,
                        int verbose) {
    int i, nactive;
    float t, s1, s2, wt, y0;
    y0 = img[0];
    s1 = 0;
    s2 = 0;
    if (verbose)
        printf("Args, img[0] %f npx %d n %d cut %f verbose %d\n", img[0], npx,
               n, cut, verbose);
#ifdef GOT_OMP_SIMD
#pragma omp parallel for simd private(t) reduction(+ : s1, s2)
#else
#pragma omp parallel for private(t) reduction(+ : s1, s2)
#endif
    for (i = 0; i < npx; i++) {
        t = img[i] - y0;
        s1 = s1 + t;
        s2 = s2 + t * t;
    }
    /* mean and std */
    *mean = (float)(s1 / npx + y0);
    *std = sqrtf((float)((s2 - (s1 * s1 / npx)) / npx));
    if (verbose > 0)
        printf("n=%d Mean %f, Std %f\n", n, *mean, *std);
    while (--n > 1) {
        y0 = *mean;
        wt = y0 + cut * (*std);
        s1 = 0;
        s2 = 0;
        nactive = 0;
#ifdef GOT_OMP_SIMD
#pragma omp parallel for simd private(t) reduction(+ : s1, s2, nactive)
#else
#pragma omp parallel for private(t) reduction(+ : s1, s2, nactive)
#endif
        for (i = 0; i < npx; i++) {
            if (img[i] < wt) {
                t = img[i] - y0;
                s1 = s1 + t;
                s2 = s2 + t * t;
                nactive++;
            }
        }
        *mean = (float)(s1 / nactive + *mean);
        *std = sqrtf(((s2 - (s1 * s1 / nactive)) / nactive));
        if (verbose > 0)
            printf("n=%d Mean %f, Std %f\n", n, *mean, *std);
    }

    /* Fill in mask */
    y0 = *mean;
    wt = y0 + cut * (*std);
    if (verbose > 0)
        printf("Cutting img > %f\n", wt);
#ifdef GOT_OMP_SIMD
#pragma omp parallel for simd
#else
#pragma omp parallel for
#endif
    for (i = 0; i < npx; i++) {
        if (img[i] < wt) {
            msk[i] = 0;
        } else {
            msk[i] = 1;
        }
    }
}

/* F2PY_WRAPPER_START

    subroutine array_stats( img, npx, minval, maxval, mean, var )
!DOC array_stats computes statistics for an image.
!DOC  img Input data array  (1D or 2D.ravel(), so sparse or dense)
!DOC npx length of data array (contiguous)
!DOC *minval minimum of the pixels
!DOC *maxval maximum of the pixels
!DOC*s1  Sum of all pixel
!DOC*s2  Sum of pixel^2
        intent(c) array_stats
        real, intent(c, in), dimension(npx) :: img
        integer, intent(c, hide), depend(img) :: npx = shape(img,0)
        ! these are intent(fortran) pointers
        real, intent(out) :: minval, maxval, mean, var
        threadsafe
    end subroutine array_stats
F2PY_WRAPPER_END */
void array_stats(float img[], int npx, float *minval, float *maxval,
                 float *mean, float *var) {
    int i;
    /* Use double to reduce rounding and subtraction errors */
    double t, s1, s2, y0, ts1, ts2;
    float mini, maxi, tmin, tmax;
    mini = FLT_MAX;
    maxi = FLT_MIN;
    s1 = 0.;
    s2 = 0.;
    y0 = img[0];
    /* Merge results - openmp 2.0 for windows has no min/max     */
#pragma omp parallel private(i, t, ts1, ts2, tmin, tmax)
    {
        tmin = FLT_MAX;
        tmax = FLT_MIN;
        ts1 = 0.;
        ts2 = 0.;
#ifdef GOT_OMP_SIMD
#pragma omp for simd
#else
#pragma omp for
#endif
        for (i = 0; i < npx; i++) {
            t = img[i] - y0;
            ts1 = ts1 + t;
            ts2 = ts2 + t * t;
            if (img[i] < tmin)
                tmin = img[i];
            if (img[i] > tmax)
                tmax = img[i];
        } // for
#pragma omp critical
        {
            s1 += ts1;
            s2 += ts2;
            if (tmin < mini)
                mini = tmin;
            if (tmax > maxi)
                maxi = tmax;
        }
    } // parallel
    /* results */
    *mean = (float)(s1 / npx + y0);
    *var = (float)((s2 - (s1 * s1 / npx)) / npx);
    *minval = mini;
    *maxval = maxi;
}

/* F2PY_WRAPPER_START
    subroutine array_histogram( img, npx, low, high, hist, nhist )
!DOC array_histogram computes the histogram for an image
!DOC Go through the data to compute a histogram of the values
!DOC  Previous call of array_stats would help to set up this call
!DOC   compare to np.bincount - this does not gain much
!DOC   better implementations can be done
!DOC
!DOC img[]  Input data array  (1D or 2D.ravel(), so sparse or dense)
!DOC npx    length of data array (contiguous)
!DOC low    Lower edge of first bin
!DOC high   Upper edge of last bin
!DOC hist[] Histogram to be output
!DOC nhist  Number of bins in the histogram
        intent(c) array_histogram
        intent(c)
        real, intent(in), dimension(npx) :: img
        integer, intent(hide), depend(img) :: npx = shape(img,0)
        real, intent(in) :: low, high
        integer*4, intent(inout), dimension(nhist) :: hist
        integer, intent(hide), depend(hist) :: nhist = shape(hist, 0 )
        threadsafe
    end subroutine array_histogram
F2PY_WRAPPER_END */
void array_histogram(float img[], int npx, float low, float high,
                     int32_t hist[], int nhist) {
    int i, ibin;
    float ostep;
    memset(hist, 0, nhist * sizeof(int32_t));
    /* Compute the multiplier to get the bin numbers */
    ostep = nhist / (high - low);
    for (i = 0; i < npx; i++) {
        ibin = (int)floorf((img[i] - low) * ostep);
        /* clip into range at ends */
        if (ibin < 0) {
            ibin = 0;
        }
        if (ibin >= nhist) {
            ibin = nhist - 1;
        }
        hist[ibin] = hist[ibin] + 1;
    }
}

/* F2PY_WRAPPER_START
    subroutine reorder_u16_a32(data, adr, out, N)
!DOC reorder_u16_a32 called in sandbox/fazit.py simple
!DOC loop with openmp saying out[adr[i]] in data[i]
!DOC e.g. semi-random writing
        intent(c) reorder_u16_a32
        intent(c)
        integer(kind=-2), dimension(N), intent(in) :: data
        integer*4, dimension(N), intent(in) :: adr
        integer(kind=-2), dimension(N), intent(inout) :: out
        integer, intent(hide), depend(data) :: N
    end subroutine reorder_u16_a32
F2PY_WRAPPER_END */
void reorder_u16_a32(const uint16_t *restrict data,
                     const uint32_t *restrict adr, uint16_t *restrict out,
                     int N) {
    int i;
    /*  printf("Hello, got N=%d\n",N);*/
#pragma omp parallel for
    for (i = 0; i < N; i++) {
        out[adr[i]] = data[i];
    }
}

/* F2PY_WRAPPER_START
    subroutine reorder_f32_a32(data, adr, out, N)
!DOC reorder_f32_a32 called in sandbox/fazit.py simple
!DOC loop with openmp saying out[adr[i]] in data[i]
!DOC e.g. semi-random writing
        intent(c) reorder_f32_a32
        intent(c)
        real, dimension(N), intent(in) :: data
        integer*4, dimension(N), intent(in) :: adr
        real, dimension(N), intent(inout) :: out
        integer, intent(hide), depend(data) :: N
    end subroutine reorder_u16_a32
F2PY_WRAPPER_END */
void reorder_f32_a32(const float *restrict data, const uint32_t *restrict adr,
                     float *restrict out, int N) {
    int i;
    /*  printf("Hello, got N=%d\n",N);*/
#pragma omp parallel for
    for (i = 0; i < N; i++) {
        out[adr[i]] = data[i];
    }
}

/* F2PY_WRAPPER_START
    subroutine reorderlut_u16_a32(data, adr, out, N)
!DOC reorderlut_u16_a32lut called in sandbox/fazit.py simple
!DOC loop with openmp saying out[i] in data[adr[i]]
!DOC e.g. semi-random reading
        intent(c) reorderlut_u16_a32
        intent(c)
        integer(kind=-2), dimension(N), intent(in) :: data
        integer*4, dimension(N), intent(in) :: adr
        integer(kind=-2), dimension(N), intent(inout) :: out
        integer, intent(hide), depend(data) :: N
    end subroutine reorderlut_u16_a32
F2PY_WRAPPER_END */
void reorderlut_u16_a32(uint16_t *restrict data, uint32_t *restrict lut,
                        uint16_t *restrict out, int N) {
    int i;
    /*  printf("Hello, got N=%d\n",N);*/
#pragma omp parallel for
    for (i = 0; i < N; i++) {
        out[i] = data[lut[i]];
    }
}

/* F2PY_WRAPPER_START
    subroutine reorderlut_f32_a32(data, adr, out, N)
!DOC reorderlut_f32_a32 lut called in sandbox/fazit.py simple
!DOC loop with openmp saying out[i] in data[adr[i]]
!DOC e.g. semi-random reading
        intent(c) reorderlut_f32_a32
        intent(c)
        real, dimension(N), intent(in) :: data
        integer*4, dimension(N), intent(in) :: adr
        real, dimension(N), intent(inout) :: out
        integer, intent(hide), depend(data) :: N
    end subroutine reorderlut_f32_a32
F2PY_WRAPPER_END */
void reorderlut_f32_a32(const float *restrict data, uint32_t *restrict lut,
                        float *restrict out, int N) {
    int i;
#pragma omp parallel for
    for (i = 0; i < N; i++) {
        out[i] = data[lut[i]];
    }
}

/* F2PY_WRAPPER_START
    subroutine reorder_u16_a32_a16(data, adr0, adr1, out, ns, nf)
!DOC reorderlut_u16_a32_a16 called in sandbox/fazit.py
!DOC data - source data read in order
!DOC adr0 - output position for the first pixel in each row
!DOC adr1 - difference offset to output next pixel in each row
        intent(c) reorder_u16_a32_a16
        intent(c)
        integer(kind=-2), dimension(ns,nf), intent(in) :: data
        integer*4, dimension(ns), intent(in) :: adr0
        integer*2, dimension(ns,nf), intent(in) :: adr1
        integer(kind=-2), dimension(ns,nf), intent(inout) :: out
        integer, intent(hide), depend(adr1) :: ns = shape(adr1,0)
        integer, intent(hide), depend(adr1) :: nf = shape(adr1,1)
    end subroutine reorder_u16_a32_a16
F2PY_WRAPPER_END */
void reorder_u16_a32_a16(uint16_t *restrict data, uint32_t *restrict a0,
                         int16_t *restrict a1, uint16_t *restrict out, int ns,
                         int nf) {
    int i, j, p;
    /*  printf("Hello, got ns=%d nf=%d\n",ns, nf);*/
#pragma omp parallel for private(p, j)
    for (i = 0; i < ns; i++) {
        p = a0[i];
        for (j = 0; j < nf; j++) {
            p += a1[i * nf + j];
            out[p] = data[i * nf + j];
        }
    }
}

/* F2PY_WRAPPER_START
    subroutine bgcalc(img, bg, msk, ns, nf, gain, sp, st)
!DOC bgcalc computes a background on a 1d signal where gain
!DOC and sp and st are defined by:
!DOC     diff = difference to neighbors or bg estimate
!DOC     sigmap = weight for abs background value
!DOC     sigmat = constant weight
!DOC     gain for b += diff * gain
!DOC img - source data
!DOC bg  - computed background
!DOC msk - mask
        intent(c) bgcalc
        intent(c)
        real, dimension(ns,nf), intent(in) :: img
        real, dimension(ns,nf), intent(inout) :: bg
        integer*1, dimension(ns,nf), intent(inout) :: msk
        integer, intent(hide), depend(img) :: ns = shape(img,0)
        integer, intent(hide), depend(img) :: nf = shape(img,1)
        real :: gain, sp, st
    end subroutine bgcalc
F2PY_WRAPPER_END */
void bgcalc(const float *restrict img, float *restrict bg,
            uint8_t *restrict msk, int ns, int nf, float gain, float sigmap,
            float sigmat) {
    int ir, i;
    float b, diff, t;
    /*
    printf("gain %f sigmap %f sigmat %f ns %d nf %d\n",
    gain, sigmap, sigmat, ns, nf); */
#pragma omp parallel for private(b, diff, t, ir, i)
    for (ir = 0; ir < ns; ir++) { // in range( 1, len(data) ):
        i = ir * nf;
        b = img[i];
        bg[i] = b;
        msk[i] = 1;
        for (i = ir * nf; i < (ir + 1) * nf; i++) {
            diff = img[i] - b;
            t = sigmap * fabsf(b) + sigmat;
            if (diff > t) {
                diff = t * diff / fabsf(diff) / 16;
                msk[i] = 1;
            } else if (diff < -t) {
                diff = t * diff / fabsf(diff) / 4;
                msk[i] = 1;
            } else {
                msk[i] = 0;
            }
            b += diff * gain;
            bg[i] = b;
        }
        i = (ir + 1) * nf - 1;
        b = img[i];
        bg[i] = b;
        msk[i] = 1;
        for (i = (ir + 1) * nf - 1; i >= ir * nf; i--) {
            diff = img[i] - b;
            t = sigmap * fabsf(b) + sigmat;
            if (diff > t) {
                diff = t * diff / fabsf(diff) / 16;
                msk[i] += 1;
            } else if (diff < -t) {
                diff = t * diff / fabsf(diff) / 4;
                msk[i] += 1;
            } else {
                if (msk[i] == 1) {
                    bg[i] = b;
                }
                if ((msk[i] == 0) || (msk[i] == 2)) {
                    bg[i] = (bg[i] + b) / 2;
                }
            }
            b += diff * gain;
        }
    }
}
