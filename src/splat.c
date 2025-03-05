
#include "cImageD11.h"
#include <stdio.h>

/* F2PY_WRAPPER_START
    subroutine splat( rgba, w, h, gve, ng, u, npx )
!DOC splat draws gvectors into an rgba image. The horror of maintaining plot3d
!DOC over the years motivated this code. See test/demo/tksplat
!DOC * set the color and markersize per peak
!DOC * perhaps also a draw order (back to front, top to bottom) ?
        intent(c) splat
        intent(c)
        integer(kind=-1), dimension(w,h,4) :: rgba
        integer,intent(hide), depend(rgba) :: h = shape(rgba, 0)
        integer,intent(hide), depend(rgba) :: w = shape(rgba, 1)
        double precision, intent(inout), dimension(ng,3) ::  gve
        integer, intent(hide), depend(gve):: ng
        double precision, intent(in), dimension(9) :: u
        integer :: npx
        ! NOT threadsafe, rgba will be recycled and on screen
    end subroutine splat

F2PY_WRAPPER_END */
void splat(uint8_t rgba[], int w, int h, double gve[][3], int ng, double u[9],
           int npx) {
    int32_t i, j, k, imx, imy, imz, w2, h2;
    double s[9];

    /* init */
    h2 = h / 2;
    w2 = w / 2;
    for (i = 0; i < 6; i++) {
        s[i] = u[i] * ((w + h) / 4);
    }
    s[6] = u[6] * 64;
    s[7] = u[7] * 64;
    s[8] = u[8] * 64;
    /* Not parallel - seems to be fast anyway and rgba is shared on write */
    for (i = 0; i < w * h * 4; i = i + 4) {
        rgba[i] = 0;
        rgba[i + 1] = 0;
        rgba[i + 2] = 0;
        rgba[i + 3] = 255;
    }
    for (i = 0; i < ng; i++) {
        imx =
            (int)(s[0] * gve[i][0] + s[1] * gve[i][1] + s[2] * gve[i][2]) + w2;
        imy =
            (int)(s[3] * gve[i][0] + s[4] * gve[i][1] + s[5] * gve[i][2]) + h2;
        imz =
            (int)(s[6] * gve[i][0] + s[7] * gve[i][1] + s[8] * gve[i][2]) + 128;
        if ((imx > npx) && (imx < w - npx) && (imy > npx) && (imy < h - npx) &&
            (imz >= 0) && (imz < 256)) {
            for (j = -npx; j <= npx; j++) {
                for (k = -npx; k <= npx; k++) {
                    rgba[w * (imy + j) * 4 + (imx + k) * 4 + 0] = 255;
                    rgba[w * (imy + j) * 4 + (imx + k) * 4 + 1] = 255;
                    rgba[w * (imy + j) * 4 + (imx + k) * 4 + 2] = 255;
                    rgba[w * (imy + j) * 4 + (imx + k) * 4 + 3] = imz;
                }
            }
        }
    }
}