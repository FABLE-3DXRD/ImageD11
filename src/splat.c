
#include "cImageD11.h"
#include <stdio.h>


/* TODO next:
 * set the color and markersize per peak
 * perhaps also a draw order (back to front, top to bottom) ?
 */

void splat( uint8_t rgba[], int w, int h,
            double gve[][3], int ng, double u[9],
            int npx
            ){
    int32_t i, j, k, imx, imy, imz, w2, h2;
    double s[9];

    /* init */
    h2 = h/2;
    w2 = w/2;
    for(i=0;i<6;i++){
        s[i] = u[i]*((w+h)/4);
    }
    s[6] = u[6]*64;
    s[7] = u[7]*64;
    s[8] = u[8]*64;
    /* Not parallel - seems to be fast anyway and rgba is shared on write */
    for(i=0;i<w*h*4;i=i+4) {
        rgba[i]  =0;
        rgba[i+1]=0;
        rgba[i+2]=0;
        rgba[i+3]=255;
    }
    for(i=0;i<ng;i++){
        imx = (s[0]*gve[i][0] + s[1]*gve[i][1] + s[2]*gve[i][2]) + w2;
        imy = (s[3]*gve[i][0] + s[4]*gve[i][1] + s[5]*gve[i][2]) + h2;
        imz = (s[6]*gve[i][0] + s[7]*gve[i][1] + s[8]*gve[i][2]) + 128;
        if( (imx > npx) && (imx < w-npx) && (imy > npx) && (imy < h-npx) && (imz >= 0) && (imz < 256) ){
            for( j = -npx; j<= npx; j++){
                for( k = -npx; k<=npx ; k++){
                    rgba[ w*(imy+j)*4 + (imx+k)*4 + 0] = 255;
                    rgba[ w*(imy+j)*4 + (imx+k)*4 + 1] = 255;
                    rgba[ w*(imy+j)*4 + (imx+k)*4 + 2] = 255;
                    rgba[ w*(imy+j)*4 + (imx+k)*4 + 3] = imz;
                }
            }
        }
    }
}