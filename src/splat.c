
#include "cImageD11.h"
#include <stdio.h>

void splat( uint8_t rgba[], int w, int h, 
            double gve[][3], int ng, double u[9],
            int npx
            ){
    int i, j,k, imx, imy, imz;
    /* init */
    for(i=0;i<w*h*4;i++) {
        rgba[i]=255;
    }
    /* not parallel as rgba is shared */
    for(i=0;i<ng;i++){
        imx = (int) (( u[0]*gve[i][0] + u[1]*gve[i][1] + u[2]*gve[i][2])*h/2. + w/2.);
        imy = (int) (( u[3]*gve[i][0] + u[4]*gve[i][1] + u[5]*gve[i][2])*h/2. + h/2.);
        imz = (int) (( u[6]*gve[i][0] + u[7]*gve[i][1] + u[8]*gve[i][2])*64. + 128.);
        if( (imx > npx) && (imx < w-npx) && (imy > npx) && (imy < h-npx) && (imz > 64) && (imz < 192) ){
            // 4*(imx*w+imy)
            for( j = -npx; j<= npx; j++){
                for( k = -npx; k<=npx ; k++){
                rgba[ w*(imy+j)*4 + (imx+k)*4 + 0] = 0;
                rgba[ w*(imy+j)*4 + (imx+k)*4 + 1] = 0;
                rgba[ w*(imy+j)*4 + (imx+k)*4 + 2] = 0;
                rgba[ w*(imy+j)*4 + (imx+k)*4 + 3] = imz;
                }
            }
        }
    }
}