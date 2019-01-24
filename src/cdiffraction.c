
#include <omp.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include "cdiffraction.h"

#define PI 3.14159265358979323846
#define RAD PI/180.0
#define DEG 180.0/PI

#define NOISY 0

void assign( double ubi[9], double gv[][3], double tol,
	     double drlv2[], int labels[], int ig, int n){
  int i;
  double dr, h[3], ttol, dh[3];
  ttol = tol * tol;
#pragma omp parallel for private(h,dr,dh)
  for(i=0;i<n;i++){
    matvec(ubi, gv[i], h)
    dh[0] = fabs(floor(h[0]+0.5) - h[0]);
    dh[1] = fabs(floor(h[1]+0.5) - h[1]);
    dh[2] = fabs(floor(h[2]+0.5) - h[2]);
    dr =  dh[0] + dh[1] + dh[2];
    dr = dr * dr;
    if ( (dr < ttol) && (dr < drlv2[i]) ) {
      drlv2[i] = dr;
      labels[i] = ig;
    }  else if (labels[i] == ig) {
      labels[i] = -1;
    } // end if
  } //    end for
}  // end function assign


void compute_gv( double xlylzl[][3], double omega[], double omegasign,
		 double wvln, double wedge, double chi, double t[3],
		 double gv[][3], int n){
  double sc, cc, sw, cw, wmat[9], cmat[9], mat[9], u[3], d[3], v[3];
  double modyz, o[3], co, so, ds, k[3];
  int i;
  // ! Fill in rotation matrix of wedge, chi
  sw = sin(wedge*RAD);
  cw = cos(wedge*RAD);
  wmat[0] = cw; wmat[1]=0.0; wmat[2]=-sw;
  wmat[3] = 0.; wmat[4]=1.0; wmat[5]= 0.;
  wmat[6] = sw; wmat[7]=0.0; wmat[8]= cw;
  sc = sin(chi*RAD);
  cc = cos(chi*RAD);
  cmat[0] = 1.; cmat[1]=0.0; cmat[2]= 0.;
  cmat[3] = 0.; cmat[4]= cc; cmat[5]=-sc;
  cmat[6] = 0.; cmat[7]= sc; cmat[8]= cc;
  // Combined mat = chi.wedge
  matmat(cmat, wmat, mat);
#pragma omp parallel for private(so,co,u,o,d,modyz,ds,v,k)
  for(i=0;i<n;i++){
    // ! Compute translation + rotation for grain origin
    so = sin(RAD*omega[i]*omegasign);
    co = cos(RAD*omega[i]*omegasign);
    // Omega matrix vector on translation
    u[0] =  co*t[0] - so*t[1];
    u[1] =  so*t[0] + co*t[1];
    u[2] = t[2];
    // grain origin, difference vec, |yz| component
    matvec( mat, u, o);
    // d is difference vector 
    vec3sub( xlylzl[i], o, d);
    modyz =  1./sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
    //     ! k-vector
    ds = 1./wvln;
    k[0] = ds*(d[0]*modyz - 1. );
    k[1] = ds*d[1]*modyz;
    k[2] = ds*d[2]*modyz;
    matTvec( mat, k, v);
    // Forwards rotation with omega finally
    gv[i][0]= co*v[0] + so*v[1];
    gv[i][1]=-so*v[0] + co*v[1];
    gv[i][2]=v[2];
  }//  enddo
} // end subroutine compute_gv



void compute_xlylzl( double s[], double f[], double p[4],
		     double r[9], double dist[3],
		     double xlylzl[][3], int n){
  double s_cen, f_cen, s_size, f_size, v[3];
  int i, j;
  s_cen = p[0];
  f_cen = p[1];
  s_size = p[2];
  f_size = p[3];
  v[0] = 0.0;
  if(NOISY){
    printf("s_cen %f f_cen %f s_size %f f_size %f\n",s_cen, f_cen, s_size, f_size);
    for(j=0;j<3;j++) printf("dist[%d]=%f ",j,dist[j]);
    for(j=0;j<9;j++) printf("r[%d]=%f ",j,r[j]);
    printf("\n");
  }
  for(i = 0; i < n; i++){
    //     ! Place on the detector plane accounting for centre and size
    //     ! subtraction of centre is done here and not later for fear of
    //     ! rounding errors
    
    v[1] = (f[i] - f_cen)*f_size;
    v[2] = (s[i] - s_cen)*s_size;
    // ! Apply the flip and rotation, python was :
    // ! fl = dot( [[o11, o12], [o21, o22]], peaks=[[z],[y]] )
    // ! vec = [0,fl[1],fl[0]]
    // ! return dist + dot(rotmat, vec)
       for(j = 0; j<3; j++){
	 //  ! Skip as v[0] is zero : r(1,j)*v(1)
	 xlylzl[i][j] = r[3*j+1]*v[1] + r[3*j+2]*v[2] + dist[j];
       } // enddo
  }// enddo
} // end subroutine compute_xlylzl

/* 
! set LDFLAGS="-static-libgfortran -static-libgcc -static -lgomp -shared"  
! f2py -m fImageD11 -c fImageD11.f90 --opt=-O3 --f90flags="-fopenmp" -lgomp -lpthread
! export OMP_NUM_THREADS=12
! python tst.py ../test/nac_demo/peaks.out_merge_t200 ../test/nac_demo/nac.prm 
! python test_xlylzl.py ../test/nac_demo/peaks.out_merge_t200 ../test/nac_demo/nac.prm

! f2py -m fImageD11 -c fImageD11.f90 --f90flags="-fopenmp" -lgomp -lpthread  --fcompiler=gnu95 --compiler=mingw32 -DF2PY_REPORT_ON_ARRAY_COPY=1


gcc -fopenmp -Wall -c cImageD11.c

*/

