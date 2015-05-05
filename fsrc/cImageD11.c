

#include <math.h>
#include <string.h>

#define PI 3.14159265358979323846
#define RAD PI/180.0
#define DEG 180.0/PI


#define vec3sub(a,b,c) \
        (c)[0] = (a)[0] - (b)[0]; \
        (c)[1] = (a)[1] - (b)[1]; \
        (c)[2] = (a)[2] - (b)[2]; 


#define matvec(a,b,c) \
    (c)[0] = (a)[0]*(b)[0] + (a)[1]*(b)[1] + (a)[2]*(b)[2]; \
    (c)[1] = (a)[3]*(b)[0] + (a)[4]*(b)[1] + (a)[5]*(b)[2]; \
    (c)[2] = (a)[6]*(b)[0] + (a)[7]*(b)[1] + (a)[8]*(b)[2]; 


#define matTvec(a,b,c) \
    (c)[0] = (a)[0]*(b)[0] + (a)[3]*(b)[1] + (a)[6]*(b)[2]; \
    (c)[1] = (a)[1]*(b)[0] + (a)[4]*(b)[1] + (a)[7]*(b)[2]; \
    (c)[2] = (a)[2]*(b)[0] + (a)[5]*(b)[1] + (a)[8]*(b)[2]; 


#define matmat(a,b,c) \
    (c)[0] = (a)[0]*(b)[0] + (a)[3]*(b)[1] + (a)[6]*(b)[2]; \
    (c)[1] = (a)[1]*(b)[0] + (a)[4]*(b)[1] + (a)[7]*(b)[2]; \
    (c)[2] = (a)[2]*(b)[0] + (a)[5]*(b)[1] + (a)[8]*(b)[2]; \
    (c)[3] = (a)[0]*(b)[3] + (a)[3]*(b)[4] + (a)[6]*(b)[5]; \
    (c)[4] = (a)[1]*(b)[3] + (a)[4]*(b)[4] + (a)[7]*(b)[5]; \
    (c)[5] = (a)[2]*(b)[3] + (a)[5]*(b)[4] + (a)[8]*(b)[5]; \
    (c)[6] = (a)[0]*(b)[6] + (a)[3]*(b)[7] + (a)[6]*(b)[8]; \
    (c)[7] = (a)[1]*(b)[6] + (a)[4]*(b)[7] + (a)[7]*(b)[8]; \
    (c)[8] = (a)[2]*(b)[6] + (a)[5]*(b)[7] + (a)[8]*(b)[8]; 





void assign( double ubi[9], double gv[][3], double tol,
	     double drlv2[], int labels[], int ig, int n);


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

/*
! compute_gv for refinegrains.py

subroutine compute_gv( xlylzl, omega, omegasign, wvln, wedge, chi, t, gv, n )
  use omp_lib
  implicit none
  real, intent(in) :: xlylzl(3,n), omega(n), wvln, wedge, chi, t(3)
  real(8), intent(inout):: gv(3,n)
  integer, intent(in) ::  n, omegasign
*/ 


void compute_gv( double xlylzl[][3], double omega[], double omegasign,
		 double wvln, double wedge, double chi, double t[3],
		 double gv[][3], int n);

void compute_gv( double xlylzl[][3], double omega[], double omegasign,
		 double wvln, double wedge, double chi, double t[3],
		 double gv[][3], int n){
  /*  real :: sc,cc,sw,cw,wmat(3,3),cmat(3,3), mat(3,3), u(3),d(3),v(3)
      real :: modyz, o(3), co, so, ds, k(3)
      real, parameter :: PI=3.141592653589793,RAD=PI/180.0,DEG=180.0/PI
      integer :: i
  */
  double sc, cc, sw, cw, wmat[9], cmat[9], mat[9], u[3], d[3], v[3];
  double modyz, o[3], co, so, ds, k[3];
  int i;
  // ! Fill in rotation matrix of wedge, chi
  sw = sin(wedge*RAD);
  cw = cos(wedge*RAD);
  sc = sin(chi*RAD);
  cc = cos(chi*RAD);
  memcpy(wmat, (double[9])  {cw, 0., -sw, 0., 1., 0., sw, 0., cw},
	 9*sizeof(double)); 
  memcpy(cmat, (double[9])  {1., 0.,  0., 0.,cc, sc , 0.,-sc, cc},
	 9*sizeof(double)); 
  //	mat = matmul(cmat, wmat)
  matmat(cmat, wmat, mat);
  //!  write(*,*)'threads',omp_get_max_threads()
  //!  write(*,*)'mat',mat
  //!$omp parallel do private(so,co,u,o,d,modyz,ds,v,k)
#pragma omp parallel for private(so,co,u,o,d,modyz,ds,v,k)
  //do i=1,n
  for(i=0;i<n;i++){
    // ! Compute translation + rotation for grain origin
    so = sin(RAD*omega[i]*omegasign);
    co = cos(RAD*omega[i]*omegasign);
    u[0] =  co*t[0] - so*t[1];
    u[1] =  so*t[0] + co*t[1];
    u[2] = t[2];
    //! grain origin, difference vec, |yz| component
    // ! o = matmul(transpose(mat),u)
    matTvec( mat, u, o);
    //     d = xlylzl(:,i) - o
    vec3sub( xlylzl[i], o, d);
    //modyz  = 1./sqrt(d(1)*d(1) + d(2)*d(2) + d(3)*d(3))
    modyz =  1./sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
    //     ! k-vector
    ds = 1./wvln;
    k[0] = ds*(d[0]*modyz - 1. );
    k[1] = ds*d[1]*modyz;
    k[2] = ds*d[2]*modyz;
    //v = matmul(mat,k)
    matmat( mat, k, v);
    gv[i][0]= co*v[0] + so*v[1];
    gv[i][1]=-so*v[0] + co*v[1];
    gv[i][2]=v[2];
  }//  enddo
} // end subroutine compute_gv


/*
subroutine compute_xlylzl( s, f, p, r, dist, xlylzl, n )
! Computes laboratory co-ordinates
!
!
  implicit none
  real(8), intent(in) :: s(n), f(n)
  real(8), intent(in) :: p(4), r(3,3), dist(3)
  real(8), intent(inout):: xlylzl(3,n)
  integer, intent(in) ::  n
  ! parameters
  real(8) :: s_cen, f_cen, s_size, f_size
  ! temporaries
  real(8) :: v(3)
  integer :: i, j
  ! unpack parameters
  s_cen = p(1)
  f_cen = p(2)
  s_size = p(3)
  f_size = p(4)
  do i = 1, n
     ! Place on the detector plane accounting for centre and size
     ! subtraction of centre is done here and not later for fear of
     ! rounding errors
     v(1) = 0.0d0
     v(2) = (f(i) - f_cen)*f_size
     v(3) = (s(i) - s_cen)*s_size
     ! Apply the flip and rotation, python was :
     ! fl = dot( [[o11, o12], [o21, o22]], peaks=[[z],[y]] )
     ! vec = [0,fl[1],fl[0]]
     ! return dist + dot(rotmat, vec)
     do j = 1, 3
        ! Skip as v(1) is zero : r(1,j)*v(1)
        xlylzl(j,i) = r(2,j)*v(2) + r(3,j)*v(3) + dist(j)
     enddo
  enddo
end subroutine compute_xlylzl


! set LDFLAGS="-static-libgfortran -static-libgcc -static -lgomp -shared"  
! f2py -m fImageD11 -c fImageD11.f90 --opt=-O3 --f90flags="-fopenmp" -lgomp -lpthread
! export OMP_NUM_THREADS=12
! python tst.py ../test/nac_demo/peaks.out_merge_t200 ../test/nac_demo/nac.prm 
! python test_xlylzl.py ../test/nac_demo/peaks.out_merge_t200 ../test/nac_demo/nac.prm

! f2py -m fImageD11 -c fImageD11.f90 --f90flags="-fopenmp" -lgomp -lpthread  --fcompiler=gnu95 --compiler=mingw32 -DF2PY_REPORT_ON_ARRAY_COPY=1


gcc -fopenmp -Wall -c cImageD11.c

*/

