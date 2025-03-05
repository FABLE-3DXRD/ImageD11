
#include "cdiffraction.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

#define PI 3.14159265358979323846
#define RAD PI / 180.0
#define DEG 180.0 / PI

#define NOISY 0

/* F2PY_WRAPPER_START
    subroutine compute_geometry(xlylzl,omega,omegasign,wvln,wedge,chi,t,out,ng)
!DOC compute_geometry is for the "updateGeometry" method of columnfiles
!DOC from xlylzl it will compute tth, eta, ds, gve into out
!DOC in the laboratory in xlylzl[npks], the omega rotation[npks], and
!DOC the rest of the parameters (wedge,wvln,chi,t[3] and omegasign)
!DOC out should contain : (tth, eta, ds, gx, gy, gz)
        intent(c) compute_geometry
        intent(c)
        integer, intent(c,hide), depend( xlylzl ) :: ng
        double precision, intent(in):: xlylzl(ng,3)
        double precision, intent(in):: omega(ng)
        double precision, intent(in):: omegasign, wvln, wedge, chi
        double precision, intent(in):: t(3)
        double precision, intent(inout):: out(ng,6)
        threadsafe
    end subroutine compute_geometry
F2PY_WRAPPER_END */
void compute_geometry(double xlylzl[][3], double omega[], double omegasign,
                      double wvln, double wedge, double chi, double t[3],
                      double out[][6], int n) {
    double sc, cc, sw, cw, wmat[9], cmat[9], mat[9], u[3], d[3], v[3];
    double modyz, o[3], co, so, ds, k[3];
    int i;
    // ! Fill in rotation matrix of wedge, chi
    sw = sin(wedge * RAD);
    cw = cos(wedge * RAD);
    wmat[0] = cw;
    wmat[1] = 0.0;
    wmat[2] = -sw;
    wmat[3] = 0.;
    wmat[4] = 1.0;
    wmat[5] = 0.;
    wmat[6] = sw;
    wmat[7] = 0.0;
    wmat[8] = cw;
    sc = sin(chi * RAD);
    cc = cos(chi * RAD);
    cmat[0] = 1.;
    cmat[1] = 0.0;
    cmat[2] = 0.;
    cmat[3] = 0.;
    cmat[4] = cc;
    cmat[5] = -sc;
    cmat[6] = 0.;
    cmat[7] = sc;
    cmat[8] = cc;
    // Combined mat = chi.wedge
    matmat(cmat, wmat, mat);
#pragma omp parallel for private(so, co, u, o, d, modyz, ds, v, k)
    for (i = 0; i < n; i++) {
        // ! Compute translation + rotation for grain origin
        so = sin(RAD * omega[i] * omegasign);
        co = cos(RAD * omega[i] * omegasign);
        // Omega matrix vector on translation
        u[0] = co * t[0] - so * t[1];
        u[1] = so * t[0] + co * t[1];
        u[2] = t[2];
        // o=grain origin
        matvec(mat, u, o);
        // d is difference vector
        vec3sub(xlylzl[i], o, d);
        modyz = 1. / sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        // two theta
        out[i][0] = DEG * atan2(sqrt(d[1] * d[1] + d[2] * d[2]), d[0]);
        //     ! k-vector
        ds = 1. / wvln;
        k[0] = ds * (d[0] * modyz - 1.);
        k[1] = ds * d[1] * modyz;
        k[2] = ds * d[2] * modyz;
        // eta
        out[i][1] = DEG * atan2(-d[1], d[2]);
        // dstar
        out[i][2] = sqrt(k[0] * k[0] + k[1] * k[1] + k[2] * k[2]);
        // g-vector
        matTvec(mat, k, v);
        // Forwards rotation with omega finally
        out[i][3] = co * v[0] + so * v[1];
        out[i][4] = -so * v[0] + co * v[1];
        out[i][5] = v[2];
    } //  enddo
} // end subroutine compute_geometry

/* F2PY_WRAPPER_START
    subroutine compute_gv(xlylzl,omega,omegasign,wvln,wedge,chi,t,gv,ng)
!DOC compute_gv computes scattering vectors given thr positions of the spot
!DOC in the laboratory in xlylzl[npks], the omega rotation[npks], and
!DOC the rest of the parameters (wedge,wvln,chi,t[3] and omegasign)
        intent(c) compute_gv
        intent(c)
        integer, intent(c,hide), depend( xlylzl ) :: ng
        double precision, intent(in):: xlylzl(ng,3)
        double precision, intent(in):: omega(ng)
        double precision, intent(in):: omegasign, wvln, wedge, chi
        double precision, intent(in):: t(3)
        double precision, intent(inout):: gv(ng,3)
        ! NOT threadsafe since gv may be shared
    end subroutine compute_gv
F2PY_WRAPPER_END */
void compute_gv(double xlylzl[][3], double omega[], double omegasign,
                double wvln, double wedge, double chi, double t[3],
                double gv[][3], int n) {
    double sc, cc, sw, cw, wmat[9], cmat[9], mat[9], u[3], d[3], v[3];
    double modyz, o[3], co, so, ds, k[3];
    int i;
    // ! Fill in rotation matrix of wedge, chi
    sw = sin(wedge * RAD);
    cw = cos(wedge * RAD);
    wmat[0] = cw;
    wmat[1] = 0.0;
    wmat[2] = -sw;
    wmat[3] = 0.;
    wmat[4] = 1.0;
    wmat[5] = 0.;
    wmat[6] = sw;
    wmat[7] = 0.0;
    wmat[8] = cw;
    sc = sin(chi * RAD);
    cc = cos(chi * RAD);
    cmat[0] = 1.;
    cmat[1] = 0.0;
    cmat[2] = 0.;
    cmat[3] = 0.;
    cmat[4] = cc;
    cmat[5] = -sc;
    cmat[6] = 0.;
    cmat[7] = sc;
    cmat[8] = cc;
    // Combined mat = chi.wedge
    matmat(cmat, wmat, mat);
#pragma omp parallel for private(so, co, u, o, d, modyz, ds, v, k)
    for (i = 0; i < n; i++) {
        // ! Compute translation + rotation for grain origin
        so = sin(RAD * omega[i] * omegasign);
        co = cos(RAD * omega[i] * omegasign);
        // Omega matrix vector on translation
        u[0] = co * t[0] - so * t[1];
        u[1] = so * t[0] + co * t[1];
        u[2] = t[2];
        // grain origin, difference vec, |yz| component
        matvec(mat, u, o);
        // d is difference vector
        vec3sub(xlylzl[i], o, d);
        modyz = 1. / sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        //     ! k-vector
        ds = 1. / wvln;
        k[0] = ds * (d[0] * modyz - 1.);
        k[1] = ds * d[1] * modyz;
        k[2] = ds * d[2] * modyz;
        matTvec(mat, k, v);
        // Forwards rotation with omega finally
        gv[i][0] = co * v[0] + so * v[1];
        gv[i][1] = -so * v[0] + co * v[1];
        gv[i][2] = v[2];
    } //  enddo
} // end subroutine compute_gv

/* F2PY_WRAPPER_START
    subroutine compute_xlylzl(s,f,p,r,dist,xlylzl,n)
!DOC computes_xlylzl finds spot positions in the laboratory frame
!DOC using packed parameters that are more general
!DOC s    = slow pixel position
!DOC f    = fast pixel position
!DOC p    = [s_cen, f_cen, s_size, f_size]
!DOC r[9] = dot( transform.detector_rotation_matrix, flipmatrix )
!DOC ... with flipmatrix = [1,0,0], [0,o22,o21], [0,o12,o11]
!DOC dist = [distancex, distancey, distancez] is 3D (beyond old model)
!DOC ... see test/test_cImageD11.py
        intent(c) compute_xlylzl
        intent(c)
        double precision, intent(in):: s(n), f(n)
        double precision, intent(inout):: xlylzl(n,3)
        double precision, intent(in):: p(4), r(9), dist(3)
        integer, intent(c,hide), depend( s ) :: n
        ! NOT threadsafe since xl may be shared
    end subroutine compute_xlylzl
F2PY_WRAPPER_END */
void compute_xlylzl(double s[], double f[], double p[4], double r[9],
                    double dist[3], double xlylzl[][3], int n) {
    double s_cen, f_cen, s_size, f_size, v[3];
    int i, j;
    s_cen = p[0];
    f_cen = p[1];
    s_size = p[2];
    f_size = p[3];
    v[0] = 0.0;
    if (NOISY) {
        printf("s_cen %f f_cen %f s_size %f f_size %f\n", s_cen, f_cen, s_size,
               f_size);
        for (j = 0; j < 3; j++)
            printf("dist[%d]=%f ", j, dist[j]);
        for (j = 0; j < 9; j++)
            printf("r[%d]=%f ", j, r[j]);
        printf("\n");
    }
    for (i = 0; i < n; i++) {
        //     ! Place on the detector plane accounting for centre and size
        //     ! subtraction of centre is done here and not later for fear of
        //     ! rounding errors

        v[1] = (f[i] - f_cen) * f_size;
        v[2] = (s[i] - s_cen) * s_size;
        // ! Apply the flip and rotation, python was :
        // ! fl = dot( [[o11, o12], [o21, o22]], peaks=[[z],[y]] )
        // ! vec = [0,fl[1],fl[0]]
        // ! return dist + dot(rotmat, vec)
        for (j = 0; j < 3; j++) {
            //  ! Skip as v[0] is zero : r(1,j)*v(1)
            xlylzl[i][j] = r[3 * j + 1] * v[1] + r[3 * j + 2] * v[2] + dist[j];
        } // enddo
    } // enddo
} // end subroutine compute_xlylzl

/* F2PY_WRAPPER_START
    subroutine quickorient( ubi, bt )
!DOC quickorient takes two g-vectors in UBI[0] and UBI[1]
!DOC and overwrites with UBI orientation using cache in bt (from h1,h2)
!DOC ... computes cross product 0x1 = ubi[0]xubi[1]
!DOC ... normalises u0=ubi0 and u2=0x1 and computes u1=u0xu2
!DOC ... returns in UBI BT.(u1,u2,u3)
!DOC algorithm was due to Busing and Levy
    intent(c) quickorient
        intent(c)
        double precision, dimension(3,3), intent (inout) :: ubi
        double precision, dimension(3,3), intent (in) :: bt
    end subroutine orient
F2PY_WRAPPER_END */
void quickorient(double UBI[9], double BT[9]) {
    /* On entry UBI[0] is g1, UBI[1] is g2, BT is made for this to work.
       0 1 2  == g1
       3 4 5  == g2
       6 7 8  <- g1xg2
     */
    double t0, t1, M[9];
    /* g2 = g0xg1 */
    M[6] = UBI[1] * UBI[5] - UBI[2] * UBI[4];
    M[7] = UBI[2] * UBI[3] - UBI[0] * UBI[5];
    M[8] = UBI[0] * UBI[4] - UBI[1] * UBI[3];
    /* u0 = norm(g0) */
    t0 = sqrt(UBI[0] * UBI[0] + UBI[1] * UBI[1] + UBI[2] * UBI[2]);
    M[0] = UBI[0] / t0;
    M[1] = UBI[1] / t0;
    M[2] = UBI[2] / t0;
    /* u3 = norm(g3) */
    t1 = sqrt(M[6] * M[6] + M[7] * M[7] + M[8] * M[8]);
    M[6] /= t1;
    M[7] /= t1;
    M[8] /= t1;
    /* u2 = u1xu3 */
    M[3] = M[1] * M[8] - M[2] * M[7];
    M[4] = M[2] * M[6] - M[0] * M[8];
    M[5] = M[0] * M[7] - M[1] * M[6];
    /* ubi = dot( BT, (u1,u2,u2) */
    UBI[0] = BT[0] * M[0] + BT[1] * M[3] + BT[2] * M[6];
    UBI[1] = BT[0] * M[1] + BT[1] * M[4] + BT[2] * M[7];
    UBI[2] = BT[0] * M[2] + BT[1] * M[5] + BT[2] * M[8];
    UBI[3] = BT[3] * M[0] + BT[4] * M[3] + BT[5] * M[6];
    UBI[4] = BT[3] * M[1] + BT[4] * M[4] + BT[5] * M[7];
    UBI[5] = BT[3] * M[2] + BT[4] * M[5] + BT[5] * M[8];
    UBI[6] = BT[6] * M[0] + BT[7] * M[3] + BT[8] * M[6];
    UBI[7] = BT[6] * M[1] + BT[7] * M[4] + BT[8] * M[7];
    UBI[8] = BT[6] * M[2] + BT[7] * M[5] + BT[8] * M[8];
}

/*
! set LDFLAGS="-static-libgfortran -static-libgcc -static -lgomp -shared"
! f2py -m fImageD11 -c fImageD11.f90 --opt=-O3 --f90flags="-fopenmp" -lgomp
-lpthread ! export OMP_NUM_THREADS=12 ! python tst.py
../test/nac_demo/peaks.out_merge_t200 ../test/nac_demo/nac.prm ! python
test_xlylzl.py ../test/nac_demo/peaks.out_merge_t200 ../test/nac_demo/nac.prm

! f2py -m fImageD11 -c fImageD11.f90 --f90flags="-fopenmp" -lgomp -lpthread
--fcompiler=gnu95 --compiler=mingw32 -DF2PY_REPORT_ON_ARRAY_COPY=1

gcc -fopenmp -Wall -c cImageD11.c

*/
