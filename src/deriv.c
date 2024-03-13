
#include "ImageD11_cmath.h"
#include <stdio.h>

struct grain {
    double UB[9];
    double t[3];
};

double norm3(double v[3]) {
    return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

/*
def anglevecs2D(a, b, c, d):
    return np.arctan2( a*d - b*c, a*c + b*d)

def derivanglevecs2d( a, b, c, d, ab2=None, cd2=None):
    if ab2 is None:
        ab2 = a*a+b*b
    if cd2 is None:
        cd2 = c*c+d*d
    return {
            'a':  b / ab2,
            'b': -a / ab2,
            'c': -d / cd2,
            'd':  c / cd2,
            }
*/

// axis_system matrix in laboratory co-ordinates is:
//    012 axis[0], axis[1], axis[2]
//    345 b[0] b[1] b[2] = cross product of c and a
//    678 c[0] c[1] c[2] = cross product of axis and beam
//    a.b = dot of axis and beam
//    ax[9]...[19] = inverse transpose matrix for conversion back

int g_to_k(double h[4], double UB[9], double ax[19], double pre[9], double wvln,
           double g[3], double k[3], double tth_calc[10], double eta_calc[10],
           double omega_calc[10]) {

    double gr[3], modg, modg2, ka[3], b, c, t;
    double gao[3], kao[3];
    int i, j;

    matvec(UB, h, g);
    modg = norm3(g);
    t = wvln * modg / 2.0;
    tth_calc[0] = DEG * asin(t);
    // d(tth_c) = DEGREES * tmp4 * dtmp3_dUB
    c = DEG * wvln / (2.0 * modg * sqrt(1 - t * t));
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            tth_calc[i * 3 + j + 1] = c * g[i] * h[j];

    matvec(pre, g, gr); // dgr_dUB
    // gr is g-vector on rotating stage in orthogonal co-ordinate axes
    // scale space to wavelength
    gr[0] = gr[0] * wvln;
    gr[1] = gr[1] * wvln;
    gr[2] = gr[2] * wvln;
    modg2 = gr[0] * gr[0] + gr[1] * gr[1] + gr[2] * gr[2];
    // Now find diffraction vector in axis co-ordinate system
    //  x along axis
    //  y along beam (not orthogonal)
    //  z makes right handed set
    // Component of gr along axis = g.axis
    ka[0] = gr[0] * ax[0] + gr[1] * ax[1] + gr[2] * ax[2]; // AXIS READ
    // Component along beam (rotated) from Laue condition |g|^2/2
    ka[1] = -modg2 / 2.;
    // Solve quadratic for ka[0]
    // |g|^2 = |k|^2
    //       = k[0]*k[0] + k[1]*k[1] + k[2]*k[2] + 2*k[0]*k[1]*ax[9]
    // ax[9] is axis.beam
    c = ka[0] * ka[0] + ka[1] * ka[1] - modg2;
    b = 2 * ka[0] * ka[1] * ax[9];
    // sign to use is in hkl[3]
    t = b * b - c * 4.;
    if (t >= 0.0) {
        ka[2] = -b + h[3] * sqrt(t) / 2.;
    } else {
        printf("Cannot find it!\n");
        return 1;
    }
    // Now we have gr and ka
    // Put ka into laboratory orthogonal system
    matvec(&ax[10], ka, k);
    k[0] = k[0] / wvln;
    k[1] = k[1] / wvln;
    k[2] = k[2] / wvln;
    // eta is found from k...
    // eta = atan2( -ky, kz )
    eta_calc[0] = DEG * atan2(-k[1], k[2]);
    // Now find the diffractometer angle
    // kao = kvector in orthogonal axis system
    kao[0] = ka[0];
    gao[0] = ka[0]; // component on axis

    kao[1] = ka[1] * (1 - ax[9]);
    kao[2] = ka[2];
    gao[1] = ax[3] * gr[0] + ax[4] * gr[1] + ax[5] * gr[2];
    gao[2] = ax[6] * gr[0] + ax[7] * gr[1] + ax[8] * gr[2];
    // return np.arctan2( a*d - b*c, a*c + b*d)
    omega_calc[0] = DEG * atan2(kao[1] * gao[2] - kao[2] * gao[1],
                                kao[1] * gao[1] - kao[2] * gao[2]);
    return 0;
}

double determinant3x3(double m[9]) {
    return m[0] * m[4] * m[8] - m[0] * m[5] * m[7] + m[1] * m[5] * m[6] -
           m[1] * m[3] * m[8] + m[2] * m[3] * m[7] - m[2] * m[4] * m[6];
}

int inverse_mat3(double m[9], double r[9]) {
    double d;
    d = determinant3x3(m);
    if (fabs(d) == 0.0) {
        return 1;
    }
    d = 1.0 / d;
    r[0] = d * (m[4] * m[8] - m[7] * m[5]);
    r[1] = d * (m[7] * m[2] - m[8] * m[1]);
    r[2] = d * (m[1] * m[5] - m[2] * m[4]);
    r[3] = d * (m[5] * m[6] - m[3] * m[8]);
    r[4] = d * (m[8] * m[0] - m[2] * m[6]);
    r[5] = d * (m[2] * m[3] - m[0] * m[5]);
    r[6] = d * (m[3] * m[7] - m[6] * m[4]);
    r[7] = d * (m[1] * m[6] - m[0] * m[7]);
    r[8] = d * (m[0] * m[4] - m[3] * m[1]);
    return 0;
}

int make_ax(double post[9], double axis[3], double beam[3], double ax[19]) {
    double rb[3];
    ax[0] = axis[0];
    ax[1] = axis[1];
    ax[2] = axis[2];
    matTvec(post, beam, rb);
    crossProduct(axis, rb, &ax[6]);
    crossProduct(axis, &ax[6], &ax[3]);
    ax[9] = ax[0] * ax[3] + ax[1] * ax[4] + ax[2] * ax[5];
    return inverse_mat3(&ax[0], &ax[10]);
}

int main() {
    double tth_calc[10], tth_calc1[10];
    double omega_calc[10], omega_calc1[10];
    double eta_calc[10], eta_calc1[10];
    int i;
    double p, wvln, testval, h[4];
    double post[9], pre[9], axis[3], beam[3], ax[19];
    double gcalc[3], kcalc[3];
    struct grain gr, gr_test;

    wvln = 0.12345;
    for (i = 0; i < 9; i++) {
        gr.UB[i] = 0;
        post[i] = 0;
        pre[i] = 0;
    }
    gr.UB[0] = 0.22;
    gr.UB[4] = 0.57;
    gr.UB[8] = 0.42;
    h[0] = 5;
    h[1] = 3;
    h[2] = 2;
    h[3] = 1;
    post[0] = 1.;
    post[4] = 1.;
    post[8] = 1.;
    pre[0] = 1.;
    pre[4] = 1.;
    pre[8] = 1.;
    axis[0] = 0.;
    axis[1] = 0.;
    axis[2] = 1.;
    beam[0] = 1.;
    beam[1] = 0.;
    beam[2] = 0.;
    make_ax(post, axis, beam, ax);
    for (i = 0; i < 19; i++)
        printf("ax[%d]=%f\n", i, ax[i]);

    g_to_k(h, gr.UB, ax, pre, wvln, gcalc, kcalc, tth_calc, eta_calc,
           omega_calc);
    printf("gcalc %f %f %f |%f|\n", gcalc[0], gcalc[1], gcalc[2], norm3(gcalc));
    printf("kcalc %f %f %f |%f|\n", kcalc[0], kcalc[1], kcalc[2], norm3(kcalc));
    printf("eta_calc %f omega_calc %f tth_calc %f\n", eta_calc[0],
           omega_calc[0], tth_calc[0]);

    p = 1e-8;
    for (i = 0; i < 9; i++) {
        // Now test
        gr_test = gr;
        gr_test.UB[i] += p;
        g_to_k(h, gr_test.UB, ax, pre, wvln, gcalc, kcalc, tth_calc1, eta_calc1,
               omega_calc1);
        testval = (tth_calc1[0] - tth_calc[0]) / p;
        printf(" dtth_dUB[%d] = %f ", i, tth_calc[i + 1]);
        printf(" %f \n", testval);
    }
}
