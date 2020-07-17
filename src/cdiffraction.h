#define PI 3.14159265358979323846
#define RAD PI / 180.0
#define DEG 180.0 / PI

#define vec3sub(a, b, c)                                                       \
    (c)[0] = (a)[0] - (b)[0];                                                  \
    (c)[1] = (a)[1] - (b)[1];                                                  \
    (c)[2] = (a)[2] - (b)[2];

#define matvec(a, b, c)                                                        \
    (c)[0] = (a)[0] * (b)[0] + (a)[1] * (b)[1] + (a)[2] * (b)[2];              \
    (c)[1] = (a)[3] * (b)[0] + (a)[4] * (b)[1] + (a)[5] * (b)[2];              \
    (c)[2] = (a)[6] * (b)[0] + (a)[7] * (b)[1] + (a)[8] * (b)[2];

#define matTvec(a, b, c)                                                       \
    (c)[0] = (a)[0] * (b)[0] + (a)[3] * (b)[1] + (a)[6] * (b)[2];              \
    (c)[1] = (a)[1] * (b)[0] + (a)[4] * (b)[1] + (a)[7] * (b)[2];              \
    (c)[2] = (a)[2] * (b)[0] + (a)[5] * (b)[1] + (a)[8] * (b)[2];

#define matmat(a, b, c)                                                        \
    (c)[0] = (a)[0] * (b)[0] + (a)[3] * (b)[1] + (a)[6] * (b)[2];              \
    (c)[1] = (a)[1] * (b)[0] + (a)[4] * (b)[1] + (a)[7] * (b)[2];              \
    (c)[2] = (a)[2] * (b)[0] + (a)[5] * (b)[1] + (a)[8] * (b)[2];              \
    (c)[3] = (a)[0] * (b)[3] + (a)[3] * (b)[4] + (a)[6] * (b)[5];              \
    (c)[4] = (a)[1] * (b)[3] + (a)[4] * (b)[4] + (a)[7] * (b)[5];              \
    (c)[5] = (a)[2] * (b)[3] + (a)[5] * (b)[4] + (a)[8] * (b)[5];              \
    (c)[6] = (a)[0] * (b)[6] + (a)[3] * (b)[7] + (a)[6] * (b)[8];              \
    (c)[7] = (a)[1] * (b)[6] + (a)[4] * (b)[7] + (a)[7] * (b)[8];              \
    (c)[8] = (a)[2] * (b)[6] + (a)[5] * (b)[7] + (a)[8] * (b)[8];

void assign(double ubi[9], double gv[][3], double tol, double drlv2[],
            int labels[], int ig, int n);

void compute_gv(double xlylzl[][3], double omega[], double omegasign,
                double wvln, double wedge, double chi, double t[3],
                double gv[][3], int n);

void compute_xlylzl(double s[], double f[], double p[4], double r[9],
                    double dist[3], double xlylzl[][3], int n);
