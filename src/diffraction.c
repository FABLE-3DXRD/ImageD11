
/* C code to do some of the diffraction geometry computations
 * a bit quicker
 *
 * For ctypes on windows (Jon's office PC) this can compile as:
 *
 * MSVC:
 * cl diffraction.c /O2 /LD 
 * cl diffraction.c /Ox /fp:fast /favor:AMD64 /Qfast_transcendentals /LD
 * mt -manifest diffraction.dll.manifest -outputresource:diffraction.dll;2
 * 0.21 s
 *
 * or with gcc from c:\mingw64\bin\gcc:
 * gcc diffraction.c -O3 -march=native -ffast-math -static-libgcc -shared -o diffraction.dll
 *
 */

#include <math.h>
#include <float.h>
#include <stdio.h>


#ifdef __cplusplus
extern "C" {
#endif


typedef double real;
#define BIGNUM DBL_MAX;
  /*
typedef float real;
#define BIGNUM FLT_MAX;
  */
#define HALF 0.5
#define ZERO 0.0
#define ONE 1.0

typedef real rmatrix[9];
typedef real vector[3];


/* If you are reading this code you may find these defines
 * confusing. It is a shame to do it like this, I agree.
 * The expand inline to avoid the function call overhead
 * and allow the compiler the maximum opportunity to 
 * optimise the code as far as it is able
 */

#define norm3( a ) \
    ( sqrt ( (a)[0]*(a)[0] + (a)[1]*(a)[1] + (a)[2]*(a)[2] ))

#define dot3( a, b )\
    (( (a)[0]*(b)[0] + (a)[1]*(b)[1] + (a)[2]*(b)[2] ))

#define round3( a, b ) \
    (b)[0] = floor( a[0] + HALF); \
    (b)[1] = floor( a[1] + HALF); \
    (b)[2] = floor( a[2] + HALF); 

#define crossProduct(a,b,c) \
	(c)[0] = (a)[1] * (b)[2] - (b)[1] * (a)[2]; \
	(c)[1] = (a)[2] * (b)[0] - (b)[2] * (a)[0]; \
	(c)[2] = (a)[0] * (b)[1] - (b)[0] * (a)[1];

#define vec3copy(a,b) \
    (b)[0] = (a)[0]; \
    (b)[1] = (a)[1]; \
    (b)[2] = (a)[2];

#define vec3add(a,b,c) \
	(c)[0] = (a)[0] + (b)[0]; \
	(c)[1] = (a)[1] + (b)[1]; \
	(c)[2] = (a)[2] + (b)[2]; 

#define vec3sub(a,b,c) \
	(c)[0] = (a)[0] - (b)[0]; \
	(c)[1] = (a)[1] - (b)[1]; \
	(c)[2] = (a)[2] - (b)[2]; 

#define vec3smul(a,b,c) \
    (c)[0] = (a)[0]*(b);\
    (c)[1] = (a)[1]*(b);\
    (c)[2] = (a)[2]*(b);

#define vec3sdiv(a,b,c) \
    (c)[0] = (a)[0]/(b);\
    (c)[1] = (a)[1]/(b);\
    (c)[2] = (a)[2]/(b);

#define matVec(a,b,c) \
    (c)[0] = (a)[0]*(b)[0] + (a)[1]*(b)[1] + (a)[2]*(b)[2]; \
    (c)[1] = (a)[3]*(b)[0] + (a)[4]*(b)[1] + (a)[5]*(b)[2]; \
    (c)[2] = (a)[6]*(b)[0] + (a)[7]*(b)[1] + (a)[8]*(b)[2]; 

#define matTVec(a,b,c) \
    (c)[0] = (a)[0]*(b)[0] + (a)[3]*(b)[1] + (a)[6]*(b)[2]; \
    (c)[1] = (a)[1]*(b)[0] + (a)[4]*(b)[1] + (a)[7]*(b)[2]; \
    (c)[2] = (a)[2]*(b)[0] + (a)[5]*(b)[1] + (a)[8]*(b)[2]; 

#define matmat(a,b,c)\
    (c)[0] = a[0]*b[0] + a[1]*b[3] + a[2]*b[6]; \
    (c)[1] = a[0]*b[1] + a[1]*b[4] + a[2]*b[7]; \
    (c)[2] = a[0]*b[2] + a[1]*b[5] + a[2]*b[8]; \
    (c)[3] = a[3]*b[0] + a[4]*b[3] + a[5]*b[6]; \
    (c)[4] = a[3]*b[1] + a[4]*b[4] + a[5]*b[7]; \
    (c)[5] = a[3]*b[2] + a[4]*b[5] + a[5]*b[8]; \
    (c)[6] = a[6]*b[0] + a[7]*b[3] + a[8]*b[6]; \
    (c)[7] = a[6]*b[1] + a[7]*b[4] + a[8]*b[7]; \
    (c)[8] = a[6]*b[2] + a[7]*b[5] + a[8]*b[8];

#define matTmatT(a,b,c)\
    (c)[0] = a[0]*b[0] + a[3]*b[1] + a[6]*b[2]; \
    (c)[1] = a[0]*b[3] + a[3]*b[4] + a[6]*b[5]; \
    (c)[2] = a[0]*b[6] + a[3]*b[7] + a[6]*b[8]; \
    (c)[3] = a[1]*b[0] + a[4]*b[1] + a[7]*b[2]; \
    (c)[4] = a[1]*b[3] + a[4]*b[4] + a[7]*b[5]; \
    (c)[5] = a[1]*b[6] + a[4]*b[7] + a[7]*b[8]; \
    (c)[6] = a[2]*b[0] + a[5]*b[1] + a[8]*b[2]; \
    (c)[7] = a[2]*b[3] + a[5]*b[4] + a[8]*b[5]; \
    (c)[8] = a[2]*b[6] + a[5]*b[7] + a[8]*b[8];

#define matTmat(a,b,c)\
    (c)[0] = a[0]*b[0] + a[3]*b[3] + a[6]*b[6]; \
    (c)[1] = a[0]*b[1] + a[3]*b[4] + a[6]*b[7]; \
    (c)[2] = a[0]*b[2] + a[3]*b[5] + a[6]*b[8]; \
    (c)[3] = a[1]*b[0] + a[4]*b[3] + a[7]*b[6]; \
    (c)[4] = a[1]*b[1] + a[4]*b[4] + a[7]*b[7]; \
    (c)[5] = a[1]*b[2] + a[4]*b[5] + a[7]*b[8]; \
    (c)[6] = a[2]*b[0] + a[5]*b[3] + a[8]*b[6]; \
    (c)[7] = a[2]*b[1] + a[5]*b[4] + a[8]*b[7]; \
    (c)[8] = a[2]*b[2] + a[5]*b[5] + a[8]*b[8];

#define matmatT(a,b,c)\
    (c)[0] = a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; \
    (c)[1] = a[0]*b[3] + a[1]*b[4] + a[2]*b[5]; \
    (c)[2] = a[0]*b[6] + a[1]*b[7] + a[2]*b[8]; \
    (c)[3] = a[3]*b[0] + a[4]*b[1] + a[5]*b[2]; \
    (c)[4] = a[3]*b[3] + a[4]*b[4] + a[5]*b[5]; \
    (c)[5] = a[3]*b[6] + a[4]*b[7] + a[5]*b[8]; \
    (c)[6] = a[6]*b[0] + a[7]*b[1] + a[8]*b[2]; \
    (c)[7] = a[6]*b[3] + a[7]*b[4] + a[8]*b[5]; \
    (c)[8] = a[6]*b[6] + a[7]*b[7] + a[8]*b[8];

#define fillomegamatrix(angle, M) \
    (M)[0] = cos( angle ); \
    (M)[1] = sin( angle ); \
    (M)[2] = ZERO; \
    (M)[3] = -sin( angle ); \
    (M)[4] = cos( angle ); \
    (M)[5] = ZERO; \
    (M)[6] = ZERO; \
    (M)[7] = ZERO; \
    (M)[8] = ONE;  

#define fillwedgematrix(angle, M) \
    (M)[0] = cos( angle ); \
    (M)[1] = ZERO; \
    (M)[2] = sin( angle ); \
    (M)[3] = ZERO; \
    (M)[4] = ONE; \
    (M)[5] = ZERO; \
    (M)[6] = -sin( angle ); \
    (M)[7] = ZERO; \
    (M)[8] = cos( angle );

#define fillchimatrix(angle, M) \
    (M)[0] = ONE; \
    (M)[1] = ZERO; \
    (M)[2] = ZERO; \
    (M)[3] = ZERO; \
    (M)[4] = cos( angle ); \
    (M)[5] = -sin( angle ); \
    (M)[6] = ZERO; \
    (M)[7] = sin( angle ); \
    (M)[8] = cos( angle );  


#ifdef DEBUG
#define pvec(s,v)\
    printf("%s %f %f %f\n",(s),(v)[0],(v)[1],(v)[2])

#define pmat(s,v)\
    printf("%s %f %f %f %f %f %f %f %f %f\n",(s),(v)[0],(v)[1],(v)[2],\
            (v)[3],(v)[4],(v)[5],(v)[6],(v)[7],(v)[8])
#define psca(s,v)\
    printf("%s %f\n",(s),(v))
#else
#define pvec(s,v)
#define pmat(s,v)
#endif 



#ifdef _WIN32
/* This line is quite upsetting...: */
__declspec(dllexport)
#endif
void c_choosegrains( vector XL[], real omega[], 
        real wedge, real chi, real wvln,
        rmatrix UB[], rmatrix UBI[], vector T[], 
        int ngrains, int npks,
        int labels[], vector hkls[]){
    int i, j, gbest;
    real MdL, scal, mod_dG, scor, scormin;
    rmatrix omat, cmat, wmat, tmpmat, gomat, diffrot, wc;
    vector origin, dL, s, beam = {ONE,ZERO,ZERO}, 
           kvec, gvec, hkl, ih, gcalc, dG, 
           hbest = {ZERO,ZERO,ZERO};




    fillwedgematrix( wedge, wmat )
    fillchimatrix( chi,     cmat )

    psca( "wedge",wedge); 
    pmat( "wmat", wmat );
    psca( "chi",chi); 
    pmat( "cmat", cmat );
    matTmatT( wmat, cmat, wc);
    pmat( "wc", wc);
    vec3smul(beam, ONE/wvln, beam);
    for(i=0; i<npks; i++){

        /* gvecs, in crystal:
         * g = (Omega)(Chi)(Wedge) k
         * So grain origins are at (Transpose is inverse)
         * O = (WedgeT)(ChiT)(OmegaT)(Translation)
         * O = gomat.(Translation)
         */
        fillomegamatrix( omega[i], omat );

        pmat( "omat", omat );

	/* wedgeT.chiT.omegaT */
        matmatT( wc, omat, gomat); /* Grain origin matrix */
        pmat( "gomat", gomat );

        matmat( omat, cmat, tmpmat ); /* omega.chi.wedge */
        matmat( tmpmat, wmat, diffrot);

        scormin = BIGNUM; /* Tolerance */
        gbest = -1;


        for(j=0; j<ngrains; j++){
            matVec( gomat, T[j], origin ); /* This grain origin */
            vec3sub( XL[i], origin, dL);   /* Ray vector in laboratory */
            MdL = norm3( dL );            /* Length of that vector */
            pvec( "T[j]",T[j]);
            pvec( "XL[i]",XL[i]);
            pvec( "origin",origin);
            pvec( "dL", dL );
	    pvec( "beam",beam);

            scal = ONE /( MdL * wvln );    /* Go to reciprocal space */
            vec3smul( dL, scal , s);
            vec3sub( s, beam, kvec);       /* Vector in lab */
            matVec( diffrot, kvec, gvec);  /* Vector in xtal */
            matVec(UBI[j], gvec, hkl );       /* hkl indices */
            pvec("hklfloat", hkl );
            round3(hkl, ih);               /* Nearest integers */

            pvec("s", s);
            pvec("trans", T[j] );
            pvec("kvec", kvec );
            pvec("gvec", gvec );
            pvec("hkl", hkl );
            pmat("UBI", UBI[j] );




            /* Error analysis */
            matVec(UB[j], ih, gcalc);         
            vec3sub( gvec, gcalc, dG );
            scor = mod_dG = norm3( dG );
            /* There is a two theta and angular part to the error etc
             * Ideally the user should configure the score function?
             */
            if (scor < scormin ){
                scormin = scor;
                gbest = j;
                vec3copy(ih, hbest);
            }
        }
	/* printf("i %d, gbest %d",i,gbest); */
        labels[i] = gbest;
        vec3copy( hbest, hkls[i] );
    }
}






#ifdef __cplusplus
}
#endif


