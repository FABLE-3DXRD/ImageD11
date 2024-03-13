


/*
 * Center in fast, slow
 * size in fast, slow 
 * distance
 * tilts
 */


/* Positions in parameter array */

#define FC 0
#define FS 1
#define SC 2
#define SS 3
#define DISTANCE 4
#define R 5
    

__kernel void tthetaf( __global float  *tth,
                       __global float  *eta,
                       __constant float* p )
{
    // data[i,j] => i == fast, j == slow
    int i;
    int j;
    float s,f,r,x,y,z;
    j = get_global_id(0);
    i = get_global_id(1);
    int address =  j * get_global_size(1) + i;

    f    = (i - p[FC]) * p[FS];
    s    = (j - p[SC]) * p[SS];
    /*  |R0  R1  R2| |0|
        |R3  R4  R5|.|f|
        |R6  R7  R8| |s| */
    z  = p[R+7]*f + p[R+8]*s;
    y  = p[R+4]*f + p[R+5]*s;
    eta[address] = atan2pi( z, y );
    x  = p[R+1]*f + p[R+2]*s + p[DISTANCE];
    r  = sqrt( y*y + z*z );
    tth[address] = atan2pi( r , x ); 
}


