


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
    
__constant float4 range4 = (float4)(0.0f, 1.0f, 2.0f, 3.0f);

__kernel void tthetaf4( __global float4 *tth,
                        __global float4 *eta,
                        __constant float* p )
{
    // data[i,j] => i == fast, j == slow
    int i;
    int j;
    float4 s,f,x,y,z;
    i = get_global_id(0);
    j = get_global_id(1);
    int address = i + j * get_global_size(0);

    f    = (4*i + range4  - p[FC]) * p[FS];
    s    = (j             - p[SC]) * p[SS];
    /*  |R0  R1  R2| |0|
        |R3  R4  R5|.|f|
        |R6  R7  R8| |s| */
    z  = p[R+7]*f + p[R+8]*s;
    y  = p[R+4]*f + p[R+5]*s;
    eta[address] = atan2pi( z, y );
    x  = p[R+1]*f + p[R+2]*s + p[DISTANCE];
    tth[address] = atan2pi( hypot( y, z) , x ); 
}


