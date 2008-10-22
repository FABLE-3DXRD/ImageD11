

#include <stdio.h>
#include "refs.h"

int main(int argc, char** argv){
    int i, nsym;    
    int hkl[] = {3,1,2,  -1,-2,-3,   -3,1,2} ;
    int res[3] ;
    int *kh, *kk, *kl;
    
    kh = &res[0];
    kk = &res[1];
    kl = &res[2];
    float sym[] = {
        1, 0, 0,
        0, 1, 0,
        0, 0, 1,
        0, 0, 0,

        0, 1, 0,
        0, 0, 1,
        1, 0, 0,
        0, 0, 0,

        0, 0, 1,
        1, 0, 0,
        0, 1, 0,
        0, 0, 0,
    };
    nsym = 3;
    
    for( i=0; i<3; i++){
        printf("%4d %4d %4d",hkl[i*3],hkl[i*3+1],hkl[i*3+2]);
        printf("->");
        *kh=0; *kk=0; *kl=0;
        get_key( hkl[i*3],hkl[i*3+1],hkl[i*3+2],
                sym, nsym,
                kh, kk, kl);
        printf("%4d %4d %4d\n",*kh,*kk,*kl);
    }
    return 0;
}
