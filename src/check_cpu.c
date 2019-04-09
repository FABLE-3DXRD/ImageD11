
/** 
 * http://newbiz.github.io/cpp/2010/12/20/Playing-with-cpuid.html
 * ... with corrections and modifications...
 */



#include <stdio.h>
#include "check_cpu_auto.c"

int main(){
    int i;
    char name[13];  
    cpuidProcessorName( name );
    printf("%s\n",name);
    printf("SSE2 %d\n", i_have_SSE2());
    printf("SSE42 %d\n", i_have_SSE42());
    printf("AVX %d\n", i_have_AVX());
    printf("AVX2 %d\n", i_have_AVX2());
    printf("AVX512F %d\n", i_have_AVX512F());
    return 0;
}
