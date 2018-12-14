
#include <stdint.h>
#include <stdio.h>

int lmw( const float* restrict im,
	 int32_t*  restrict l,
	 float threshold,
	 int dim0,
	 int dim1){
  float mx;
  int npk, iq, p, q, o[8]={ -1-dim1, -dim1, +1-dim1,
			    -1     ,          +1,
			    -1+dim1, +dim1, +1+dim1 };
  // Set edges to zero
  for( int i=0; i<dim0; i++){
    l[i]=0;
    l[dim1*i]=0;
    l[dim1*i+dim0-1]=0;
    l[dim1*dim0-i]=0;
  }
for( int i=1; i<dim0-1; i++){
    for( int j=1; j<dim1-1; j++){
      p = i*dim1+j;
      mx = im[p];
      if( mx < threshold ){
	l[p]=0;
	continue;
      }
      iq = p;
      for( int k=0; k<8; k++){
	q = p+o[k];
	if( im[q] > mx ){
	  mx = im[q];
	  iq = q;
	}
      }
      l[p]=iq;
    }
  } // i
  
  // Now make all point to their max
  int k, n, t;
  npk = 0;
  for( int i=1; i<dim0-1; i++){
    for( int j=1; j<dim1-1; j++){
      p = i*dim1+j;
      if( l[p] == 0 ) continue;
      if( l[p] == p ) { // on max
	npk += 1;
	continue;
      } else {
	k = l[p]; // next higher
	while( l[k] != k ){ // look for max
	  k = l[k]; // found it at k
	}
        n = l[p]; // now walk path renaming to k
	l[p] = k;
	while( l[n] != k ){
	  t = l[n];
	  l[n] = k;
	  n = t;
	}
      }
    }
  }
  // At this point we should renumber them with
  // new labels and/or sort by height, size etc
  return npk;
}

// gcc  -shared -fPIC -O3 -o lmw.so lmwold.c -std=c99
// python localmaxwshed.py



