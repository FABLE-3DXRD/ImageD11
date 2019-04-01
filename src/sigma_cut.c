
#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <string.h>


/**
 * Go through the data to compute mean and standard deviation and 
 * then repeat this masking out the values above c*std
 * 
 * @param img Input data array 
 * @param npx length of data array (contiguous)
 * @param cut cutoff for img < cut to be used
 * @param *s0  Number of pixels 
 * @param *s1  Sum of selected pixels 
 * @param *s2  Sum of selected pixel^2
 */
void sigma_cut( float img[], int npx,
		float cut, 
		int *s0, double *s1, double *s2){
  int i;
  double t;
  *s0 = npx;
  *s1 = 0.;
  *s2 = 0.;
  for( i = 0; i < npx; i++ ){
    if( img[i] > cut){      
      *s0 = *s0 - 1;
    } else {
      t = img[i];
      *s1 = *s1 + t;
      *s2 = *s2 + (t * t);
    }
  }
}


void threshold_image( float img[], int npx,
		      float nsigma,
		      char  msk[] ){
  double s1, s2, cut, mean, std;
  int i, s0;
  s0=npx;
  s1=0.;
  s2=0.;
  /* First guess is all data */
  for( i = 0; i < npx; i++ ){
      s1 = s1 + img[i];
      s2 = s2 + (img[i] * img[i]);
  }
  mean = s1/s0;
  std  = sqrt(s2/s0 - mean*mean);
  cut  = mean + nsigma * std;
  printf("mean %f std %f cut %f %d %f %f \n", mean, std, cut,s0,s1,s2);
  for( i = 0; i < 5; i++ ){
    cut  = mean + nsigma * std;
    sigma_cut( img, npx, cut, &s0, &s1, &s2);
    mean = s1/s0;
    std  = sqrt(s2*s0 - s1*s1)/s0;
    cut  = mean + nsigma * std;
    printf("mean %f std %f cut %f %f %d %f %f\n", mean, std, cut, cut,s0,s1,s2);    
  }
  for( i = 0; i < npx; i++ ){
    if( img[i] > cut){
      msk[i] = 1;
    } else {
      msk[i] = 0;
    }
  }
}
		 

void histogram_image( float img[],
		      int npx,
		      float low,
		      float step,
		      int32_t hist[],
		      int nhist){
  int i, ibin, s;
  float ostep, rbin, w;
  memset( hist, 0, nhist * sizeof( int32_t ) );
  ostep = 1./step;
  for( i=0; i<npx; i++){
    rbin = (img[i] - low)*ostep;
    if( rbin < (nhist-1) ){
      ibin = (int) rbin;
      hist[ibin]++; 
    }
  }
  /* 1: Locate median in the histogram */
  i=0;
  s=0;
  w=0;
  rbin = low;
  while( s < (npx/2) ){
    s = s + hist[i];
    w = w + hist[i]*rbin;
    i = i + 1;
    rbin = rbin + step;
  }
  /* 2: Compute mean and std of values less than the median
   *    based on the values from the histogram
   *
  mean = w/s; 
std = ???
  */
  
  
     
  
}
