
#include <math.h>
#include <stdio.h>


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
		float cut, float mean,
		int *s0, float *s1, float *s2){
  int i;
  float t;
  *s0 = npx;
  *s1 = 0.;
  *s2 = 0.;
  for( i = 0; i < npx; i++ ){
    if( img[i] > cut){      
      *s0 = *s0 - 1;
    } else {
      t = img[i] - mean;
      *s1 = *s1 + t;
      *s2 = *s2 + (t * t);
    }
  }
}


void threshold_image( float img[], int npx,
		      float nsigma,
		      char  msk[] ){
  float s1, s2, cut, mean, std;
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
  std  = sqrt(s2*s0 - s1*s1)/s0;
  cut  = mean + nsigma * std;
  printf("mean %f std %f cut %f %d %f %f\n", mean, std, cut,s0,s1,s2);
  for( i = 0; i < 5; i++ ){
    sigma_cut( img, npx, cut, mean, &s0, &s1, &s2);
    mean = s1/s0 + mean;
    std  = sqrt(s2*s0 - s1*s1)/s0;
    cut  = mean + nsigma * std;
    printf("mean %f std %f cut %f %d %f %f\n", mean, std, cut,s0,s1,s2);    
  }
  for( i = 0; i < npx; i++ ){
    if( img[i] > cut){
      msk[i] = 1;
    } else {
      msk[i] = 0;
    }
  }
}
		 
