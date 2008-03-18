

/*
# ImageD11_v1.0 Software for beamline ID11
# Copyright (C) 2005-2007  Jon Wright
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include "blobs.h"
#include <stdlib.h> /* malloc */
#include <stdio.h>
#include <math.h> /* sqrt */

void compute_moments( double b[], int nb  ){
  /* Convert the blob (summing) representation to a human readable form */
  int i, off;
  double ts,to,tf,us,uo,uf,vs,vo,vf, tc;

  /* loop over a set of blobs in a flat array */

  for(i=0; i < nb ; i++){
    /* avg_i */
    
    off = i*NPROPERTY;

    if( b[s_1 + off] == 0){
      /* Zero pixels - skip this fella! */
      /* printf("Skipped %d\n",i);*/
      continue;
    }

    tc = b[ s_I + off ]; /* tc = total counts */

    b[avg_i + off] =  b[ s_I  + off ]/b[s_1 + off]; 

    /* First and second moments */


    uf =  b[ s_fI  + off ]/tc;
    tf =  b[ s_ffI + off ]/tc;
    b[f_raw + off] = uf;
    vf = tf - uf*uf;
    if (( vf + 1.0) > 0.0) b[m_ff + off] = sqrt(vf+1.);
    else b[m_ff + off] = 1.0;

    
    us =  b[ s_sI  + off ]/tc;
    ts =  b[ s_ssI + off ]/tc;
    b[s_raw + off] = us;
    vs = ts - us*us;
    if (( vs + 1.0) > 0.0) b[m_ss + off] = sqrt(vs+1.);
    else b[m_ss + off] = 1.0;

    uo =  b[ s_oI  + off ]/tc;
    to =  b[ s_ooI + off ]/tc;
    b[o_raw + off] = uo;
    vo = to - uo*uo;
    if (( vo + 1.0) > 0.0) b[m_oo + off] = sqrt(vo+1.);
    else b[m_oo + off] = 1.;

    b[m_so + off] = (b[s_soI + off] / tc - us*uo )/b[m_ss + off]/b[m_oo + off] ;

    b[m_fo + off] = (b[s_foI + off] / tc - uf*uo )/b[m_ff + off]/b[m_oo + off] ;

    b[m_sf + off] = (b[s_sfI + off] / tc - us*uf )/b[m_ss + off]/b[m_ff + off] ;

  }

  
}


void add_pixel( double b[], int s, int f, double I, double o){


  b[s_1]   += 1;       /* Npix*/
  b[s_I]   += I;       /* Sum intensity */
  b[s_I2]  += I*I;     /* Sum intensity^2 */
  b[s_fI]  += f*I;     /* Sum f * intensity */
  b[s_ffI] += f*f*I;   /* Sum f * f* intensity */
  b[s_sI]  += s*I;     /* Sum s * intensity */
  b[s_ssI] += s*s*I;   /* Sum s * s * intensity */
  b[s_sfI] += s*f*I;   /* Sum f * s * intensity */
  b[s_oI]  += o*I;     /* Sum o * intensity */
  b[s_ooI] += o*o*I;   /* Sum o * intensity */
  b[s_soI] += s*o*I;   /* Sum s * o * intensity */
  b[s_foI] += f*o*I;   /* Sum f * o * intensity */



  if(I > b[mx_I]){
    b[mx_I]   = I;     /* Max intensity */
    b[mx_I_f] = f;     /* fast at Max intensity */
    b[mx_I_s] = s;     /* slow at Max intensity */
    b[mx_I_o] = o;
  }

  
  /* Bounding box */
  b[bb_mx_f] = ((f > b[bb_mx_f])? f :  b[bb_mx_f] );  
  b[bb_mx_s] = ((s > b[bb_mx_s])? s :  b[bb_mx_s] );  
  b[bb_mx_o] = ((o > b[bb_mx_o])? o :  b[bb_mx_o] );  
  b[bb_mn_f] = ((f < b[bb_mn_f])? f :  b[bb_mn_f] );  
  b[bb_mn_s] = ((s < b[bb_mn_s])? s :  b[bb_mn_s] );  
  b[bb_mn_o] = ((o < b[bb_mn_o])? o :  b[bb_mn_o] );  

}


void merge(double b1[], double b2[]){
  /* b2 is killed, b1 is kept */ 
  int i;

  b1[s_1]   += b2[s_1];         /* Npix*/
  b1[s_I]   += b2[s_I];         /* Sum intensity */
  b1[s_I2]  += b2[s_I2];        /* Sum intensity^2 */
  b1[s_fI]  += b2[s_fI];        /* Sum f * intensity */
  b1[s_ffI] += b2[s_ffI];       /* Sum f * f* intensity */
  b1[s_sI]  += b2[s_sI];        /* Sum s * intensity */
  b1[s_ssI] += b2[s_ssI];       /* Sum s * s * intensity */
  b1[s_sfI] += b2[s_sfI];       /* Sum f * s * intensity */
  b1[s_oI]  += b2[s_oI];        /* Sum o * intensity */
  b1[s_ooI] += b2[s_ooI];       /* Sum o * intensity */
  b1[s_soI] += b2[s_soI];       /* Sum s * o * intensity */
  b1[s_foI] += b2[s_foI];       /* Sum f * o * intensity */
  
  if(b2[mx_I] > b1[mx_I]){
    b1[mx_I]   = b2[mx_I];      /* Max intensity */
    b1[mx_I_f] = b2[mx_I_f];    /* fast at Max intensity */
    b1[mx_I_s] = b2[mx_I_s];    /* slow at Max intensity */
    b1[mx_I_o] = b2[mx_I_o];    /* slow at Max intensity */
  }

  /* Bounding box */
  b1[bb_mx_f] = ((b2[bb_mx_f] > b1[bb_mx_f])? b2[bb_mx_f] :  b1[bb_mx_f] );  
  b1[bb_mx_s] = ((b2[bb_mx_s] > b1[bb_mx_s])? b2[bb_mx_s] :  b1[bb_mx_s] );  
  b1[bb_mx_o] = ((b2[bb_mx_o] > b1[bb_mx_o])? b2[bb_mx_o] :  b1[bb_mx_o] );  
  b1[bb_mn_o] = ((b2[bb_mn_o] < b1[bb_mn_o])? b2[bb_mn_o] :  b1[bb_mn_o] );  
  b1[bb_mn_f] = ((b2[bb_mn_f] < b1[bb_mn_f])? b2[bb_mn_f] :  b1[bb_mn_f] );  
  b1[bb_mn_s] = ((b2[bb_mn_s] < b1[bb_mn_s])? b2[bb_mn_s] :  b1[bb_mn_s] );  

  /* Trash b2 to be on the safe side */
  for(i=0;i<NPROPERTY;i++)b2[i]=0;

}


/* ********************************************************************
 * Disjoint set
 *
 * See a good book on algorithms (eg)
 *    "Introduction to Algorithms"
 *    Cormen, Leiserson, Rivest, Stein
 *    MIT Press
 *
 * Stored in array S
 *
 * S[0] = malloc'ed length of the array
 * S[:-1] = current end element 
 *
 * dset_initialise(int size) => S[:] = 0
 * dset_new(S) => new item at end of working part of array
 * dset_makeunion(S,i,j) => make S[i]==S[j]
 * dset_find(i,S) => return set of i
 *
 * eg: sets are 1,2,3,4,5
 *     Now say 1 and 3 are actually the same thing:
 *     1,2,3=1,4,5
 *     Ask for set of 3 => S[3]=1 => S[1]=1 => return 1
 *
 * ********************************************************************  */


int * dset_initialise(int size){
   int i, *S;
   S=(int *) ( malloc(size*sizeof(int) ) );
   if (S==NULL){
   printf("Memory allocation error in dset_initialise\n");
   exit(1);
   }
   for(i=0;i<size;i++)S[i]=0;
   S[0]=size;
   return S;
}



int * dset_new( int ** pS){
   /* S[0] will always hold the length of the array */
   /* S[:-1] holds the current element */
   int length, current, i;
   int *S;
   S = (int *) * pS;
   length=S[0];
   current=(++S[S[0]-1]);
   if(current+3>length){
      S=(int *)(realloc(S,length*2*sizeof(int))); 
      if(S==NULL){
      printf("Memory allocation error in dset_new\n");
      exit(1);
      }
      /*      printf("Realloced S to %d in dset_new\n",length*2); */
      /* Fails on gcc 4 but worked on gcc 3.x
	 dont actually what it does or why it is there
	 hence commented it out
      (int *) pS = S; */
      S[0]=length*2;
      S[length-1]=0;
      for(i=length-1;i<length*2;i++)S[i]=0;
      S[length*2-1]=current;
   }
   S[current]=current;
   return S;
}



void dset_makeunion(int * S, int r1, int r2){
   int a,b;
   a=dset_find(r1,S);
   b=dset_find(r2,S);
   dset_link(S, a, b);
}



void dset_link(int * S, int r2, int r1){
   if (r1 > r2){
     /* The higher # r1 is changed to point to the lower */
     S[r1] = r2;
   } 
   if(r1 < r2){
     /* Again, the higher number, r2 is the one to point to lower */
     S[r2] = r1;
   } 
   /* if r1==r2 then they are already a union */
}



int dset_find(int x, int *S){
   if (x==0){
    /* oh shit */
    printf("Oh dear, you tried to find zero in your "\
	   "disjoint set, and it is not there!\n");
    return 0;
   }
   if (S[x] != x){
      S[x]=dset_find(S[x],S);
      }
   return S[x];
}
