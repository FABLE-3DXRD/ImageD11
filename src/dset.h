

/*
# ImageD11_v0.4 Software for beamline ID11
# Copyright (C) 2005  Jon Wright
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



#ifndef _dset_h
#define _dset_h
int * dset_initialise(int size); /* array to hold real values of each */
int * dset_new(int ** S);
void dset_makeunion(int * S, int r1, int r2);
void dset_link(int *S, int r1, int r2);
int dset_find(int x, int * S);

int *S;

int * dset_initialise(int size){
   int i;
   /* printf("Initialising the disjoint set at size %d, S=%x\n",size,S); */
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
      (int *) pS = S;
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
   if(r1>r2){S[r1]=r2;}
   if(r1<r2){S[r2]=r1;}
   /* if r1==r2 then they are already a union */
}

int dset_find(int x, int *S){
   if (x==0){
    /* oh shit */
    printf("Oh shit, you tried to find zero in your disjoint set, and it is not there!\n");
    return 0;
   }
   if (S[x] != x){
      S[x]=dset_find(S[x],S);
      }
   return S[x];
}

#endif /* _dset_h */

