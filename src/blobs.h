




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

#ifndef _blobs_h
#define _blobs_h
int * dset_initialise(int size); /* array to hold real values of each */
int * dset_new(int ** S);
void dset_makeunion(int *S, int r1, int r2);
void dset_link(int *S, int r1, int r2);
int dset_find(int x, int *S);


  /* Spot_ID - to be generated when writing out */
enum {
  s_1=0,       /* 1 Npix */
  s_I,       /* 2 Sum intensity */
  s_I2,      /* 3 Sum intensity^2 */
  s_fI,      /* 4 Sum f * intensity */
  s_ffI,     /* 5 Sum f * f* intensity */
  s_sI,      /* 6 Sum s * intensity */
  s_ssI,     /* 7 Sum s * s * intensity */
  s_sfI,     /* 8 Sum f * s * intensity */
  s_oI,         /* 9 sum omega * intensity */ 
  s_ooI,        /* 10 sum omega2 * intensity */ 
  s_soI,        /* 11 sum omega * s * intensity */
  s_foI,        /* 12 sum omega * f * intensity */

  mx_I,      /* 13  Max intensity */
  mx_I_f,    /* 14 fast at Max intensity */
  mx_I_s,    /* 15 slow at Max intensity */
  mx_I_o,    /* 16 omega at max I */

  bb_mx_f,      /* 17 max of f */
  bb_mx_s,      /* 18 max of s */
  bb_mx_o,      /* 19 max of omega */
  bb_mn_f,      /* 20 min of f */
  bb_mn_s,      /* 21 min of s */
  bb_mn_o,      /* 22 min of o */  


  avg_i,   /* Average intensity */   
  f_raw,       /* fast centre of mass */
  s_raw,       /* slow centre of mass */
  o_raw,       /* omega centre of mass */
  m_ss,   /* moments */
  m_ff,
  m_oo,
  m_sf,
  m_so,
  m_fo,

  f_cen,  /* Filled in elsewhere  - spatial correction */
  s_cen,  /* ditto */

  dety, /*Filled in elsewhere  - flip to HFPO book */
  detz, /*Filled in elsewhere  - flip to HFPO book */

  NPROPERTY     /* Number of properties if starting at 0 */
} ;





/*void new_blob(double blob[], int i, int j, double val);*/

void add_pixel(double blob[], int i, int j, double val, double omega);

void merge(double blob1[], double blob2[]);

void compute_moments(double blobs[], int nblobs);
#endif /* _blobs_h */

