



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
void dset_makeunion(int * S, int r1, int r2);
void dset_link(int *S, int r1, int r2);
int dset_find(int x, int * S);

int *S;

typedef struct {

  /* Spot_ID - to be generated when writing out */

  double s_1;       /* 1 Npix */
  double s_I;       /* 2 Sum intensity */
  double s_I2;      /* 3 Sum intensity^2 */
  double s_fI;      /* 4 Sum f * intensity */
  double s_ffI;     /* 5 Sum f * f* intensity */
  double s_sI;      /* 6 Sum s * intensity */
  double s_ssI;     /* 7 Sum s * s * intensity */
  double s_sfI;     /* 8 Sum f * s * intensity */

  double mx_I;      /* 9  Max intensity */
  double mx_I_f;    /* 10 fast at Max intensity */
  double mx_I_s;    /* 11 slow at Max intensity */

  double bb_mx_f;      /* 12 max of f */
  double bb_mx_s;      /* 13 max of s */
  double bb_mn_f;      /* 14 min of f */
  double bb_mn_s;      /* 15 min of s */ 
} blob2d;

#define NPROPERTY 15;

void new_blob(blob2d*, int , int , double);

void add_pixel(blob2d*, int , int , double);

void merge(blob2d*, blob2d *);

#endif /* _blobs_h */

