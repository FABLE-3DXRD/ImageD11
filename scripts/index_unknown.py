#!/usr/bin/env python

from __future__ import print_function


# ImageD11_v1 Software for beamline ID11
# Copyright (C) 2008  Jon Wright
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  0211-1307  USA

import os, sys, logging
from argparse import ArgumentParser
from ImageD11 import ImageD11options
from ImageD11 import fft_index_refac, rc_array,\
    lattice_reduction, indexing
import numpy as np
from ImageD11.fft_index_refac import grid


def get_options(parser):
    parser.add_argument('-g', '--gve',
                      action = 'store',
                      type=ImageD11options.GvectorFileType(mode='r'),
                      dest = 'gvefilename',
                      default = None,
                      help = "Filename for g-vectors")
    parser.add_argument('-k', '--ngrains',
                      action = 'store',
                      dest = 'ngrains',
                      type = int,
                      default = 1,
                      help = "number of grains to try to find")
    parser.add_argument('-o', '--output',
                      action = 'store',
                      type=ImageD11options.UbiFileType(mode='w'),
                      default = 'grains.ubi',
                      dest = 'outfile',
                      help = "Name of ubi file to save grains in")

    parser = lattice_reduction.get_options(parser)

    parser.add_argument('--fft',
                      action = 'store_true',
                      dest = 'use_fft',
                      default = False,
                      help = "Use fft to generate lattice vectors [False]")

    parser.add_argument('--score_fft',
                      action = 'store_true',
                      dest = 'score_fft',
                      default = False,
                      help = "Score fft peaks using fft peaks first [True]")

    parser.add_argument('--no_sort',
                      action = 'store_false',
                      dest = 'sort_gve',
                      default = True,
                      help = "Sorting the gvector by length before indexing [True]")

    parser.add_argument('--noisy',
                      action = 'store_true',
                      dest = 'noisy',
                      default = False,
                      help = "Print more output")
    fft_index_refac.get_options( parser )
    return parser






if __name__=="__main__":

    logging.basicConfig( level=logging.INFO )
    parser = get_options( ArgumentParser() )

    options = parser.parse_args()
    
    
    o = indexing.indexer()
    if options.gvefilename is None or \
            not os.path.exists(options.gvefilename):
        print("You need to supply a gvector file with the -g option")
        sys.exit()

    o.readgvfile(options.gvefilename)


    mags = [ [np.dot( g, g ), tuple(g)] for g in o.gv]
    mags.sort()
    sorted = np.array( [ np.array(g[1]) for g in mags ] )
    
    all_gvecs = rc_array.rc_array(sorted.copy(), direction = 'row' )
    cur_gvecs = rc_array.rc_array(sorted.copy(), direction = 'row' )

    ubis = []
    
    try:
        for i in range(options.ngrains):
            if len(cur_gvecs) < 3:
                print("Ran out of unindexed peaks")
                break

            # print "Peak remaining",len(cur_gvecs)
            if options.use_fft:
                # do fft
                g = grid( npx = options.npx,
                          mr = options.mr,
                          nsig = options.nsig)
                g.gv_to_grid_new(cur_gvecs)
                g.fft()
                g.props()
                g.peaksearch(open("eu.patterson_pks","w"))
                g.read_peaks("eu.patterson_pks")
                vecs = rc_array.rc_array(g.UBIALL.T , direction='col')
                assert vecs.shape == (3, len(g.UBIALL))
                order = np.argsort( g.colfile.sum_intensity )[::-1]
                vecs = rc_array.rc_array( np.take( vecs, order, axis = 1),
                                          direction = 'col')
                assert g.gv.shape[1] == 3
                print("Finding lattice from patterson")
                # Go through and make a shortlist of lattices by scoring fft only
                if options.score_fft:
                    test_set = vecs
                else:
                    test_set = all_gvecs
                    #for v in vecs[:options.n_try]:
                    #     print v,np.sqrt(np.dot(v,v))
                try:
                    # mv2 = 1/options.min_vec2
                    l = lattice_reduction.find_lattice(
                        vecs,
                        min_vec2 = options.min_vec2,
                        n_try = options.n_try,
                        test_vecs = test_set,
                        tol = options.tol,
                        fraction_indexed = options.fraction_indexed,
                        noisy = options.noisy,
                        )
                except IndexError:
                    print(vecs, options.n_try)
                    print(vecs.shape)
                    raise
            else:
                # Test vectors are the g-vectors
                # Sort by length before searching??
                logging.info("Doing lattice search")
                l = lattice_reduction.find_lattice(
                    cur_gvecs, 
                    min_vec2=options.min_vec2,
                    n_try=options.n_try,
                    test_vecs = all_gvecs, 
                    tol = options.tol,
                    fraction_indexed = options.fraction_indexed,
                    noisy = options.noisy,
                    )
            if l is None:
                break
            l.r2c = indexing.refine( l.r2c, all_gvecs, tol = options.tol)
            l.c2r = np.linalg.inv(l.r2c)

            #print "all_gv 3,4", all_gvecs[3:5]
            #print np.dot(l.r2c, all_gvecs[3:5].T)
            #print l.score(all_gvecs[3:5], tol = options.tol)
            #print l.remainders( all_gvecs[3:5])
                        
            all = len(all_gvecs)
            npks = l.score( all_gvecs , tol = options.tol )
            if (1.0*npks/len(all_gvecs)) < options.fraction_indexed:
                break
            dr = indexing.calc_drlv2( l.r2c, all_gvecs) 
            t2 = options.tol * options.tol
            # print "tol, dr",t2,dr.shape, [d for d in dr[:10]] # print it !!!
            n2 = np.sum(np.where(dr < t2 , 1, 0 ))
            n3 = np.sum(np.where(dr > t2 , 1, 0 ) )
            assert npks == n2 , "debug scoring"
            assert n3 == len(all_gvecs) - npks, "debug scoring"
            print(l.r2c)
            print("Unit cell:",(6*"%.6f ")%indexing.ubitocellpars(l.r2c))
            
            # Put this grain in the output list
            ubis.append(l.r2c)

            # Remove from the gvectors
            drlv2 = indexing.calc_drlv2( l.r2c, cur_gvecs )
            # print drlv2[:20]
            # print cur_gvecs.shape
            # print drlv2.shape,drlv2[:10],options.tol*options.tol
            
            cur_gvecs = rc_array.rc_array(
                np.compress( drlv2 > options.tol*options.tol,
                             cur_gvecs, axis = 0),
                direction = 'row')
            print("Lattice found, indexes", npks, "from all",all)
            print("Number of unindexed peaks remaining %d"%(len(cur_gvecs)))
            print("Current vector shape",cur_gvecs.shape)
        if len(ubis)>0:
            indexing.write_ubi_file(options.outfile,ubis)
            print("Wrote to file",options.outfile)
        else:
            print("No unit cell found, sorry, please try again")
    except:
        if len(ubis)>0:
            indexing.write_ubi_file(options.outfile,ubis)
            print("Wrote to file",options.outfile)
        raise 
