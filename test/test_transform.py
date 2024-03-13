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


import unittest

import numpy as np

from ImageD11 import transform

#import sys
#sys.path.append(".")
#import transform

class testtransform(unittest.TestCase):
    def setUp(self):
        self.peaks = np.array([ [ 10, 20 , 50, 100, 512, 2000 ],
                                [ 100, 10, 50, 99, 518, 2010  ] ], float)

    def test_compute_tth_eta1(self):
        # Check translation of 0,0,0 has no effect
        yc = 49. ; ys=0.05 ; ty = 0.001
        zc = 51. ; zs=0.04 ; tz =-0.002
        dist = 112.345
        not_trans = transform.compute_tth_eta(self.peaks,
                                              y_center=yc, y_size=ys, tilt_y=ty,
                                              z_center=zc, z_size=zs, tilt_z=tz,
                                              distance=dist)
        om = np.ones(self.peaks.shape[1], float)
        trans =  transform.compute_tth_eta(self.peaks,
                                           y_center=yc, y_size=ys, tilt_y=ty,
                                           z_center=zc, z_size=zs, tilt_z=tz,
                                           distance=dist,
                                           t_x=0.,t_y=0.,t_z=0.,
                                           omega=om,
                                           wedge=10.,
                                           chi=-11.)

        diff = not_trans[0] - trans[0]
        self.assertAlmostEqual(np.sum(diff*diff), 0, 5 )
        diff = not_trans[1] - trans[1]
        self.assertAlmostEqual(np.sum(diff*diff), 0, 5 )



    def test_xyz_from_tth_eta(self):
        # Check the compute_xyz_from_tth_eta works
        # mm
        yc = 1032. ; ys=0.05  ; ty = 0.001
        zc = 1024. ; zs=0.047 ; tz =-0.002
        dist = 112.345
        omega = np.ones( len(self.peaks[0]) )
        tth, eta =  transform.compute_tth_eta(self.peaks,
                                              y_center=yc, y_size=ys, tilt_y=ty,
                                              z_center=zc, z_size=zs, tilt_z=tz,
                                              distance=dist)
        fc, sc = transform.compute_xyz_from_tth_eta(tth, eta, omega,
                                                    y_center=yc, y_size=ys, tilt_y=ty,
                                                    z_center=zc, z_size=zs, tilt_z=tz,
                                                    distance=dist)
        self.assertAlmostEqual( np.abs(sc - self.peaks[0]).sum(), 0, 5 )
        self.assertAlmostEqual( np.abs(fc - self.peaks[1]).sum(), 0, 5 )

    def test_xyz_from_tth_eta_trans(self):
        # Check the compute_xyz_from_tth_eta works
        # mm
        yc = 1032. ; ys=0.05  ; ty = 0.001
        zc = 1024. ; zs=0.047 ; tz =-0.002
        t_x = 0.1  ; t_y = 0.2 ; t_z = 0.3
        dist = 112.345
        omega = np.linspace(-720, 720, len(self.peaks[0]) )
        # With translation
        tth, eta =  transform.compute_tth_eta(self.peaks,
                                              omega= omega,
                                              t_x = t_x, t_y = t_y, t_z = t_z,
                                              y_center=yc, y_size=ys, tilt_y=ty,
                                              z_center=zc, z_size=zs, tilt_z=tz,
                                              distance=dist)
        fc, sc = transform.compute_xyz_from_tth_eta(tth, eta, omega,
                                                    t_x = t_x, t_y = t_y, t_z = t_z,
                                                    y_center=yc, y_size=ys, tilt_y=ty,
                                                    z_center=zc, z_size=zs, tilt_z=tz,
                                                    distance=dist)
        self.assertAlmostEqual( np.abs(sc - self.peaks[0]).sum(), 0, 5 )
        self.assertAlmostEqual( np.abs(fc - self.peaks[1]).sum(), 0, 5 )



if __name__=="__main__":
    unittest.main()
