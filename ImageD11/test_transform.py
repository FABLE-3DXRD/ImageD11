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

class testtransform(unittest.TestCase):
    def setUp(self):
        self.peaks = np.array([ [ 10, 20 , 50, 100 ],
                                [ 100, 10, 50, 99  ] ], np.float)

    def test_compute_tth_eta1(self):
        # Check translation of 0,0,0 has no effect
        yc = 49. ; ys=0.05 ; ty = 0.001
        zc = 51. ; zs=0.04 ; tz =-0.002
        dist = 112.345
        not_trans = transform.compute_tth_eta(self.peaks,
                                              y_center=yc, y_size=ys, tilt_y=ty,
                                              z_center=zc, z_size=zs, tilt_z=tz,
                                              distance=dist)
        om = np.ones(self.peaks.shape[1], np.float)
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


if __name__=="__main__":
    unittest.main()
