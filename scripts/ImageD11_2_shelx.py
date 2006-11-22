#!/bliss/users/blissadm/python/bliss_python/suse82/bin/python


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



"""
Script for getting integrated intensities from .flt + .ubi file
from the command line
"""





if __name__=="__main__":
    import sys
    try:
        flt = sys.argv[1]
        ubi = sys.argv[2]
        par = sys.argv[3]
        out = sys.argv[4]
        from ImageD11.refinegrains import refinegrains
        o=refinegrains()
        o.loadparameters(par)
        o.loadfiltered(flt)
        o.readubis(ubi)
        o.generate_grains()
        o.tolerance = 0.1  # to command line arg
        o.refineubis(False)


    except:
        print "Usage: %s fltfile ubifile parfile outfile"
        raise
