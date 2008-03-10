## Automatically adapted for numpy.oldnumeric Sep 06, 2007 by alter_code1.py


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

import numpy


class grain:
    def __init__(self,ubi,translation=None):
        self.ubi = numpy.array(ubi,numpy.float)
        self.ub = numpy.linalg.inv(ubi)
        if translation==None:
            self.translation = numpy.zeros(3, numpy.float)
        else:
            self.translation = numpy.array(translation,numpy.float)


    
def write_grain_file(filename, list_of_grains):
    f = open(filename, "w")
    for g in list_of_grains:
        t = g.translation
        f.write("#translation: %f %f %f\n"%(t[0],t[1],t[2]))
        if hasattr(g,"name"):
            f.write("#name %s\n"%(g.name))
        if hasattr(g,"x"):
            f.write("#npks %d\n"%(len(g.x)))
        f.write("#UBI:\n")
        u = g.ubi
        f.write("%f %f %f\n"  %(u[0,0],u[0,1],u[0,2]))
        f.write("%f %f %f\n"  %(u[1,0],u[1,1],u[1,2]))
        f.write("%f %f %f\n\n"%(u[2,0],u[2,1],u[2,2]))
    f.close()

def read_grain_file(filename):
    """read ubifile and return a list of ubi arrays """
    f = open(filename, "r")
    grainsread = []
    u = []
    t = [0,0,0]
    for line in f:
        if line.find("#translation:")==0:
            t = [ float(x) for x in line.split()[1:]]
            continue
        if line[0] == "#":
            continue
        vals = [ float(x) for x in line.split() ]
        if len(vals) == 3:
            u = u + [vals]
        if len(u)==3:
            grainsread.append( grain(u, t) )
            u = []
            t = [0,0,0]
    f.close()
    return grainsread
