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

import numpy, math
import ImageD11.indexing


class grain:
    def __init__(self,ubi,translation=None, **kwds):
        self.ubi = numpy.array(ubi,numpy.float)
        self.ub = numpy.linalg.inv(ubi)
        self.u = ImageD11.indexing.ubitoU(self.ubi)
        try:
            self.Rod = ImageD11.indexing.ubitoRod(self.ubi)
        except:
            print self.ubi
            raise
        self.mt = numpy.dot(self.ubi, self.ubi.T)
        self.rmt = numpy.linalg.inv(self.mt)
        if translation is None:
            # If translation has not been read from ubi file make it 
            # be None to avoid confusion with a grain which is known
            # to be at [0,0,0]
            self.translation = None
        else:
            self.translation = numpy.array(translation,numpy.float)
    def set_ubi(self, ubi):
        self.ubi = numpy.array(ubi,numpy.float)
        self.ub = numpy.linalg.inv(ubi)
        self.u = ImageD11.indexing.ubitoU(self.ubi)
        self.Rod = ImageD11.indexing.ubitoRod(self.ubi)
        self.mt = numpy.dot(self.ubi, self.ubi.T)
        self.rmt = numpy.linalg.inv(self.mt)
        


    
def write_grain_file(filename, list_of_grains):
    f = open(filename, "w")
    for g in list_of_grains:
        t = g.translation
        f.write("#translation: %g %g %g\n"%(t[0],t[1],t[2]))
        if hasattr(g,"name"):
            f.write("#name %s\n"%(g.name.rstrip()))
        if hasattr(g,"intensity_info"):
            f.write("#intensity_info %s\n"%(g.intensity_info.rstrip()))
        if hasattr(g,"npks"):
            f.write("#npks %d\n"%(int(g.npks)))
        if hasattr(g,"Rod"):
            try:
                f.write("#Rod %f %f %f\n"%tuple([float(r) for r in g.Rod]))
            except:
                f.write("#Rod %s"%(g.Rod))
        f.write("#UBI:\n")
        u = g.ubi
        # More than float32 precision
        f.write("%.9g %.9g %.9g\n"  %(u[0,0],u[0,1],u[0,2]))
        f.write("%.9g %.9g %.9g\n"  %(u[1,0],u[1,1],u[1,2]))
        f.write("%.9g %.9g %.9g\n\n"%(u[2,0],u[2,1],u[2,2]))
    f.close()

def read_grain_file(filename):
    """read ubifile and return a list of ubi arrays """
    f = open(filename, "r")
    grainsread = []
    u = []
    t = None
    p = {}
    for line in f:
        if line.find("#translation:")==0:
            t = [ float(x) for x in line.split()[1:]]
            continue
        if line[0] == "#" and line.find("UBI")<0:
            k,v=line[1:].split(" ",1)
            p[k]=v
            continue
        if line[0] == "#" and line.find("intensity_info")>-1:
            p["intensity_info"] = line.split("intensity_info")[1].rstrip()
        if line.find("#")==0: continue
        vals = [ float(x) for x in line.split() ]
        if len(vals) == 3:
            u = u + [vals]
        if len(u)==3:
            grainsread.append( grain(u, t) )
            for k in ["name","npks","Rod","intensity_info"]:
                if p.has_key(k):
                    setattr(grainsread[-1], k, p[k])
            p={}
            u = []
            t = None
    f.close()
    return grainsread
