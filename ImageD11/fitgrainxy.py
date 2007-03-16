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

import Numeric

from ImageD11 import transform, indexing, closest

from ImageD11.grain import grain

class fitgrainxy:
    def __init__(self):
        self.parameters={'wedge':0.0 , 'chi':0.0}
        self.ubisread=[]
    
    def loadparameters(self,filename):
        lines = open(filename,"r").readlines()
        for line in lines:
            name,value=line.split()
            try:
                self.parameters[name]=float(value)
            except:
                self.parameters[name]=value

    def saveparameters(self,filename):
        out = open(filename,"w")
        keys=self.parameters.keys()
        keys.sort()
        for k in keys:
            try:
                out.write("%s %f\n"%(k,self.parameters[k]))
            except:
                out.write("%s %s\n"%(k,self.parameters[k]))
        out.close()
        
    def readubis(self,filename):
        """
        Save the generated ubi matrices into a text file
        """
        f=open(filename,"r")
        u = []
        i=0
        for line in f:
            vals = [ float(x) for x in line.split() ]
            if len(vals) == 3:
                u = u + [vals]
            if len(u)==3:
                # name=filename + " " + str(i)
                self.ubisread.append(Numeric.array(u))
                i=i+1
                u = []
        f.close()

    def loadfiltered(self,filename):
        """
        Read in file containing filtered and merged peak positions
        """
        f=open(filename,"r")
        #                   0123456789012
        line=f.readline()
        if line[0:12] !="# xc yc omega"[0:12]:
            print line
            raise Exception("Sorry That does not seem to be a filter peaks file, output from the peaksearching menu option")
        self.scantitles = line.replace("#","").split()
        bigarray=[]
        for line in f.readlines():
            v=[float(z) for z in line.split()]
            bigarray.append(v)
        f.close()
        self.scandata=Numeric.array(bigarray)
        self.allscandata=self.scandata.copy()

    def compute_gv(self,g):
        """
        Makes self.gv refer be g-vectors computed for this grain in this scan
        """

        try:
            xc      = self.scantitles.index("xc")
        except:
            print self.scantitles
            raise
#        print self.scandata[scanname].shape
        x  = self.scandata[:,xc]
        yc      = self.scantitles.index("yc")
        y  = self.scandata[:,yc]
        om      = self.scantitles.index("omega")
        om = self.scandata[:,om]

        tth,eta = transform.compute_tth_eta( Numeric.array([x, y]) ,
                                             self.parameters['y-center'],
                                             self.parameters['y-size'],
                                             self.parameters['tilt-y'],
                                             self.parameters['z-center'],
                                             self.parameters['z-size'],
                                             self.parameters['tilt-z'],
                                             self.parameters['distance']*1.0e3,
                                             crystal_translation = g.translation,
                                             omega = om,
                                             axis_orientation1 = self.parameters['wedge'],
                                             axis_orientation2 = self.parameters['chi'])
        self.gv = transform.compute_g_vectors(tth,eta,om,float(self.parameters['wavelength']), self.parameters['wedge'])
        self.gv = Numeric.transpose(self.gv)



    def map(self):
        for m in self.ubisread:
           g = grain(m)
           self.scandata=self.allscandata.copy()
           self.compute_gv(g)
           h=Numeric.matrixmultiply(g.ubi,Numeric.transpose(self.gv))
           hint=Numeric.floor(h+0.5).astype(Numeric.Int) # rounds down
           diff=h-hint
           drlv=Numeric.sqrt(Numeric.sum(diff*diff,0))
           indices = Numeric.compress(drlv < 0.05, range(self.scandata.shape[0]))
           print indices.shape,"hello"
           self.scandata = Numeric.take(self.allscandata,indices)
           npts = 10
           for i in range(npts):
               for j in range(npts):
                   x=(i-npts*0.5)*1000/npts
                   y=(j-npts*0.5)*1000/npts
                   g.translation[0]=x
                   g.translation[1]=y
                   self.compute_gv(g)
                   for tol in [ 0.1]:
                       mat = g.ubi.copy()
                       npks = closest.score_and_refine(mat, self.gv, tol)
                       # 2nd time with refined
                       npks = closest.score_and_refine(mat, self.gv, tol)
                       h=Numeric.matrixmultiply(mat,Numeric.transpose(self.gv))
                       hint=Numeric.floor(h+0.5).astype(Numeric.Int) # rounds down
                       diff=h-hint
                       drlv=Numeric.sqrt(Numeric.sum(diff*diff,0))
                       tthscore = Numeric.sum(Numeric.sum(hint*diff)*Numeric.sum(hint*diff)/Numeric.sum(h*h))
                       print x,y, tol, "%5d"%(npks),sum(drlv)/drlv.shape[0],tthscore/drlv.shape[0],indexing.ubitocellpars(mat),
                   print
               sys.stdout.flush()
           print
           print
           print

if __name__=="__main__":
    import sys
    o=fitgrainxy()
    o.loadparameters(sys.argv[1])
    o.readubis(sys.argv[2])
    o.loadfiltered(sys.argv[3])
    o.compute_gv(grain(o.ubisread[0]))
    o.map()
