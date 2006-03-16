
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

import transform
import indexing
from grain import grain

class refinegrains:
    def __init__(self,tolerance=0.01):
        self.tolerance=tolerance
        # list of ubi matrices (1 for each grain in each scan)
        self.grainnames = []
        self.ubisread = {}
        # list crystal translations (1 for each matrix in each scans)
        # list of scans and corresponding data
        self.scannames=[]
        self.scantitles={}
        self.scandata={}
        # grains in each scan
        self.grains = {}
        self.wedge = 0.0
        self.chi = 0.0
        self.parameters = {}
        self.varyingtranslations=False

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
                out.write("%s %f\n"%(key,self.parameters[key]))
            except:
                out.write("%s %s\n"%(key,self.parameters[key]))
        f.close()
                    
        
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
                name=filename + " " + str(i)
                self.grainnames.append(name)
                self.ubisread[name]=Numeric.array(u)
                i=i+1
                u = []
        f.close()
        #for key in self.ubisread.keys():
        #   print self.ubisread[key]


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
        titles = line.replace("#","").split()
        bigarray=[]
        for line in f.readlines():
            v=[float(z) for z in line.split()]
            bigarray.append(v)
        f.close()
        self.scannames.append(filename)
        self.scantitles[filename]=titles
        self.scandata[filename]=Numeric.array(bigarray)
        
                         
    def generate_grains(self):
        for grainname in self.grainnames:
            for scanname in self.scannames:
                try:
                    gr = self.grains[(grainname,scanname)]
                except KeyError:
                    self.grains[(grainname,scanname)] = grain(self.ubisread[grainname])
        #for key in self.grains.keys():
        #    print key,":",self.grains[key].ubi

    def compute_gv(self,grainname,scanname):
        """
        Makes self.gv refer be g-vectors computed for this grain in this scan
        """
        try:
            xc      = self.scantitles[scanname].index("xc")
        except:
            print self.scantitles[scanname]
#        print self.scandata[scanname].shape
        x  = self.scandata[scanname][:,xc]
        yc      = self.scantitles[scanname].index("yc")
        y  = self.scandata[scanname][:,yc]
        om      = self.scantitles[scanname].index("omega")
        om = self.scandata[scanname][:,om]
        g = self.grains[(grainname,scanname)]

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
        self.gv = transform.compute_g_vectors(tth,eta,om,float(self.parameters['wavelength']), self.wedge)
        self.gv = Numeric.transpose(self.gv)


    def refine(self,ubi,quiet=True):
        from ImageD11 import closest
        mat=ubi.copy()
        npks = closest.score_and_refine(mat, self.gv, self.tolerance)
        if not quiet:
            print npks 
        #mat = indexing.refine(ubi,self.gv,self.tolerance,quiet=quiet)
        #print ubi, testmatrix,ubi-testmatrix,mat-testmatrix
        return mat
    
    def gof(self,args):
        """
        <drlv> for all of the grains in all of the scans
        """

        self.applyargs(args)
        diffs = 0.
        contribs = 0.
        sumdrlv = 0.
        goodpks = 0
        
        
        first = True
        for key in self.grains.keys():
            g = self.grains[key]
            grainname = key[0]
            scanname = key[1]
            if self.varyingtranslations or first:
                first = False
                # Compute gv using current parameters
                self.compute_gv(grainname,scanname)
                #print self.gv.shape
                #print self.gv[0:10,:]

            g.ubi = self.refine(g.ubi)

            h=Numeric.matrixmultiply(g.ubi,Numeric.transpose(self.gv))
            hint=Numeric.floor(h+0.5).astype(Numeric.Int) # rounds down
            diff=h-hint
            drlv=Numeric.sqrt(Numeric.sum(diff*diff,0))
            t=self.tolerance
            n_ind = Numeric.sum(Numeric.where(drlv<t,1,0)) # drlv
            drlv = Numeric.where(drlv<t,drlv,0)
            diffs += Numeric.sum(drlv)
            contribs+= n_ind
            sumdrlv -= Numeric.sum(drlv)
            goodpks += n_ind
            

            if n_ind>300:
                print g.ubi
                print diff.shape
                print diff[:,:40]
                print drlv[:40]
                raise Exception("")
        #print "indexed %10d with mean error %f"%(goodpks,sumdrlv/goodpks)

        return 1e6*diffs/contribs

    stepsizes = {
        "wavelength" : 0.001,
        'y-center'   : 0.1,
        'z-center'   : 0.1,
        'distance'   : 0.1,
        'tilt-y'     : transform.radians(0.1),
        'tilt-z'     : transform.radians(0.1),
        'wedge'      : transform.radians(0.01),
        'chi'        : transform.radians(0.01),
        'translation' : 0.1
        }

    def applyargs(self,arg):
        for i in range(len(self.varylist)):
            item = self.varylist[i]
            value = arg[i]
            try:
                self.parameters[item]=value
            except:
                # Hopefully a crystal translation
                pass

    def printresult(self,arg):
        for i in range(len(self.varylist)):
            item = self.varylist[i]
            value = arg[i]
            try:
                self.parameters[item]=value
                print item,value
            except:
                # Hopefully a crystal translation
                pass
        

    def fit(self):
        import simplex
        guess = []
        inc = []
        for item in self.varylist:
            try:
                guess.append( float(self.parameters[item]) )
            except:
                # Hopefully a crystal translation
                print item
            try:
                inc.append(self.stepsizes[item])
            except:
                pass

        print "Start of simplex:"
        self.printresult(guess)
        print
            
        s=simplex.Simplex(self.gof,guess,inc)
      
        newguess,error,iter=s.minimize(maxiters=100)

        print
        print "End of simplex"
        self.printresult(newguess)
    
    def refineubis(self,quiet=True):
        
        for key in self.grains.keys():
            g = self.grains[key]
            grainname = key[0]
            scanname = key[1]
            if not quiet:
                print "%10s %10s"%(grainname,scanname),
            # Compute gv using current parameters
            self.compute_gv(grainname,scanname)
            #print self.gv.shape
            #print self.gv[0:10,:]
            g.ubi = self.refine(g.ubi,quiet=quiet)
        
            

        
        

if __name__ == "__main__":
    o=refinegrains()
    o.loadparameters("start.prm")
    scans  = ["0N.flt" ,"811N.flt" ,  "2211N.flt" , "3311N.flt" , "4811N.flt" ]
    for name in scans:
        o.loadfiltered(name)
    o.readubis("ubitest")
    o.generate_grains()
    o.tolerance = 0.1
    #o.varytranslations()
    for o.tolerance in [0.1,0.2,0.15,0.1,0.075,0.05]:
        o.refineubis(quiet=False)
    print "***"
    o.refineubis(quiet=False)
    o.varylist = ['y-center','z-center','distance','tilt-y','tilt-z','wedge','chi']
    o.fit()
    o.refineubis(quiet=False)
