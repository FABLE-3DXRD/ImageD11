

# Automatically adapted for numpy.oldnumeric Sep 06, 2007 by alter_code1.py


# ImageD11_v1.0.6 Software for beamline ID11
# Copyright (C) 2008  Jon Wright
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

from ImageD11 import transform, indexing, parameters
from ImageD11 import grain, columnfile, closest

import simplex

class refinegrains:
    
    """ 
    Class to refine the position and orientation of grains
    and also the detector parameters
    """

    # Default parameters
    pars = {"cell__a" : 4.1569162,
        "cell__b" :  4.1569162,
        "cell__c" : 4.1569162,
        "cell_alpha" : 90.0,
        "cell_beta" : 90.0,
        "cell_gamma" : 90.0,
        "cell_lattice_[P,A,B,C,I,F,R]" : 'P',
        "chi" : 0.0,
        "distance" : 7367.8452302,
        "fit_tolerance" : 0.5,
        "o11" : 1,
        "o12" : 0,
        "o21" : 0,
        "o22" : 1,
        "omegasign" : 1,
        "t_x" : 33.0198146824,
        "t_y" : 14.6384893741,
        "t_z" : 0.0 ,
        "tilt_x" : 0.0,
        "tilt_y" : 0.0623952920101,
        "tilt_z" : 0.00995011461696,
        "wavelength" : 0.2646,
        "wedge" : 0.0,
        "y_center" : 732.950204632,
        "y_size" : 4.6,
        "z_center" : 517.007049626,
        "z_size" : 4.6 }

    # Default stepsizes for the 
    stepsizes = {
        "wavelength" : 0.001,
        'y_center'   : 0.1,
        'z_center'   : 0.1,
        'distance'   : 0.1,
        'tilt_y'     : transform.radians(0.1),
        'tilt_z'     : transform.radians(0.1),
        'tilt_x'     : transform.radians(0.1),
        'wedge'      : transform.radians(0.1),
        'chi'        : transform.radians(0.1),
        't_x' : 0.1,
        't_y' : 0.1,
        't_z' : 0.1,
        }

    def __init__(self, tolerance = 0.01):
        """

        """
        self.tolerance=tolerance
        # list of ubi matrices (1 for each grain in each scan)
        self.grainnames = []
        self.ubisread = {}
        self.translationsread = {}
        # list of scans and corresponding data
        self.scannames=[]
        self.scantitles={}
        self.scandata={}
        # grains in each scan
        self.grains = {}
        self.grains_to_refine = []
        # ?
        self.drlv = None
        self.parameterobj = parameters.parameters(**self.pars)
        for k,s in self.stepsizes.items():
            self.parameterobj.stepsizes[k]=s

    def loadparameters(self, filename):
        self.parameterobj.loadparameters(filename)
        
    def saveparameters(self, filename):
        self.parameterobj.saveparameters(filename)

    def readubis(self, filename):
        """
        Read ubi matrices from a text file
        """
        try:
            ul = grain.read_grain_file(filename)
        except:
            print filename, type(filename)
            raise
        for i, g in enumerate(ul):
            name = filename + "_" + str(i)
            # Hmmm .... multiple grain files?
            name = i
            self.grainnames.append(i)
            self.ubisread[name] = g.ubi
            self.translationsread[name] = g.translation
        #print "Grain names",self.grainnames

    def savegrains(self, filename, sort_npks=True):
        """
        Save the refined grains
        
        """
        ks = self.grains.keys()
        # sort by number of peaks indexed to write out
        if sort_npks:
        #      npks
            gl = [ (len(self.grains[k].x), self.grains[k]) for k in ks ]
            gl.sort()
            gl = [ g[1] for g in gl[::-1] ]
        else:
            ks.sort()
            gl = [ self.grains[k] for k in ks ]
        grain.write_grain_file(filename, gl)


    def makeuniq(self, symmetry):
        """
        Flip orientations to a particular choice
        
        Be aware that if the unit cell is distorted and you do this,
        you might have problems...
        """
        from ImageD11.sym_u import find_uniq_u, getgroup
        g = getgroup( symmetry )()
        for k  in self.ubisread.keys():
            self.ubisread[k] = find_uniq_u(self.ubisread[k], g)
        for k in self.grains.keys():
            self.grains[k].ubi = find_uniq_u(self.grains[k].ubi, g)
            

    def loadfiltered(self, filename):
        """
        Read in file containing filtered and merged peak positions
        """
        col = columnfile.columnfile(filename)
        self.scannames.append(filename)
        self.scantitles[filename] = col.titles
        if not "drlv2" in col.titles:
            col.addcolumn( numpy.ones(col.nrows, numpy.float),
                           "drlv2" )
        if not "labels" in col.titles:
            col.addcolumn( numpy.ones(col.nrows, numpy.float)-2,
                           "labels" )
        self.scandata[filename] = col
        
        


    def generate_grains(self):
        t = numpy.array([ self.parameterobj.parameters[s]
                          for s in ['t_x', 't_y','t_z']] )
        for grainname in self.grainnames:
            for scanname in self.scannames:
                try:
                    gr = self.grains[(grainname,scanname)]
                except KeyError:
                    if (self.translationsread[grainname] == 0).all():
                        self.grains[(grainname,scanname)] = grain.grain(
                            self.ubisread[grainname],
                            translation=t)
                        self.grains[(grainname,scanname)].name = \
                            (str(grainname)+":"+scanname).replace(" ","_")
                    else:
                        self.grains[(grainname,scanname)] = grain.grain(
                            self.ubisread[grainname],
                            translation = self.translationsread[grainname] )
                        self.grains[(grainname,scanname)].name = \
                            (str(grainname)+":"+scanname).replace(" ","_")

                    

        for scanname in self.scannames:
            self.reset_labels(scanname)

    def reset_labels(self, scanname):
        
        try:
            x = self.scandata[scanname].xc
            y = self.scandata[scanname].yc
        except AttributeError:
            x = self.scandata[scanname].sc
            y = self.scandata[scanname].fc
        om = self.scandata[scanname].omega
        # only for this grain
        self.scandata[scanname].labels = self.scandata[scanname].labels*0 - 2 
        self.scandata[scanname].drlv2 = self.scandata[scanname].drlv2*0 + 1
        for g in self.grainnames:
            self.grains[(g,scanname)].x = x
            self.grains[(g,scanname)].y = y
            self.grains[(g,scanname)].om = om

            

    def compute_gv(self, grainname, scanname):
        """
        Makes self.gv refer be g-vectors computed for this grain in this scan
        """
        x  = self.grains[ ( grainname, scanname ) ].x
        y  = self.grains[ ( grainname, scanname ) ].y
        om = self.grains[ ( grainname, scanname ) ].om
        try:
            sign = self.parameterobj.parameters['omegasign']
        except:
            sign = 1.0
                
        self.tth,self.eta = transform.compute_tth_eta( numpy.array([x, y]),
                                      omega = om * sign,
                                      **self.parameterobj.parameters)
        self.gv = transform.compute_g_vectors(self.tth, self.eta, om*sign,
                                              float(self.parameterobj.parameters['wavelength']),
                                              self.parameterobj.parameters['wedge'],
                                              self.parameterobj.parameters['chi'])
        self.gv = self.gv.T
        # print "in compute_gv",self.gv[0]

    def refine(self, ubi, quiet=True):
        """
        Fit the matrix without changing the peak assignments
        
        """
        mat=ubi.copy()
        # print "In refine",self.tolerance, self.gv.shape
        self.npks, self.avg_drlv2 = closest.score_and_refine(mat, self.gv, 
                                                             self.tolerance)

        if not quiet:
            import math
            print "%-8d %.6f"%(self.npks,math.sqrt(self.avg_drlv2))

        #print self.tolerance
        # self.npks = closest.score_and_refine(mat, self.gv, self.tolerance)
        #tm = indexing.refine(ubi,self.gv,self.tolerance,quiet=quiet)
        #print ubi, tm,ubi-tm,mat-tm

        return mat

    def applyargs(self,args):
        self.parameterobj.set_variable_values(args)

    def printresult(self,arg):
        # print self.parameterobj.parameters
        # return
        for i in range(len(self.parameterobj.varylist)):
            item = self.parameterobj.varylist[i]
            value = arg[i]
            try:
                self.parameterobj.parameters[item]=value
                print item,value
            except:
                # Hopefully a crystal translation
                pass


    def gof(self,args):
        """
        <drlv> for all of the grains in all of the scans
        """
        self.applyargs(args)
        diffs = 0.
        contribs = 0.

        # defaulting to fitting all grains
        # print "refineing grains",self.grains_to_refine
        for key in self.grains_to_refine:

            g = self.grains[key]
            grainname = key[0]
            scanname = key[1]

            # Compute gv using current parameters
            # Keep labels fixed
            self.compute_gv(grainname,scanname)
            #print self.gv.shape
            #print self.gv[0:10,:]

            g.ubi = self.refine(g.ubi)
            self.npks # number of peaks it got

            diffs +=  self.npks*self.avg_drlv2
            contribs+= self.npks
        if contribs > 0:
            return 1e6*diffs/contribs
        else:
            print "No contribs???"
            return 1e6


    def fit(self, maxiters=100):
        """
        Fit the global parameters
        """
        self.assignlabels()
        guess = self.parameterobj.get_variable_values()
        inc = self.parameterobj.get_variable_stepsizes()
        self.printresult(guess)
        self.grains_to_refine = self.grains.keys()
        s=simplex.Simplex(self.gof, guess, inc)
        newguess,error,iter=s.minimize(maxiters=maxiters , monitor=1)
        print
        self.printresult(newguess)

    def getgrains(self):
        return self.grains.keys()


    def set_translation(self,gr,sc):
        self.parameterobj.parameters['t_x'] = self.grains[(gr,sc)].translation[0]
        self.parameterobj.parameters['t_y'] = self.grains[(gr,sc)].translation[1]
        self.parameterobj.parameters['t_z'] = self.grains[(gr,sc)].translation[2]

        
    def refinepositions(self, quiet=True,maxiters=100):
        self.assignlabels()
        ks  = self.grains.keys()
        ks.sort()
        # assignments are now fixed
        tolcache = self.tolerance
        self.tolerance = 1.0
        for key in ks:
            g = key[0]
            self.grains_to_refine = [key]
            self.parameterobj.varylist = [ 't_x', 't_y', 't_z' ]
            self.set_translation(key[0],key[1])
            guess = self.parameterobj.get_variable_values()
            inc =   self.parameterobj.get_variable_stepsizes()

            s =     simplex.Simplex(self.gof, guess, inc)

            newguess, error, iter = s.minimize(maxiters=maxiters,monitor=1)

            self.grains[key].translation[0] = self.parameterobj.parameters['t_x']
            self.grains[key].translation[1] = self.parameterobj.parameters['t_y'] 
            self.grains[key].translation[2] = self.parameterobj.parameters['t_z'] 
            print key,self.grains[key].translation,
            self.refine(self.grains[key].ubi,quiet=False)
        self.tolerance = tolcache

 
    def refineubis(self, quiet=True, scoreonly=False):
        #print quiet
        ks = self.grains.keys()
        ks.sort()
        if not quiet:
            print "%10s %10s"%("grainname","scanname"),
            print "npeak   <drlv> "
        for key in ks:
            g = self.grains[key]
            grainname = key[0]
            scanname = key[1]
            if not quiet:
                print "%10s %10s"%(grainname,scanname),
            # Compute gv using current parameters
            self.compute_gv(grainname, scanname)
            res = self.refine(g.ubi , quiet=quiet)
            if not scoreonly:
                g.ubi = res

    def assignlabels(self):
        """
        Fill out the appropriate labels for the spots
        """
        # print "Assigning labels"
        for s in self.scannames:
            self.scandata[s].labels = self.scandata[s].labels * 0 - 2 # == -1
            self.scandata[s].drlv2 = self.scandata[s].drlv2*0 + 1  # == 1
            nr = self.scandata[s].nrows
            int_tmp = numpy.zeros(nr , numpy.int32 )-1
            for g in self.grainnames:
                gr = self.grains[ ( g, s) ]
                self.set_translation( g, s)
                self.compute_gv( g, s )
                #print "about to assign"
                closest.score_and_assign( gr.ubi,
                                          self.gv,
                                          self.tolerance,
                                          self.scandata[s].drlv2,
                                          int_tmp,
                                          int(g) )
                #print "assigned"
            # Second loop after checking all grains
            self.scandata[s].labels = int_tmp * 1.
            for g in self.grainnames:
                gr = self.grains[ ( g, s) ]
                ind = numpy.compress(int_tmp == g,
                                     range(nr) )
                #print 'x',gr.x[:10]
                #print 'ind',ind[:10]
                try:
                    gr.x = numpy.take(self.scandata[s].xc , ind)
                    gr.y = numpy.take(self.scandata[s].yc , ind)
                except AttributeError:
                    gr.x = numpy.take(self.scandata[s].sc , ind)
                    gr.y = numpy.take(self.scandata[s].fc , ind)
                gr.om = numpy.take(self.scandata[s].omega , ind)
                self.grains[ ( g, s) ] = gr
                print "Grain",g,"Scan",s,"npks=",len(ind)
                #print 'x',gr.x[:10]



def test_benoit():
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


def test_nac():
    import sys
    o = refinegrains()
    o.loadparameters(sys.argv[1])
    print "got pars"
    o.loadfiltered(sys.argv[2])
    print "got filtered"
    o.readubis(sys.argv[3])
    print "got ubis"
    o.tolerance = float(sys.argv[4])
    print "generating"
    o.generate_grains()
    print "Refining posi too"
    o.refineubis(quiet = False , scoreonly = True)
    print "Refining positions too"
    o.refinepositions()
    print "Refining positions too"    
    o.refineubis(quiet = False , scoreonly = True)
    
    


if __name__ == "__main__":
    test_nac()
