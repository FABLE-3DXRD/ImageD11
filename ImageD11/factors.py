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


"""
Some routines for factor analysis, a bit slow, via the SVD
"""
import glob, numpy.oldnumeric as Numeric, LinearAlgebra, struct, logging

from ImageD11 import opendata

class factors:
    """
    Class for factor analysis
    """
    def __init__(self):
        self.obsdata = None
        self.svd = None
        self.gendata = None
        self.x = None 
        self.nfactors=0

    def generatedata(self):
        """
        compute the data from our svd
        """
        if self.svd is not None:
            l,s,r=self.svd # left , singularvals, right
            nf=self.nfactors
            I=Numeric.identity(nf,Numeric.Float)
            suse=I*s[:nf]
            ls = Numeric.matrixmultiply(l[:,:nf], suse)
            self.gendata=Numeric.matrixmultiply(ls,r[:nf,:])
            logging.debug("generated data.shape"+str(self.gendata.shape))


    def factorsvd(self):
        """
        compute the svd
        """
        logging.debug("In svd")
        if self.obsdata is not None:
            logging.debug("Calling svd")
            self.svd = LinearAlgebra.singular_value_decomposition(self.obsdata)

    def savesvd(self,filename):
        """
        Save the time consuming svd step
        """
        if self.svd is not None:
            l,s,r=self.svd
            out=open(filename,"wb")
            out.write(struct.pack("lllll",
                                  l.shape[0],
                                  l.shape[1],
                                  s.shape[0],
                                  r.shape[0],
                                  r.shape[1]))
            print len(l.astype(Numeric.Float).tostring())
            out.write(l.astype(Numeric.Float).tostring())
            out.write(s.astype(Numeric.Float).tostring())
            out.write(r.astype(Numeric.Float).tostring())
            out.close()

    def loadsvd(self,filename):
        """
        Save the time consuming svd step
        """
        out=open(filename,"rb")
        dims=struct.unpack("lllll",out.read(struct.calcsize("lllll")))
        print dims,8*dims[0]*8*dims[1]
        l=Numeric.fromstring(out.read(8*dims[0]*dims[1]),Numeric.Float)
        s=Numeric.fromstring(out.read(8*dims[2]          ),Numeric.Float)
        r=Numeric.fromstring(out.read(8*dims[3]*dims[4]),Numeric.Float)
        print l.shape,s.shape,r.shape,dims
        l=Numeric.reshape(l,(dims[0],dims[1]))
        r=Numeric.reshape(r,(dims[3],dims[4]))
        self.svd=(l,s,r)



    def setnfactors(self,n):
        """Decide on the number of factors in the data"""
        self.nfactors=n
        print "Number of factors set to",self.nfactors

    def loadchis(self,filename):
        """
        Glob for filenames
        """
        fl=glob.glob(filename[:-8]+"????"+".chi")
        fl.sort()
        print "Number of chi files is:",len(fl)
        dl=[opendata.openchi(f).data[:,1] for f in fl]
        self.obsdata=Numeric.array(dl)
        self.x=opendata.openchi(fl[0]).data[:,0]
        print self.x.shape,self.obsdata.shape

    def saveobsdata(self,filename):
        """
        Save in binary format?
        """
        out=open(filename,"wb")
        out.write(struct.pack("lll",self.x.shape[0],
                                    self.obsdata.shape[0],
                                    self.obsdata.shape[1]))
        out.write(self.x.astype(Numeric.Float).tostring())
        out.write(self.obsdata.astype(Numeric.Float).tostring())
        out.close()

    def readobsdata(self,filename):
        """ Reads observed data from a binary format """
        infile=open(filename,"rb")
        sizes=struct.unpack("lll",infile.read(struct.calcsize("lll")))
        # Type is Float therefore 8 bytes per item
        print sizes
        self.x=Numeric.fromstring(infile.read(sizes[0]*8),Numeric.Float)
        self.obsdata=Numeric.fromstring(infile.read(8*sizes[1]*sizes[2]),
                                        Numeric.Float)
        print self.obsdata.shape
        self.obsdata=Numeric.reshape(self.obsdata,(sizes[1],sizes[2]))
        print self.x.shape,self.obsdata.shape


if __name__=="__main__":
    o=factors()
    import sys
    if sys.argv[1].find("chi")>-1:
        o.loadchis(sys.argv[1])
    else:
        o.readobsdata(sys.argv[1])
    print dir(o)
    o.setnfactors(int(sys.argv[2]))
    o.factorsvd()
    o.generatedata()
