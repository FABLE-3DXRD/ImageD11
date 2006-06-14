
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

from Numeric import transpose, matrixmultiply, Float, Int, floor, sum, zeros

class grain:
    def __init__(self,ubi,translation=None):
        self.ubi=ubi
        if translation==None:
            self.translation = zeros(3,Float)
        else:
            self.translation = array(translation)

    def ubitocellpars(self):
        """
        convert ubi matrix to unit cell
        """
        g=matrixmultiply(self.ubi,transpose(self.ubi))
        a=sqrt(g[0,0])
        b=sqrt(g[1,1])
        c=sqrt(g[2,2])
        from math import acos, degrees
        alpha=degrees(acos(g[1,2]/b/c))
        beta =degrees(acos(g[0,2]/a/c))
        gamma=degrees(acos(g[0,1]/a/b))
        return a,b,c,alpha,beta,gamma
 
    def refine(self,gv,tol,quiet=True):
        """
        Refine an orientation matrix and rescore it.
        
        From Paciorek et al Acta A55 543 (1999)
        UB = R H-1
        where:
        R = sum_n r_n h_n^t
        H = sum_n h_n h_n^t
        r = g-vectors
        h = hkl indices
        """
        #      print "Orientation and unit cell refinement of"
        #      print "UBI\n",UBI
        #      print "Scores before",self.score(UBI)
        # Need to find hkl indices for all of the peaks which are indexed
        h=matrixmultiply(self.ubi,transpose(gv))
        hint=floor(h+0.5).astype(Int) # rounds down
        diff=h-hint
        drlv2=sum(diff*diff,0)
        tol = float(tol)
        tol = tol*tol
        # Only use peaks which are assigned to rings for refinement
        ind = compress( less(drlv2,tol) , arange(gv.shape[0]) )
        scoreb4=ind.shape[0]
        contribs = take(drlv2,ind)
        try:
            fitb4=math.sqrt(sum(contribs)/contribs.shape[0])
            if not quiet:
                print "Fit before refinement %.8f %5d"%(fitb4,contribs.shape[0]),
        except:
            print "No contributing reflections for\n",UBI
            raise
        drlv2_old=drlv2
        R=zeros((3,3),Float)
        H=zeros((3,3),Float)
        for i in ind:
            r = gv[i,:]
            k = hint[:,i].astype(Float)
            #           print r,k
            R = R + outerproduct(r,k)
            H = H + outerproduct(k,k)
        from LinearAlgebra import inverse
        try:
            HI=inverse(H)
            UBoptimal=matrixmultiply(R,HI)
            UBIo=inverse(UBoptimal)
        except:
            # A singular matrix - this sucks.
            UBIo=UBI
        h=matrixmultiply(UBIo,transpose(gv))
        hint=floor(h+0.5).astype(Int) # rounds down
        diff=h-hint
        drlv2=sum(diff*diff,0)
        ind = compress( less(drlv2,tol), arange(gv.shape[0]) )
        scorelastrefined=ind.shape[0]
        contribs = take(drlv2,ind)
        try:
            fitlastrefined=math.sqrt(sum(contribs)/contribs.shape[0])
            if not quiet:
                print "after %.8f %5d"%(fitlastrefined,contribs.shape[0])
        except:
            print "\n\n\n"
            print "No contributing reflections for\n",UBI
            print "After refinement, it was OK before ???"
            print "\n\n\n"
            return self.ubi
        #raise
        #      for i in ind:
        #         print "( %-6.4f %-6.4f %-6.4f ) %12.8f %12.8f"%(h[0,i],h[1,i],h[2,i],sqrt(drlv2[i]),sqrt(drlv2_old[i]))
        #      print UBIo
        #      print "Scores after", self.score(UBIo,self.hkl_tol)
        #      print "diff\n",UBI-UBIo
        #      print "Mean drlv now",sum(sqrt(drlv2))/drlv2.shape[0],
        #      print "Mean drlv old",sum(sqrt(drlv2_old))/drlv2_old.shape[0]
        return UBIo

