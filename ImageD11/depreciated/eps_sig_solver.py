
from __future__ import print_function

# Get Strain/Stress from ImageD11 UBI/map files
# Copyright (C) 2015 Younes ELHACHI
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

import numpy as np
import math
from ImageD11 import transform, unitcell, columnfile
from ImageD11.parameters import par, parameters
from ImageD11.grain import read_grain_file
try:
    from FitAllB.conversion import grain2sample, strain2stress,formStiffnessMV
except:
    print( "You need to install FitAllB")
    print( "You will get error messages if you try to compute strain!")
from xfab.tools import ubi_to_u_and_eps,ubi_to_cell

def readubis(ubifile):
    """read ubifile and return a list of ubi arrays """
    f = open(ubifile, "r")
    ubisread = []
    u = []
    for line in f:
        if line[0]=="#":
            continue
        vals = [ float(x) for x in line.split() ]
        if len(vals) == 3:
            u = u + [vals]
        if len(u)==3:
            ubisread.append(np.array(u))
            u = []
    f.close()
    return ubisread

def write_ubi_file(filename, ubilist):
    """ save 3x3 matrices into file """
    f=open(filename,"w")
    for u in ubilist:
        f.write("%f %f %f\n"  %(u[0][0],u[0][1],u[0][2]))
        f.write("%f %f %f\n"  %(u[1][0],u[1][1],u[1][2]))
        f.write("%f %f %f\n\n"%(u[2][0],u[2][1],u[2][2]))
    f.close() 



class solver:
    """
    A class for getting strain and stress tensors
    """
    def __init__(self,
            unitcell=None,
            ubis=None,
            crystal_symmetry=None,
            c11=None,c12=None,c13=None,c14=None,c15=None,c16=None,
            c22=None,c23=None,c24=None,c25=None,c26=None,
            c33=None,c34=None,c35=None,c36=None,
            c44=None,c45=None,c46=None,
            c55=None,c56=None,
            c66=None):
        """
        unitcell would be a list of six elements [a, b, c, alpha, beta, gamma]
        ubis would be a list of orientation matrices as by ImageD11 convention
        symmetry would be isotropic, cubic, tetragonal_high... (see FitAllB.conversion) to form the stiffness tensor C
        The rest of the arguments are parameters.
        """
        
        self.cell__a = None
        self.cell__b = None
        self.cell__c = None
        self.cell_alpha = None
        self.cell_beta = None
        self.cell_gamma = None
        if unitcell is not None:
            if len(unitcell)==6:
                self.cell__a = unitcell[0]
                self.cell__b = unitcell[1]
                self.cell__c = unitcell[2]
                self.cell_alpha = unitcell[3]
                self.cell_beta = unitcell[4]
                self.cell_gamma = unitcell[5]
            else:
                raise Exception("The unit cell must be defined by six parameters!")
        self.ubis=ubis
        self.crystal_symmetry=crystal_symmetry
        self.c11=c11
        self.c12=c12
        self.c13=c13
        self.c14=c14
        self.c15=c15
        self.c16=c16
        self.c22=c22
        self.c23=c23
        self.c24=c24
        self.c25=c25
        self.c26=c26
        self.c33=c33
        self.c34=c34
        self.c35=c35
        self.c36=c36
        self.c44=c44
        self.c45=c45
        self.c46=c46
        self.c55=c55
        self.c56=c56
        self.c66=c66
        self.parameterobj = parameters(cell__a=self.cell__a,
                                       cell__b=self.cell__b,
                                       cell__c=self.cell__c,
                                       cell_alpha=self.cell_alpha,
                                       cell_beta=self.cell_beta,
                                       cell_gamma=self.cell_gamma,                                                 
                                       crystal_symmetry=self.crystal_symmetry,
                                       c11=self.c11,c12=self.c12,c13=self.c13,c14=self.c14,c15=self.c15,c16=self.c16,
                                       c22=self.c22,c23=self.c23,c24=self.c24,c25=self.c25,c26=self.c26,
                                       c33=self.c33,c34=self.c34,c35=self.c35,c36=self.c36,
                                       c44=self.c44,c45=self.c45,c46=self.c46,
                                       c55=self.c55,c56=self.c56,
                                       c66=self.c66)
        self.epsilon=[]
        self.sigma=[]

        
    def loadmap(self,filename):
        try:
            self.map=read_grain_file(filename)
            self.ubis=[x.ubi for x in self.map]
        except:
            print("error when reading %s\n",filename)
            raise
    
    
    def loadpars(self,filename=None):
        if filename is not None:
            self.parameterobj.loadparameters(filename)
        self.parameterobj.update_other(self)
        #update also the unitcell list (because the element are included in parameterobj but not the list):
        #self.unitcell=[self.cell__a, self.cell__b, self.cell__c, self.cell_alpha, self.cell_beta, self.cell_gamma]

    def updateparameters(self):
        self.savepars()
        self.pars=self.parameterobj.parameters
        #update also the unitcell list (because the element are included in parameterobj but not the list):
        #self.unitcell=[self.cell__a, self.cell__b, self.cell__c, self.cell_alpha, self.cell_beta, self.cell_gamma]

    def savepars(self,filename=None):
        self.parameterobj.update_yourself(self)
        if filename is not None:
            self.parameterobj.saveparameters(filename)
    
        
    def unitcell(self):
        return [self.cell__a, self.cell__b, self.cell__c, self.cell_alpha, self.cell_beta, self.cell_gamma]
    
    def setunitcell(self,uc):
        """    this is used to set all the unit cell elements in one shot    """
        self.cell__a = uc[0]
        self.cell__b = uc[1]
        self.cell__c = uc[2]
        self.cell_alpha = uc[3]
        self.cell_beta = uc[4]
        self.cell_gamma = uc[5]
        
    
    def MVStiffness(self):
        return formStiffnessMV(crystal_system=self.crystal_symmetry,
                               c11=self.c11,c12=self.c12,c13=self.c13,c14=self.c14,c15=self.c15,c16=self.c16,
                               c22=self.c22,c23=self.c23,c24=self.c24,c25=self.c25,c26=self.c26,
                               c33=self.c33,c34=self.c34,c35=self.c35,c36=self.c36,
                               c44=self.c44,c45=self.c45,c46=self.c46,
                               c55=self.c55,c56=self.c56,
                               c66=self.c66)
          

    def compute_write_eps_sig(self,outputfile):
        """    Compute strain and stress in crystal and sample co-ordinates system    """
        
        if self.ubis is not None:
            
            writestress = True
            
            f = open(outputfile,'w')
            ''' the used parameters will be the header of the output file'''
            for k,v in sorted(self.parameterobj.parameters.items()):
                f.write(("%s %s\n")%(k,v))
            ''' write titles'''
            f.write("##############################################\n")
            f.write("cell__a cell__b cell__c cell_alpha cell_beta cell_gamma u11 u12 u13 u21 u22 u23 u31 u32 u33 ")
            f.write("eps11_c eps22_c eps33_c eps12_c eps13_c eps23_c eps11_s eps22_s eps33_s eps12_s eps13_s eps23_s ")
            f.write("sig11_c sig22_c sig33_c sig12_c sig13_c sig23_c sig11_s sig22_s sig33_s sig12_s sig13_s sig23_s\n")
            
            '''this is the part where we compute strain and stress'''
            for ubi in self.ubis:
                U, eps = ubi_to_u_and_eps(ubi, self.unitcell())
                #U, eps = ubi_to_u_and_eps(ubi, ubi_to_cell(self.ubis[0]))
                epsM = [ [ eps[0], eps[1], eps[2] ],        #write the strain tensor list as a matrix
                         [ eps[1], eps[3], eps[4] ],
                         [ eps[2], eps[4], eps[5] ] ]
                epsS = np.dot( U, np.dot( epsM, U.T ) )     #epsilon in sample co-ordinates
                sigM= np.empty((3,3))
                sigS= np.empty((3,3))
                try:
                    sigM = strain2stress( np.array(epsM), self.MVStiffness() )      #write the stress tensor as a symmetric matrix in crystal co-ordinates
                    sigS = np.dot( U, np.dot( sigM, U.T ) )     #sigma in sample co-ordinates
                except:
                    print("couldn't compute stress! please check the crystal_symmetry parameters and elastic constants")
                    writestress = False
                
                
                ''' writing down the results'''
                
                f.write(("%f "*6)%tuple(ubi_to_cell(ubi)))
                f.write(("%f "*9)%tuple(U.ravel()))
                ligne = ""
                for i,j in [(0,0),(1,1),(2,2),(0,1),(0,2),(1,2)]:
                    ligne = ligne + " " + str(100.*np.array(epsM)[i,j])
                for i,j in [(0,0),(1,1),(2,2),(0,1),(0,2),(1,2)]:
                    ligne = ligne + " " + str(100.*epsS[i,j])
                if writestress==True:
                    for i,j in [(0,0),(1,1),(2,2),(0,1),(0,2),(1,2)]:
                        ligne = ligne + " " + str(np.array(sigM)[i,j])
                    for i,j in [(0,0),(1,1),(2,2),(0,1),(0,2),(1,2)]:
                        ligne = ligne + " " + str(sigS[i,j])
                
                f.write(ligne[1:])
                f.write("\n")
            f.close()
