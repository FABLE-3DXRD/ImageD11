
# coding: utf-8
from __future__ import print_function, division

import numpy as np
import ImageD11.finite_strain, ImageD11.unitcell, ImageD11.parameters
from ImageD11.grain import read_grain_file


# Compute strain and stress from a list of ImageD11 UBI matrices 
# Copyright (C) 2023 Jean-Baptiste Jacob
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



class EpsSigSolver:
    """ A class to handle strain and stress computations on a list of ubis. """
    def __init__(self,
                 name = 'default',
                 symmetry = None,
                 stress_unit = 'GPa',
                 UBI_list = [],
                 unitcell = None,
                 c11=None,c12=None,c13=None,c14=None,c15=None,c16=None,
                 c22=None,c23=None,c24=None,c25=None,c26=None,
                 c33=None,c34=None,c35=None,c36=None,
                 c44=None,c45=None,c46=None,
                 c55=None,c56=None,
                 c66=None):

        """
        Initialize the EpsSigSolver class.

        Parameters:
        -----------
        cij (float)        : elastic constants of the crystal
        Cij_symmetry (str) : symmetry considered for the Stiffness and Compliance matrices. Should be one of the following:
                               'cubic', 'trigonal_high', 'trigonal_low', 'tetragonal_high', 'tetragonal_low', 'hexagonal', 'orthorhombic', 'monoclinic', 'triclinic'
                              
        unitcell (array_like)   : Unstrained unit cell parameters [a, b, c, alpha,beta, gamma]
        UBI_list (list of 3x3 arrays) : List of real-space unit cell vectors (ubi in ImageD11).
        """
        
        assert symmetry in Cij_symmetry.keys(), 'symmetry not recognized!'
        
        self.phase_name = name
        self.unitcell = unitcell

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
                
        self.parameterobj = ImageD11.parameters.parameters(cell__a=self.cell__a,
                                               cell__b=self.cell__b,
                                               cell__c=self.cell__c,
                                               cell_alpha=self.cell_alpha,
                                               cell_beta=self.cell_beta,
                                               cell_gamma=self.cell_gamma,                                                 
                                               c11=c11,c12=c12,c13=c13,c14=c14,c15=c15,c16=c16,
                                               c22=c22,c23=c23,c24=c24,c25=c25,c26=c26,
                                               c33=c33,c34=c34,c35=c35,c36=c36,
                                               c44=c44,c45=c45,c46=c46,
                                               c55=c55,c56=c56,
                                               c66=c66)
        
        self.symmetry = symmetry
        self.stress_unit = 'GPa'
        self.Cij_symmetry = Cij_symmetry[self.symmetry]
        self.Cij = self.form_stiffness_tensor()
        if not np.alltrue(self.Cij == 0):
            self.Sij = np.linalg.inv(self.Cij)
        
        self.UBIs = UBI_list
        self.F_list = None

    @property
    def nUBIs(self):
        return len(self.UBIs)
        
    def __str__(self):
        return ("EpsSigSolver:\n phase name: {phase_name}\n reference unitcell: {unitcell}\n symmetry:" +\
               "{symmetry}\n unit:{stress_unit}\n Stiffness:\n {Cij}\n  n ubis: {nUBIs)}").format(self)

    # Load / save / update functions for parameters (from former eps_sig_solver.py)
    ########################################
    def loadmap(self,filename):
        try:
            self.glist = ImageD11.grain.read_grain_file(filename)
            self.UBIs=[g.ubi for g in self.glist]
        except:
            print("error when reading %s\n",filename)
            raise
            
            
    def loadpars(self,filename=None):
        if filename is not None:
            self.parameterobj.loadparameters(filename)
        self.parameterobj.update_other(self)
        #update also the unitcell list (because the element are included in parameterobj but not the list):
        self.unitcell=[self.cell__a, self.cell__b, self.cell__c, self.cell_alpha, self.cell_beta, self.cell_gamma]

    def updateparameters(self):
        self.savepars()
        self.pars=self.parameterobj.parameters
        #update also the unitcell list (because the element are included in parameterobj but not the list):
        self.unitcell=[self.cell__a, self.cell__b, self.cell__c, self.cell_alpha, self.cell_beta, self.cell_gamma]

    def savepars(self,filename=None):
        self.parameterobj.update_yourself(self)
        if filename is not None:
            self.parameterobj.saveparameters(filename)
            
    
    
    
    # Methods for Strain & stress computation
    #########################################
    
    def form_stiffness_tensor(self):
        """
        Form the Voigt 6x6 stiffness tensor (Cij) from elastic constants (cij), using the symmetry specified in self.symmetry.
        By definition, this is only valid for strain in grain coordinates, which must be provided as a 6-component vector 
        in voigt notation (e11, e22, e33, 2.e23, 2.e13, 2.e12). Stress in lab coordinates can be obtained by rotating the 3x3 stress
        tensor from the grain to the lab coordinate system (see strain2stress_lab for details).
        
        Return
        ----------
        Cij : 6x6 matrix containing the elastic components
        """
        Cij = np.zeros((6,6))
        pattern = self.Cij_symmetry  # pattern for the stiffness matrix. From Mouhat & Coudert (2014). Necessary and sufficient elastic stability conditions in various crystal systems. Physical Review

        
        parlist = 'c11,c22,c33,c44,c55,c66,c12,c13,c14,c15,c16,c23,c24,c25,c26,c34,c35,c36,c45,c46,c56'.split(',')
        
        for i,parname in enumerate(parlist):
            v = self.parameterobj.parameters[parname]
            if v is None:
                continue
            if self.symmetry in ['hexagonal','trigonal_high','trigonal_low'] and parname == 'c66':
                v = 0.5 * (self.parameterobj.parameters['c11'] - self.parameterobj.parameters['c12'])
                
            c_ij = np.where( np.abs(pattern) == i+1, np.sign(pattern) * v, 0 )
            Cij += c_ij
        return Cij 
    
    
    
    def compute_Deformation_Gradient_Tensors(self):
        """ 
        compute the deformation gradient tensor F for all ubis
        
        Return
        --------
        New instance F_list added to EpsSiSolver, containing a list of F tensors
        """
        F_list = []
        B0 = ImageD11.unitcell.unitcell(self.unitcell).B
        
        for ubi in self.UBIs:
            F = ImageD11.finite_strain.DeformationGradientTensor(ubi, B0)
            F_list.append(F)
        self.F_list = F_list

    

    def strain2stress_Ref(self, m=1):
        """
        Compute elastic strain and stress in grain reference coordinates for all ubis, using the stiffness matrix in self.Cij.
        Returns strain and stress as two lists of 3x3 symmetric tensors 'eps_Ref' and 'sigma_Ref'
        
        Parameters
        ----------
        m (float)    : exponent for the Seth-Hill finite strain tensors E = 1/2m (U^2m - I). see documentation in ImageD11.finite_strain.py
        for more detail
        
        Return
        ----------
        New instances 'eps_Ref' and 'sigma_Ref' added to  EpsSigSolver, containing respectively strain and stress as lists of 3x3 tensors
        """
        
        assert self.F_list is not None, 'No deformation gradient tensors to process. Run "self.compute_Deformation_Gradient_Tensors" first' 
        
        self.eps_Ref = []
        self.sigma_Ref = []
        
        for F in self.F_list:
            
            # compute strain and stress
            eRef = F.finite_strain_ref(m)
            eRef_v = full_3x3_to_vector(eRef, 'voigt', is_strain=True)
            sRef_v = np.dot(self.Cij, eRef_v)
            sRef = vector_to_full_3x3(sRef_v, 'voigt', is_strain=False)
            
            # append strain and stress to lists
            self.eps_Ref.append(eRef)
            self.sigma_Ref.append(sRef)
                
            
            
    def strain2stress_Lab(self, m=1, debug=0):
        """
        Compute elastic strain and stress in Lab coordinates for all ubis, using the stiffness matrix in self.Cij. 
        Computation is done first in the grain coordinate system, and then stress in lab coordinates is obtained by
        rotating the 3x3 stress tensor from the grain to the lab coordinate system using the following transormation
        σ' = RT.σ.R where R is the rotation matrix yielded by the polar decomposition of the finite deformation
        gradient tensor F.
        
        Returns strain and stress as two lists of 3x3 symmetric tensors 'eps_Lab' and 'sigma_Lab'. 
        
        Parameters
        ----------
        m (float)    : exponent for the Seth-Hill finite strain tensors E = 1/2m (U^2m - I). see documentation in ImageD11.finite_strain.py
        for more detail
         
        Return
        ----------
        New instances 'eps_Lab' and 'sigma_Lab' added to  EpsSigSolver, containing respectively strain and stress as lists of 3x3 tensors
        """
        
        self.eps_Lab = []
        self.sigma_Lab = []
        
        for F in self.F_list:
            
            # compute strain in lab coordinates
            eLab = F.finite_strain_lab(m)
            
            # compute strain and stress in grain coordinates
            eRef = F.finite_strain_ref(m)
            eRef_v = full_3x3_to_vector(eRef, 'voigt', is_strain=True)
            sRef =  vector_to_full_3x3( np.dot(self.Cij, eRef_v), 'voigt', is_strain=False )
            
            # Rotate stress to lab coordinates           
            sLab = F.U.dot(sRef).dot(F.U.T)
            
            # check rotation is working ok for strain
            assert np.allclose(eLab, F.U.dot(eRef).dot(F.U.T)), 'eRef and eLab are inconsistent. Rotation matrix might be dodgy'
            
            # Check stress tensor has the same eigenvalues in lab and grain coordinates
            eigv_Ref,_ = eigen_decomposition(sRef)
            eigv_Lab,_ = eigen_decomposition(sLab)
            assert np.allclose(eigv_Ref, eigv_Lab), 'sRef and sLab are inconsistent'
            
            # append strain and stress to lists
            self.eps_Lab.append(eLab)
            self.sigma_Lab.append(sLab)     
            
            
    def convert_tensors_to_vecs(self, output_format = 'voigt', debug=0):
        """
        convert strain and stress tensors to 6-component vectors in different possible formats
        
        output_format (str) : select one of the output formats below:
        
        default : e11, e22, e33, e23, e13, e12          | s11, s22, s33, s23, s13, s12
        xfab    : e11, e12, e13, e22, e23, e33          | s11, s12, s13, s22, s23, s33
        mandel  : e11, e22, e33, √2.e22, √2.e23, √2.e33 | s11, s22, s33, √2.s22, √2.s23, √2.s33
        voigt   : e11, e22, e33, 2.e23, 2.e13, 2.e12    | s11, s22, s33, s23, s13, s12
        
        Return
        ----------
        New instances 'eps_xx_fmt' and 'sigma_xx_fmt' added to  EpsSigSolver, containing respectively strain and stress as 6-component vectors
        instances names refer explicitely to the format specified in output_format, e.g. esp_Ref_voigt, etc.
        """
        formats = 'default,voigt,mandel,xfab'.split(',')
        assert output_format in formats, 'output format not recognized!'
        
        # select all strain and stress tensors list
        dnames = [attr for attr in dir(self) if any([attr.startswith(s) for s in ['eps','sigma']]) ]
         # filter out all data that begins with 'eps' or 'sigma' but are not strain or stress tensors  
        dnames = [d for d in dnames if any([d.endswith(s) for s in ['Lab','Ref','_d']]) ] 
        
        # stop if no strain / stress tensor list fond
        if len(dnames) == 0:
            print('No strain or stress tensors list found')
            return
        
        # convert all tensors lists to selected vector format
        for d in dnames:
            if debug:
                print(d)
            tensor_list = self.__getattribute__(d)
            assert np.all([T.shape == (3,3) for T in tensor_list]), 'data in input list are not 3x3 tensors' 
            
            if 'eps' in d:
                strain = True
            else:
                strain = False
                
            vector_list = [full_3x3_to_vector(T, output_format, is_strain=strain) for T in tensor_list]
            outname = '_'.join([d,output_format])
            
            setattr(self, outname, vector_list)
            
            
    def deviatoric_tensors(self, debug=0):
        """
        compute deviatoric component for all strain and stress tensors lists
        
        Return
        ----------
        New instances 'eps_xx_d' and 'sigma_xx_d' added to  EpsSigSolver, containing lists of deviatoric 3x3 tensors 
        """  
        
        # select all strain and stress tensors list that are not deviatoric
        dnames = [attr for attr in dir(self) if any([attr.startswith(s) for s in ['eps','sigma']]) ]
        dnames = [d for d in dnames if any([d.endswith(s) for s in ['Lab','Ref']]) ] 
        
        # stop if no strain / stress tensor list fond
        assert len(dnames) > 0, 'No strain or stress tensors list found'
        
        # compute deviatoric tensors
        for d in dnames:
            if debug:
                print(d)
            
            tensor_list = self.__getattribute__(d)
            assert np.all([T.shape == (3,3) for T in tensor_list]), 'data in input list are not 3x3 tensors' 
            tensors_dev = [deviatoric(T) for T in tensor_list]
            
            setattr(self, d+'_d', tensors_dev)
            
            
    
    def invariant_props(self, dname):
        # NOTE : not sure about the expression of von Mises strain. In any case it is related to √J2 by a multiplication factor k,
        # but it seems to be different from the definition of von Mises stress √(3.J2).
        # see https://www.continuummechanics.org/vonmisesstress.html
        """
        compute invariant properties for selected data column
        compute invariant properties for selected data column: volumetric strain / pressure (-I1/3) and von mises strain /stress (√3.J2)
        
        dname (str) : name of the input data column. Must be a non-deviatoric 3x3 tensors 
        
        Returns
        ----------
        New instances added to EpsSigSolver, containing list of floats
        if strain tensor in input
        dname+'_vol' : volumetric strain 
        dname+'_vM'  : von Mises strain (√2.J2)
        if stress tensor in input:
        dname+'_P_hyd'   : hydrostatic Pressure (if stress tensor in input)
        dname+'_vM'  : von Mises stress (√3.J2)
        """
        assert dname in dir(self), 'dname not recognized'
        assert '_d_' not in dname, 'tensor is deviatoric. Please use the non-deviatoric tensor'

        tensor_list = self.__getattribute__(dname)
        assert np.all([T.shape == (3,3) for T in tensor_list])

        tensor_list_dev = [deviatoric(T) for T in tensor_list]
        Invts = [invariants(T) for T in tensor_list]
        Invts_dev = [invariants(T) for T in tensor_list_dev]

        Inv1 = [-i[0]/3 for i in Invts]
        Inv2 = [np.sqrt(3*i[1]) for i in Invts_dev]

        if 'eps' in dname:
            setattr(self, dname+'_vol', Inv1)
        else:
            setattr(self, dname+'_P_hyd', Inv1)
        setattr(self, dname+'_vM', Inv2)
        
        
        
    def compute_principal_components(self, dname):
        """
        compute principal components resulting from eigen decomposition of the strain / stress tensor for selected data column
        
        dname (str) : name of the input data column 
        
        Returns
        ----------
        New instances dname+'_eigvals' , dname+'_eigvecs' added to EpsSigSolver, containing respectively a list of eigenvalues and eigenvectors
        
        eigvals (1x3 array) : Principal components in decreasing order from the largest positive to largest negative
                              for strain : ε1 > ε2 > ε3 (positive strain = elongation)
                              for stress σ3 > σ2 > σ1 (positive stress  = tension)
        eigvecs (3x3 array) : Normalized principal component vectors by columns sorted accordingly to eigvals
        """
        assert dname in dir(self), 'dname not recognized'
        
        tensor_list = self.__getattribute__(dname)
        assert np.all([T.shape == (3,3) for T in tensor_list]), 'data in input list are not 3x3 tensors' 
        
        PC = [eigen_decomposition(T) for T in tensor_list]
        eigvals = [p[0] for p in PC]
        eigvecs = [p[1] for p in PC]
        
        setattr(self, dname+'_eigvals', eigvals)
        setattr(self, dname+'_eigvecs', eigvecs)
            
            
            
            
# GENERAL FUNCTIONS
################################################################################

def full_3x3_to_vector(T, output_format = 'default', is_strain=True):
    """
    Form a 6 component strain or stress vector from a 3x3 tensor, using specified convention
    
    Parameters
    -----------
    T (array_like) : strain or stress tensor of the form [e11 e12 e13]
                                                         [e12 e22 e23]
                                                         [e13 e23 e33]
                                                         
    is_strain (bool) : specify whether T is a strain or stress tensor
    
    output_format (str) : specify convention for the 6-component vector: must be one of the following: 'default', 'xfab', 'mandel', 'voigt'

    default : e11, e22, e33, e23, e13, e12          | s11, s22, s33, s23, s13, s12
    xfab    : e11, e12, e13, e22, e23, e33          | s11, s12, s13, s22, s23, s33
    mandel  : e11, e22, e33, √2.e22, √2.e23, √2.e33 | s11, s22, s33, √2.s22, √2.s23, √2.s33
    voigt   : e11, e22, e33, 2.e23, 2.e13, 2.e12    | s11, s22, s33, s23, s13, s12
    """
 
    assert output_format in 'default,voigt,mandel,xfab'.split(','), 'Format not recognized'
    if output_format == 'default' or not is_strain:
        return np.array([T[0, 0], T[1, 1], T[2, 2], T[1, 2], T[0, 2], T[0, 1]])
    
    elif output_format == 'xfab':
        return np.array([T[0, 0], T[0, 1], T[0, 2], T[1, 1], T[1, 2], T[2, 2]])
    
    elif output_format == 'mandel':
        sqr2 = np.sqrt(2)
        return np.array([T[0, 0], T[1, 1], T[2, 2], sqr2*T[1, 2], sqr2*T[0, 2], sqr2*T[0, 1]])
    
    else:
        return np.array([T[0, 0], T[1, 1], T[2, 2], 2*T[1, 2], 2*T[0, 2], 2*T[0, 1]])
        


def vector_to_full_3x3(vec, input_format='default', is_strain=True):
    """
    Reconstruct a 3x3 strain or stress tensor from a 6-component vector, using specified convention
    
    Parameters
    -----------
    vec (array_like) : 6-component strain or stress vector
    
    is_strain (bool) : specify whether the intput vector is a strain or stress vector
    
    input_format (str) : specify convention for the 6-component vector: must be one of the following: 'default', 'xfab', 'mandel', 'voigt'
    
    default : e11, e22, e33, e23, e13, e12          | s11, s22, s33, s23, s13, s12
    xfab    : e11, e12, e13, e22, e23, e33          | s11, s12, s13, s22, s23, s33
    mandel  : e11, e22, e33, √2.e22, √2.e23, √2.e33 | s11, s22, s33, √2.s22, √2.s23, √2.s33
    voigt   : e11, e22, e33, 2.e23, 2.e13, 2.e12    | s11, s22, s33, s23, s13, s12
    """
 
    assert input_format in 'default,voigt,mandel,xfab'.split(','), 'Format not recognized'
    
    if input_format == 'default' or not is_strain:
        return np.array([[vec[0], vec[5], vec[4]],
                         [vec[5], vec[1], vec[3]],
                         [vec[4], vec[3], vec[2]]])
    
    elif input_format == 'xfab':
        return np.array([[vec[0], vec[1], vec[2]],
                         [vec[1], vec[3], vec[4]],
                         [vec[2], vec[4], vec[5]]])
    
    elif input_format == 'mandel':
        sqr2 = np.sqrt(2)
        return np.array([[vec[0]     , vec[5]/sqr2, vec[4]/sqr2],
                         [vec[5]/sqr2, vec[1]     , vec[3]/sqr2],
                         [vec[4]/sqr2, vec[3]/sqr2, vec[2]     ]])
    
    else:
        return np.array([[vec[0]  , vec[5]/2, vec[4]/2],
                         [vec[5]/2, vec[1]  , vec[3]/2],
                         [vec[4]/2, vec[3]/2, vec[2]]])
    


    
def rotate_3x3_tensor(S, R, tol = 1e-6):
        """Return 3x3 matrix in rotated coordinate system
        
        Parameters
        -----------
        S (array_like) : 3x3 strain / stress matrix
        R (array_like) : 3x3 Rotation matrix
        tol (float)    : tolerance to consider R as a rotation matrix
        
        Returns
        -----------
        S' (array_like) : S in rotated coordinate system
        """
        # check R is a rotation matrix
        assert np.allclose(np.dot(R,R.T), np.eye(3)) and abs(np.linalg.det(R) - 1) < tol, 'R is not a rotation matrix'
        # return rotated coordinates
        return np.linalg.multi_dot([R,S,R.T])

    
    
    
def build_6x6_rot_mat(R, tol):
    """
    Return 6x6 transormation matrix corresponding to rotation R for a Voigt 6x6 stiffness tensor
    
    Parameters
    -----------
    R (array_like) : 3x3 rotation matrix
    tol (float)    : tolerance to consider R as a rotation matrix
    
    Returns
    -----------
    K (array_like) : 6x6 transformation matrix for the 6x6 stiffness tensor
    NOTE : K is not a rotation matrix, and unlike R it does not satisfy K.KT = I
    """
    
    assert np.allclose(np.dot(R,R.T), np.eye(3)) and abs(np.linalg.det(R) - 1) < tol, 'R is not a rotation matrix'
    
    q11,q12,q13,q21,q22,q23,q31,q32,q33 = R.flatten()
    
    k1 = np.array( [q11**2,    q12**2,    q13**2,    2*q12*q13,         2*q11*q13,         2*q11*q12] )
    k2 = np.array( [q21**2,    q22**2,    q23**2,    2*q22*q23,         2*q21*q23,         2*q21*q22] )
    k3 = np.array( [q31**2,    q32**2,    q33**2,    2*q32*q33,         2*q31*q33,         2*q31*q32] )
    k4 = np.array( [q21*q31,   q22*q32,   q23*q33,   q22*q33 + q23*q32, q21*q33 + q23*q31, q21*q32 + q22*q31] )
    k5 = np.array( [q11*q31,   q12*q32,   q13*q33,   q12*q33 + q13*q32, q11*q33 + q13*q31, q11*q32 + q12*q31] )
    k6 = np.array( [q11*q21,   q12*q22,   q13*q23,   q12*q23 + q13*q22, q11*q23 + q13*q21, q11*q22 + q12*q21] )
    
    return np.array([k1,k2,k3,k4,k5,k6])


def eigen_decomposition(T):
    """
    Compute the norm and orientation of principal components resulting from eigen decomposition of a strain / stress tensor

    Parameters
    ----------
    T (array_like) : 3x3 symmetric tensor

    Returns
    ----------
    eigvals_s (1x3 array) : Principal components values in decreasing order
    eigvecs_s (3x3 array) : Normalized principal component vectors by columns sorted accordingly to eigenvalues
    """
    eigenvalues, eigenvectors = np.linalg.eig(T)
    
    # Sort eigenvalues and eigenvectors in descending order
    idx = np.argsort(eigenvalues)
    eigvals_s = eigenvalues[idx][::-1]
    eigvecs_s = eigenvectors[:, idx[::-1]]
    
    return eigvals_s, eigvecs_s



def deviatoric(T, return_iso = False):
    """
    Decompose the strain / stress tensor into an isotropic and a deviatoric component : T = H + D, where H = tr(T) / 3
    
    Parameters
    ----------
    T (array_like)    : 3x3 symmetric tensor
    return_iso (bool) : if True, the isotropic component will be returned in addition to the deviatoric tensor. Default is False
    
    Returns
    ----------
    pi (float)    : Hydrostatic strain / stress. Hydrostatic pressure then equals to -Pi
    D (array_like): 3x3 deviatoric tensor.
    """
    pi = np.trace(T) / 3.0
    D = T - pi * np.identity(3)
    
    if return_iso:
        return pi, D
    else:
        return D


def invariants(T):
    """
    Returns the three invariants I1, I2, I3 of a 3x3 tensor
    """
    I1 = np.trace(T)
    I2 = 0.5 * (np.trace(np.dot(T, T)) - np.trace(T) ** 2)
    I3 = np.linalg.det(T)
    
    return I1, I2, I3


    
# Cij_symmetry for stiffness tensor (from Mouhat & Coudert 2014)
########################

Cij_symmetry = {
   'cubic':           np.array([[1, 7, 7, 0, 0, 0],
                                [7, 1, 7, 0, 0, 0],
                                [7, 7, 1, 0, 0, 0],
                                [0, 0, 0, 4, 0, 0],
                                [0, 0, 0, 0, 4, 0],
                                [0, 0, 0, 0, 0, 4]]),

   'trigonal_high':   np.array([[1, 7, 8, 9,  0, 0],
                                [7, 1, 8, -9, 0, 0],
                                [8, 8, 3, 0,  0, 0],
                                [9, -9, 0, 4, 0, 0],
                                [0,  0, 0, 0, 4, 9],
                                [0, 0, 0, 0,  9, 6]]),

   'trigonal_low':    np.array([[1,  7,  8,  9,  10,  0 ],
                                [7,  1,  8, -9, -10,  0 ],
                                [8,  8,  3,  0,   0,  0 ],
                                [9, -9,  0,  4,   0, -10],
                                [10,-10, 0,  0,   4,  9 ],
                                [0,  0,  0, -10 , 9,  6 ]]),

   'tetragonal_high': np.array([[1, 7, 8, 0, 0, 0],
                                [7, 1, 8, 0, 0, 0],
                                [8, 8, 3, 0, 0, 0],
                                [0, 0, 0, 4, 0, 0],
                                [0, 0, 0, 0, 4, 0],
                                [0, 0, 0, 0, 0, 6]]),

   'tetragonal_low':  np.array([[1, 7, 8, 0, 0, 11],
                                [7, 1, 8, 0, 0, -11],
                                [8, 8, 3, 0, 0, 0],
                                [0, 0, 0, 4, 0, 0],
                                [0, 0, 0, 0, 4, 0],
                                [11, -11, 0, 0, 0, 6]]),

   'orthorhombic':    np.array([[ 1,  7,  8,  0,  0,  0],
                                [ 7,  2, 12,  0,  0,  0],
                                [ 8, 12,  3,  0,  0,  0],
                                [ 0,  0,  0,  4,  0,  0],
                                [ 0,  0,  0,  0,  5,  0],
                                [ 0,  0,  0,  0,  0,  6]]),

   'monoclinic':      np.array([[ 1,  7,  8,  0,  10,  0],
                                [ 7,  2, 12,  0, 14,  0],
                                [ 8, 12,  3,  0, 17,  0],
                                [ 0,  0,  0,  4,  0,  20],
                                [10, 14, 17,  0,  5,  0],
                                [ 0,  0,  0, 20,  0,  6]]),

    'triclinic':       np.array([[ 1,  7,  8,  9,  10, 11],
                                 [ 7,  2, 12,  13, 14, 15],
                                 [ 8, 12,  3,  16, 17, 18],
                                 [ 9, 13, 16,  4,  19, 20],
                                 [10, 14, 17, 19,  5,  21],
                                 [11, 15, 18, 20,  21, 6 ]]),
   }


Cij_symmetry['hexagonal'] = Cij_symmetry['tetragonal_high']
Cij_symmetry['tetragonal'] = Cij_symmetry['tetragonal_high']
Cij_symmetry[None] = Cij_symmetry['triclinic']
