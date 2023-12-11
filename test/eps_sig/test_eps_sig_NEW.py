
from __future__ import print_function, division


import unittest, os
import numpy as np
from ImageD11 import eps_sig_solver_NEW as ess
import ImageD11.parameters

def p(x):
    return os.path.join( os.path.dirname( __file__ ), x )

def test_create(parfile=p('mypar.par'), grainfile=p('CuAlBe_scan10.map')):    
    pars = ImageD11.parameters.read_par_file(parfile)
    ubis = [ g.ubi for g in ImageD11.grain.read_grain_file( grainfile ) ]
    unitcell = [pars.get('cell__' + a) for a in 'abc'] +\
               [pars.get('cell_' + a) for a in ('alpha','beta','gamma')]
    args = { 'unitcell' : unitcell, 'UBI_list' : ubis }
    
    for newname, oldname in [('name', 'name'),
                             ('stress_unit','stress_unit'),
                             ('symmetry', 'crystal_symmetry')]:
        if newname in pars.parameters:
            args[newname] = pars.get(newname)
        elif oldname in pars.parameters:
            args[newname] = pars.get(oldname)
    for i in range(1,7):
        for j in range(i,7):
            name = 'c%d%d'%(i,j)
            val = pars.get( name )
            if val == 'None':
                val = None
            args[name] = val
    # todo : finite strain should accept grains to get ubi/ub
    # todo : eps_sig_solver should be able to get loaded from a parfile
    e = ess.EpsSigSolver( **args )
    assert e is not None
    


