



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
Class to handle groups of parameters to be saved in and out
of files and edited in guis with the fixed/varied info etc
"""

class par:
    """
    Represents a thing which can vary
    """
    def __init__(self, name, value, helpstring = None ,
                 vary=False, can_vary=False , stepsize = None):
        """
        name : unique key used as keyword arg to some functions
        value : value of the parameter
        helpstring : optional string to help user
        vary : value should be optimised
        can_vary : value is not fixed
        stepsize : guessestimated ball park step size (eg not 1e99!)
        """
        self.name = name
        self.value = value
        if helpstring is None:
            self.helpstring = "parameter : "+name
        else:
            self.helpstring = helpstring
        self.vary = vary
        self.can_vary = can_vary
        self.stepsize = stepsize
        
    def fromstringlist(self, sl):
        """ to send to Java """
        [self.name,
         self.value,
         self.helpstring,
         self.vary,
         self.can_vary,
         self.stepsize] = sl
        
    def tostringlist(self):
        """ to catch from Java """
        return [self.name ,
                self.value ,
                self.helpstring ,
                self.vary ,
                self.can_vary,
                self.stepsize ]

import logging 
class parameters:
    """
    Class to hold a set of named parameters
    """
    def __init__(self,**kwds):
        """
        name=value style arg list
        """
        self.parameters = kwds
        self.varylist = []
        self.can_vary = {}
        self.variable_list = []
        self.stepsizes = {}
        self.par_objs = {}
        for k,v in self.parameters.items():
            self.addpar(par(k,v))
            
    def addpar(self,par):
        """
        add a parameter object
        """
        self.parameters[par.name] = par.value
        self.can_vary[par.name] = par.can_vary
        if par.vary and par.name not in self.varylist:
            self.varylist.append(par.name)
        if par.can_vary and par.name not in self.variable_list:
            self.variable_list.append(par.name)
            self.stepsizes[par.name]=par.stepsize
        self.par_objs[par.name]=par

    def get_variable_list(self):
        return self.variable_list

    def get_variable_values(self):
        """ values of the parameters """
        return [self.parameters[name] for name in self.varylist]

    def get_variable_stepsizes(self):
        """ stepsizes for optimisers """
        return [self.stepsizes[name] for name in self.varylist]

    def set_varylist(self, vl):
        ks = self.parameters.keys()
        for v in vl:
            assert v in ks
            assert v in self.variable_list
        self.varylist = vl

    def set_variable_values(self,values):
        """ set values of the parameters"""
        assert len(values)==len(self.varylist)
        for name, value in zip(self.varylist,values):
            self.parameters[name]=value

    def set_parameters(self,d):
        """
        Updates the values of parameters
        """
        self.parameters.update(d)
        self.dumbtypecheck()

    def get_parameters(self):
        """
        Returns a dictionary of parameters
        """
        return self.parameters

    def get(self,name):
        return self.parameters[name]
    
    def set(self,name,value):
        self.parameters[name] = value

    def update_yourself(self,other):
        """
        Sychronise this parameter objects list of values with another object
        """
        for k,v in self.parameters.items():
            if hasattr(other,k):
                var = getattr(other,k)
                print "setting: pars[%s] from %s to %s"%(k,v,var)
                self.parameters[k]=var
            else:
                print "error: %s has no attribute %s, ignoring"%(other,k)

    def update_other(self,other):
        """
        Synchronise an object with the values in this object
        """
        for k,v in self.parameters.items():
            if hasattr(other,k):
                var = getattr(other,k)
                logging.debug("setting: %s.%s from %s to %s"%(other,k,var,v))
                setattr(other,k,v)
            else:
                logging.debug("error: %s has no attribute %s, ignoring"%
                              (other,k))

    def saveparameters(self,filename):
        """
        Write parameters to a file
        """
        f=open(filename,"w")
        keys=self.parameters.keys()
        keys.sort()
        for key in keys:
            f.write("%s %s\n"%(key,str(self.parameters[key])))
        f.close()

    def loadparameters(self,filename):
        """
        Load parameters from a file
        """
        lines = open(filename,"r").readlines()
        for line in lines:
            try:
                [name, value] = line.split(" ") 
                name=name.replace("-","_")
                self.parameters[name]=value
            except ValueError:
                print "Failed to read:",line
        self.dumbtypecheck()

    def dumbtypecheck(self):
        """
        Eventually parameter types (and units and fixed/varied to be
        specifieable
        For now it just tries to coerce to float, then does nothing
        """
        for name, value in self.parameters.items():
            if type(value) == type("string"):
                try:
                    vf = float(value)
                except ValueError:
                    # it really is a string
                    self.parameters[name] = value.lstrip().rstrip()
                    continue
                # here if float worked
                try:
                    vi = int(value)
                except ValueError:
                    # it really is a float
                    self.parameters[name] = vf
                    continue
                    
                # here if float and int worked
                # should not be needed, depends on int valueerror
                if abs(vi - vf) < 1e-9:
                    # use int
                    self.parameters[name] = vi
                    continue
                else:
                    self.parameters[name] = vf
                    continue
            else:
                # int/float preserve type
                self.parameters[name] = value

def read_par_file( filename ):
    p = parameters()
    p.loadparameters( filename )
    return p
