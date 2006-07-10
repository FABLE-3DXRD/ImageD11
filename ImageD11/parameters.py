



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

class parameters:
    def __init__(self,**kwds):
        """
        name=value style arg list
        """
        self.parameters = kwds

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
                print "setting: %s.%s from %s to %s"%(other,k,var,v)
                setattr(other,k,v)
            else:
                print "error: %s has no attribute %s, ignoring"%(other,k)

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
                self.parameters[name]=value
            except:
                print "Failed to read:",line
        self.dumbtypecheck()

    def dumbtypecheck(self):
        """
        Eventually parameter types (and units and fixed/varied to be
        specifieable
        For now it just tries to coerce to int, then float, then does nothing
        """
        for name, value in self.parameters.items():
            try:
                value = int(value)
                self.parameters[name]=value
            except:
                try:
                    value = float(value)
                    self.parameters[name]=value
                except:
                    # Hope it is a string?
                    pass
