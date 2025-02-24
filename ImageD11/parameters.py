# temporary comment out for new pars development:
# from xfab.parameters import *

from __future__ import print_function

import os
import sys

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
#
# Moved from ImageD11 to xfab 07/11/2019

"""
Class to handle groups of parameters to be saved in and out
of files and edited in guis with the fixed/varied info etc
"""


from xfab import xfab_logging

logger = xfab_logging.get_module_level_logger(__name__)


class AnalysisSchema(object):
    """Class to handle more complicated data analysis parameters."""

    def __init__(self, filename=None):
        # the dictionary of everything we found from the json file
        self.pars_dict = dict()

        # a parameters object for the geometry
        self.geometry_pars_obj = None

        # a dict of parameters object for each phase
        self.phase_pars_obj_dict = dict()
        
        # If a filename is supplied, load the filename and initialise the pars objects
        if filename is not None:
            self.load_json(filename)
            self.get_pars_objects()
    
    def load_json(self, filename):
        """
        Load json from file
        """
        import json
        with open(filename, 'r') as json_string:
            self.pars_dict = json.load(json_string)

        # store the path to the json file in 'parfile' key
        self.pars_dict['parfile'] = filename
        
    @classmethod
    def from_json(cls, filename):
        return cls(filename)
    
    def get_pars_objects(self):
        """Parses self.pars_dict, reads the .par files, makes parameter objects for them"""
        # get geometric parameters as an ImageD11 parameter object
        geometry_dict = self.pars_dict["geometry"]
        geometry_file = self.rel_to_absolute(geometry_dict["file"])
        geometry_pars_obj, _ = self.split_parameters_objs(parameters.from_file(geometry_file))
        self.geometry_pars_obj = geometry_pars_obj
        # remember the geometry.par filename
        self.geometry_pars_obj.set('filename', geometry_file)

        # get the phases entry from the json file
        phase_filenames_dict = self.pars_dict["phases"]
        for phase_name, phase_entry in phase_filenames_dict.items():
            # get the absolute phase.par filename
            phase_file = self.rel_to_absolute(phase_entry["file"])
            # read the phase pars from disk
            phase_pars_obj = parameters.from_file(phase_file)
            # add parameters to keep track of the name of the phase and its file path
            phase_pars_obj.set('filename', phase_file)
            # add the phase
            self.add_phase_from_pars_obj(phase_name, phase_pars_obj)
    
    @staticmethod
    def split_parameters_dicts(pars_dict):
        """
        Split pars dict parameters into (geometry, phase) dicts
        """
        # split the cell parameters from the geometry parameters
        geometry_pars_dict = {key:value for key,value in pars_dict.items() if not (('cell' in key) or ('phase_name' in key))}
        phase_pars_dict = {key:value for key,value in pars_dict.items() if (('cell' in key) or ('phase_name' in key))}
        return geometry_pars_dict, phase_pars_dict
    
    @staticmethod
    def split_parameters_objs(pars_obj):
        """
        Split pars object into (geometry, phase) objects
        """
        pars_dict = dict(pars_obj.get_parameters())
        geometry_pars_dict, phase_pars_dict = AnalysisSchema.split_parameters_dicts(pars_dict)
        geometry_pars_obj = parameters.from_dict(geometry_pars_dict)
        phase_pars_obj = parameters.from_dict(phase_pars_dict)
        return geometry_pars_obj, phase_pars_obj
        
    def add_phase_from_pars_obj(self, phase_name, phase_pars_obj):
        """
        Add phase from ImageD11 pars object
        """
        # the base case
        # need to strip geometry pars out
        _, phase_pars_obj = self.split_parameters_objs(phase_pars_obj)
        phase_pars_obj.set('phase_name', phase_name)
        # put this pars object in self.phase_pars_obj_dict
        self.phase_pars_obj_dict[phase_name] = phase_pars_obj
    
    def add_phase_from_dict(self, phase_name, phase_dict):
        """Add phase from pars dict"""
        phase_obj = parameters.from_dict(phase_dict)
        self.add_phase_from_pars_obj(phase_name, phase_obj)
        
    def add_phase_from_pars_file(self, phase_name, phase_pars_file):
        """
        Add phase from ImageD11 pars file
        """
        # get pars object from ImageD11 pars file
        phase_obj = parameters.from_file(phase_pars_file)
        self.add_phase_from_pars_obj(phase_name, phase_obj)
    
    def add_phase_from_unitcell(self, phase_name, unitcell):
        """
        Add phase from unitcell object
        
        phase_name: the name of the phase
        unitcell: the unitcell object
        """
        # get parameter object from unitcell
        phase_pars_obj = unitcell.to_par_obj()
        self.add_phase_from_pars_obj(phase_name, phase_pars_obj)

    def rel_to_absolute(self, rel_file):
        """Get absolute path from (could be relative) filename"""
        parent_folder = os.path.dirname(os.path.abspath(self.pars_dict['parfile']))
        return os.path.join(parent_folder, rel_file)
    
    def update_parfiles(self):
        # updates geometry and phase par files on disk
        self.update_geometryfile()
        for phase_name in self.phase_pars_obj_dict.keys():
            self.update_phasefile(phase_name)
    
    def update_geometryfile(self):
        """Write the geometry.par file to disk"""
        if self.geometry_pars_obj is None:
            raise ValueError('No pars object in self.geometry_pars_obj!')
        self.geometry_pars_obj.saveparameters(self.geometry_pars_obj.get('filename'))

    def update_phasefile(self, phase_name, guess_filename=False):
        """Write the phase with phase_name to disk.
           guess_filename: If no filename present for this phase, try to guess it from the geometry filename"""
        # update phase file on disk from memory
        try:
            phase_pars = self.phase_pars_obj_dict[phase_name]
        except KeyError:
             raise KeyError("Phase " + phase_name + " not found")
        
        # find out where to put the phase file
        try:
            filename = phase_pars.get('filename')
        except KeyError:
            if guess_filename:
                # guess the filename from where we put geometry.par
                dest_folder = os.path.dirname(self.geometry_pars_obj.get('filename'))
                filename = os.path.join(dest_folder, phase_name + ".par")
                # save it to the phase file
                self.phase_pars_obj_dict[phase_name].set('filename', filename)
            else:
                raise ValueError("Don't know where to save this phase! Try setting guess_filename=True")
        
        self.phase_pars_obj_dict[phase_name].saveparameters(filename)
    
    def get_any_phase_pars_obj(self):
        """Returns any parameters object from self.phase_pars_obj_dict"""
        return next(iter(self.phase_pars_obj_dict.values()))
    
    def to_old_pars_dict(self, phase_name=None):
        """Produce an old-style ImageD11 parameters dict"""
        # get geometry pars as a dict
        geometry_pars_dict = self.geometry_pars_obj.get_parameters()
        if phase_name is not None:
            # get parameters for a specific phase
            phase_pars_dict = self.phase_pars_obj_dict[phase_name].get_parameters()
            # copy the geometry dict to a new dict
            pars_dict = geometry_pars_dict.copy()
            # add in the phase pars
            pars_dict.update(phase_pars_dict)
            # re-add the geometry dict to override things like filename if present
            pars_dict.update(geometry_pars_dict)
            if 'phase_name' in pars_dict.keys():
                del pars_dict['phase_name']
        else:
            # just copy the geometry dict
            pars_dict = geometry_pars_dict.copy()
        if 'filename' in pars_dict.keys():
            del pars_dict['filename']
        return pars_dict
    
    def to_old_pars_object(self, phase_name=None):
        """Produce an old-style ImageD11 parameters object"""
        pars_dict = self.to_old_pars_dict(phase_name=phase_name)
        pars_object = parameters.from_dict(pars_dict)
        return pars_object
    
    def to_old_pars_file(self, filename, phase_name=None):
        """Write an old-style ImageD11 .par file"""
        pars_object = self.to_old_pars_object(phase_name=phase_name)
        pars_object.saveparameters(filename)
    
    def get_xfab_pars_dict(self, phase_name=None):
        """Produce an old-style parameters dict"""
        return self.to_old_pars_dict(phase_name=phase_name)
    
    def update_pars_dict(self):
        """Update pars_dict from current geometry and phases"""
        if ('phases' not in self.pars_dict.keys()) and (len(self.phase_pars_obj_dict) > 0):
            self.pars_dict['phases'] = dict()
        for phase_name, phase_pars_obj in self.phase_pars_obj_dict.items():
            if phase_name not in self.pars_dict['phases'].keys():
                self.pars_dict['phases'][phase_name] = {}
            if 'filename' not in phase_pars_obj.get_parameters().keys():
                dest_folder = os.path.dirname(self.geometry_pars_obj.get('filename'))
                filename = os.path.join(dest_folder, phase_name + ".par")
                # save it to the phase file
                self.phase_pars_obj_dict[phase_name].set('filename', filename)
            self.pars_dict['phases'][phase_name]['file'] = phase_pars_obj.get('filename')
    
    def save(self, filename=None):
        """
        Save current state - both to json and to individual .par files. Looks for filename in self.pars_dict['parfile'], or can be supplied.
        """
        if filename is None:
            try:
                filename = self.pars_dict['parfile']
            except KeyError:
                raise ValueError("Don't know where to save to! Supply a filename")
        
        self.update_pars_dict()
        
        import json
        json_object = json.dumps(self.pars_dict, indent=2)
        with open(filename, 'w') as json_file:
            json_file.write(json_object)
        
        # write the geometry file
        self.update_geometryfile()
        
        # write the phases
        for phase_name in self.phase_pars_obj_dict.keys():
            self.update_phasefile(phase_name, guess_filename=True)
    
    @classmethod
    def from_geom_and_phase_dict(cls, geom_dict, phase_dict):
        """Create a new AnalysisSchema object from a geometry dict and a phase dict"""
        # set up empty schema object
        schema_obj = cls()
        # see if we have a filename for the geometry.par file
        geom_file = geom_dict.get('filename', 'geometry.par')
        gd = {
            "file": geom_file
        }
        schema_obj.pars_dict["geometry"] = gd
        schema_obj.geometry_pars_obj = parameters.from_dict(geom_dict)
        schema_obj.geometry_pars_obj.set('filename', geom_file)
        
        # set up the phases
        # add parameters to keep track of the name of the phase and its file path
        phase_name = phase_dict.get('phase_name', 'phase')
        phase_file = phase_dict.get('filename', phase_name + '.par')
        phase_dict['filename'] = phase_file
        schema_obj.add_phase_from_dict(phase_name, phase_dict)
        
        return schema_obj
    
    @classmethod
    def from_old_pars_dict(cls, pars_dict):
        """Create a new AnalysisSchema object from an old parameters dict"""
        geometry_pars_dict, phase_pars_dict = cls.split_parameters_dicts(pars_dict)
        print(phase_pars_dict)
        return cls.from_geom_and_phase_dict(geometry_pars_dict, phase_pars_dict)
    
    @classmethod
    def from_old_pars_object(cls, pars_obj):
        """Create a new AnalysisSchema object from an old parameters object"""
        pars_dict = pars_obj.get_parameters()
        return cls.from_old_pars_dict(pars_dict)
    
    @classmethod
    def from_old_pars_file(cls, filename):
        """Create a new AnalysisSchema object from an old parameters file"""
        pars_obj = parameters.from_file(filename)
        return cls.from_old_pars_object(pars_obj)
    
    @classmethod
    def from_default(cls, detector='eiger'):
        """Load default detector parameters from disk for either 'eiger' or 'frelon' detecor"""
        # get path to either eiger or frelon default geometric parameters
        # import pkg_resources
        # geom_par_path = pkg_resources.resource_filename("ImageD11","data/{det}_example_geometry.par".format(det=detector))
        # phase_par_path = pkg_resources.resource_filename("ImageD11","data/CeO2.par")
        # geom_par_path = os.path.join(os.path.dirname(sys.modules['ImageD11'].__file__), '..', 'data', '{det}_example_geometry.par'.format(det=detector))
        # phase_par_path = os.path.join(os.path.dirname(sys.modules['ImageD11'].__file__), '..', 'data', 'CeO2.par')
        # geom_obj = parameters.from_file(geom_par_path)
        # phase_obj = parameters.from_file(phase_par_path)
        
        if detector == 'eiger':
            geom_dict = {'chi': 0.0,
                         'distance': 152736.55305695778,
                         'fit_tolerance': 0.05,
                         'min_bin_prob': 1e-05,
                         'no_bins': 10000,
                         'o11': -1,
                         'o12': 0,
                         'o21': 0,
                         'o22': -1,
                         'omegasign': 1.0,
                         't_x': 0,
                         't_y': 0,
                         't_z': 0,
                         'tilt_x': 0.0,
                         'tilt_y': 0.0,
                         'tilt_z': 0.0,
                         'wavelength': 0.2845704,
                         'wedge': 0.0,
                         'weight_hist_intensities': 0,
                         'y_center': 1049.9295061162281,
                         'y_size': 75.0,
                         'z_center': 1116.4472483389864,
                         'z_size': 75.0}
        elif detector == 'frelon':
            geom_dict = {'chi': 0.0,
                         'distance': 136062.54166078364,
                         'fit_tolerance': 0.05,
                         'min_bin_prob': 1e-05,
                         'no_bins': 10000,
                         'o11': 1,
                         'o12': 0,
                         'o21': 0,
                         'o22': -1,
                         'omegasign': 1.0,
                         't_x': 0.0,
                         't_y': 0.0,
                         't_z': 0.0,
                         'tilt_x': 0.0,
                         'tilt_y': 0.0,
                         'tilt_z': 0.0,
                         'wavelength': 0.28457041,
                         'wedge': 0.0,
                         'weight_hist_intensities': 0,
                         'y_center': 1081.849909211387,
                         'y_size': 47.0,
                         'z_center': 1017.0452628833956,
                         'z_size': 47.0}
        else:
            raise ValueError('Invalid detector! Options are frelon, eiger')
        
        phase_dict = {'cell__a': 5.41143,
                     'cell__b': 5.41143,
                     'cell__c': 5.41143,
                     'cell_alpha': 90.0,
                     'cell_beta': 90.0,
                     'cell_gamma': 90.0,
                     'cell_lattice_[P,A,B,C,I,F,R]': 225}

        phase_obj = parameters.from_dict(phase_dict)
        geom_obj = parameters.from_dict(geom_dict)
        
        phase_obj.set('phase_name', 'CeO2')
        phase_obj.set('filename', 'CeO2.par')
        asc = cls.from_geom_and_phase_dict(geom_obj.get_parameters(), phase_obj.get_parameters())
        return asc

class par:
    """
    Represents a thing which can vary
    """

    def __init__(self, name, value, helpstring=None,
                 vary=False, can_vary=False, stepsize=None):
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
            self.helpstring = "parameter : " + name
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
        return [self.name,
                self.value,
                self.helpstring,
                self.vary,
                self.can_vary,
                self.stepsize]


class parameters:
    """
    Class to hold a set of named parameters
    """

    def __init__(self, **kwds):
        """
        name=value style arg list
        """
        self.parameters = kwds
        self.varylist = []
        self.can_vary = {}
        self.variable_list = []
        self.stepsizes = {}
        self.par_objs = {}
        for k, v in list(self.parameters.items()):
            self.addpar(par(k, v))

    def addpar(self, par):
        """
        add a parameter object
        """
        self.parameters[par.name] = par.value
        self.can_vary[par.name] = par.can_vary
        if par.vary and par.name not in self.varylist:
            self.varylist.append(par.name)
        if par.can_vary and par.name not in self.variable_list:
            self.variable_list.append(par.name)
            self.stepsizes[par.name] = par.stepsize
        self.par_objs[par.name] = par

    def get_variable_list(self):
        return self.variable_list

    def get_variable_values(self):
        """ values of the parameters """
        return [self.parameters[name] for name in self.varylist]

    def get_variable_stepsizes(self):
        """ stepsizes for optimisers """
        return [self.stepsizes[name] for name in self.varylist]

    def set_varylist(self, vl):
        ks = list(self.parameters.keys())
        for v in vl:
            assert v in ks
            assert v in self.variable_list
        self.varylist = vl

    def set_variable_values(self, values):
        """ set values of the parameters"""
        assert len(values) == len(self.varylist)
        for name, value in zip(self.varylist, values):
            self.parameters[name] = value

    def set_parameters(self, d):
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

    def get(self, name):
        return self.parameters[name]

    def set(self, name, value):
        self.parameters[name] = value

    def update_yourself(self, other):
        """
        Sychronise this parameter objects list of values with another object
        """
        for k, v in list(self.parameters.items()):
            if hasattr(other, k):
                var = getattr(other, k)
                logger.debug("setting: pars[%s] from %s to %s" % (k, v, var))
                self.parameters[k] = var
            else:
                logger.debug("error: %s has no attribute %s, ignoring" % (other, k))

    def update_other(self, other):
        """
        Synchronise an object with the values in this object
        """
        for k, v in list(self.parameters.items()):
            if hasattr(other, k):
                var = getattr(other, k)
                logger.debug("setting: %s.%s from %s to %s" % (other, k, var, v))
                setattr(other, k, v)
            else:
                logger.debug("error: %s has no attribute %s, ignoring" %
                             (other, k))

    def saveparameters(self, filename):
        """
        Write parameters to a file
        """
        with open(filename, "w") as f:
            keys = list(self.parameters.keys())
            keys.sort()
            for key in keys:
                f.write("%s %s\n" % (key, str(self.parameters[key])))

    def loadparameters(self, filename, phase_name=None):
        """
        Load parameters from a file
        """
        # is it a json file?
        if filename.endswith('json'):
            # create a JsonPars object and get the pars dict from it
            pars_dict = AnalysisSchema(filename=filename).get_xfab_pars_dict(phase_name)
            self.parameters.update(pars_dict)
        else:
            with open(filename, "r") as f:
                lines = f.readlines()
            for line in lines:
                try:
                    [name, value] = line.split(" ")
                    name = name.replace("-", "_")
                    self.parameters[name] = value
                except ValueError:
                    logger.error("Failed to read:%s" % (line))
        self.dumbtypecheck()

    def dumbtypecheck(self):
        """
        Eventually parameter types (and units and fixed/varied to be
        specifieable
        For now it just tries to coerce to float, then does nothing
        """
        for name, value in list(self.parameters.items()):
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
    
    @classmethod
    def from_dict(cls, pars_dict):
        p = parameters()
        p.set_parameters(pars_dict)
        return p
        
    @classmethod
    def from_file(cls, filename, phase_name=None):
        return read_par_file(filename, phase_name=None)


def read_par_file(filename, phase_name=None):
    p = parameters()
    p.loadparameters(filename, phase_name=phase_name)
    return p

