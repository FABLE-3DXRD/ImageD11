

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
Interface between Tkinter gui and the actual useful code.

There should be no scientific algorithms (eventually) on
the gui side of this class.

This class will eventually offer macro recording capability.
"""

import logging, sys
# Things to offer from gui
from ImageD11 import peakmerge, indexing, transformer

# import imaged11 # breaks the code for esrf/fable


# To autoconvert arrays to lists for Java XMLRPC
RETURN_NUMERICS = False
import Numeric
TYPE_NUMERIC = type(Numeric.zeros(1)) 

class guicommand:
    """
    Keeps a log of all commands issued - separates gui code from
    algorithmical code
    """
    def __init__(self):
        self.objects = { "peakmerger" : peakmerge.peakmerger(),
                         "transformer": transformer.transformer(),
                         "indexer"    : indexing.indexer()
                         }

        self.commandscript = """# Create objects to manipulate - they hold your data
#
from ImageD11 import peakmerge, indexing, transformer
mypeakmerger = peakmerge.peakmerger()
mytransformer = transformer.transformer()
myindexer = indexing.indexer()
#
# Your work starts here:
#    
"""

    def execute(self,object,command,*args,**kwds):
        """
        Pass in object as string [peakmerger|transformer|indexer]
        Pass in command as string, getattr(command) will be used
        Returns the return value of the function....

        TODO : change this interface???
             eg : works - returns True
                          you look for self.lastreturned
                  fails - returns False
                          you look for self.lasttraceback
        """
        if object not in self.objects.keys():
            raise Exception("ERROR! Unknown command object")
        o = self.objects[object]
        ran = "my%s.%s("%(object,command)
        if command.find("."):
            subobjs = command.split(".")[:-1]
            for s in subobjs:
                o = getattr(o,s)
            command = command.split(".")[-1]
        func = getattr(o,command)
        try:
            addedcomma = ""
            for a in args:
                ran="%s %s %s"%(ran,addedcomma,repr(a))
                addedcomma=","
            for k,v in kwds.items():
                ran="%s %s %s=%s "%(ran,addedcomma,k,v)
                addedcomma=","
            ran+=" )\n"
            print "Running:",ran,
            sys.stdout.flush()
            ret = func(*args, **kwds)

        except:
            print self
            print object
            print command
            print func
            print args
            print kwds
            print "Exception occured"
            import traceback
            traceback.print_exc()
            return "Exception occured in the python " + ran
        print " OK!"
        self.commandscript+=ran
        return ret

    def getdata(self,object,name):
        """
        Allows access to "live" data in the objects wrapped

        By passing references back you can circumvent the
        cleanliness of the interface. Please dont.

        Returns object.name
        """
        if object not in self.objects.keys():
            raise Exception("ERROR! Unknown command object")
        attribute = getattr(self.objects[object],name)
        if RETURN_NUMERICS:
            # Normally python will get this
            print "python return array"
            return attribute
        if type(attribute) == TYPE_NUMERIC:
            # Java gets this for arrays
            print "Java return list"
            return attribute.tolist()
        return attribute

    def gethistory(self):
        return self.commandscript
