
import os, logging

# Try to manage the options coming and going to the different scripts
#
# Load / save / log command lines to reproduce options set when driving scripts
# Offer gui helpers to launch programs
#
# Previously we had "optparse", this becomes "argparse" for later pythons
#
# The idea is to use something like gooey or argparseui which hook to 
# the parsers _actions list to see what can be offered.
#
# We define a few "types" corresponding to ImageD11 known filetypes to 
# help when identifying what an optiong is so we could hook up plotting
# and editing to the different entities

class FileType(object):
    def __init__(self, mode='r'):
        assert mode in 'rw'
        self._mode = mode
    def __call__(self, string):
        if  'r' in self._mode:
            if not os.path.exists( string ):
                logging.error("File %s not found", string)
        return string

class ParameterFileType(FileType):
    pass

class ColumnFileType(FileType):
    pass

class UbiFileType(FileType):
    pass

class ImageFileType(FileType):
    pass

class SplineFileType(FileType):
    pass

class GvectorFileType(FileType):
    pass

class HdfFileType(FileType):
    pass


