
"""
Module(s) for scaling peaks

- Read in a series of columnfile type objects
    collection of columnfiles

- Read in some symmetry information: 
    there are 32 crystallographic point groups
    11 contain a centre of symmetry and are also laue groups 

- For each peak 
    assign it a unique "key" in a list of reflections
    typically h>=k>=l of the symmetry equivalents
    Give it a scale factor (or two)
    Compute the corrected intensity and derivative w.r.t scale factors

- For each uniq key compute the average group intensity
    with or without median style filtering

- Compute the overall R(int) merging statistic and minimise that w.r.t scaling 
    variables

- Optimise scale factor variables
"""

import columnfile

class data_table(object):
    """
    Wraps a columnfile/hklf4/brukerraw file
    Gives h,k,l
      -    variables [scales per frame, scan, grain etc]
      -    derivatives [of all variables]
    To constrain grain scale factors across different scans imply
    we need to have grains in same datafile.
    """
    def __init__(self):
        """
        Initialise some variables
        """
        self.varnames = [ "scale" ]
        self.varvals  = [ 1.0 ]
        self.varlimits = [ (0.0, None) ]
        self.npks = 0
    def get_hkls(self):
        """
        Return the hkl indices of the peaks
        """
        raise Exception("Write me")
    def set_keys(self, keys):
        """
        Store a list of uniq keys which go with each peak
        These are used to group peaks into symmetry equivalents
        """
        self.keys = keys
    def calc(self):
        pass
        

class scalem(object):
    """
    This class holds a series of tables (assumes everything fits in memory)
    ... and does some fitting
    """
    def __init__(self):
        """
        Holds a list of unmerged data tables
        These correspond to a single scans containing N grains

        """
        self.data_objects = []
    def add_data(self, filename, format="ImageD11"):
        """
        Read a columnfile object. Allowed to have format == ImageD11|bruker|hklf4
        """
        assert fo
        if format == "ImageD11":
            c = columnfile.columnfile(filename)
            assert "h" in c.titles
            assert "k" in c.titles
            assert "l" in c.titles
            # Need to compute the Intensity using the Lorentz factor
            assert "Iobs" in c.titles
            # Need to know which grain is which in multigrain data...
            # etc ....
        if format == "bruker":
            # Columns are from a bruker raw file
            raise Exception("Write that module")
        if format == "hklf4":
            # Data are 3I4,2F,I
            # h k l F2 sF2 ibatch
            raise Exception("Write that module")
        # Per columnfile we have an overall scale factor
        # Also scale factors as a function of some other columns
        # These need to be generalised
    def set_sym(self, sym):
        """
        Set the symmetry to sym
         - list of matrices (3x3)
         - list of operators (4x4)
         - string (eg m3m etc)
        """
        raise Exception("Write me")
    def calc(self):
        """
        Compute for each data object:
            Columns: I_corr, sigma_corr
            Foreach scaling variable: dI_corr/dvar, dsigma_corr/dvar
        """
        pass
        # for o in self.





# Sheldrick Sienna tutorial coding
#
# 
#
# h,k,l integer arrays (n reflections)
# ip,iq (?)
# FF,SI - F*F, sig_F
# SY - 24 symmetry matrices
#
# Fill in symmetry matrices
#
# LOOP over -h, h
#    LOOP over symops
#        pick one
# read next
# == generate key
#
# Call INSORT routine. Sorts reflections. iP, iQ used here
# Sorts on h, then k, then l
# == sort on key
#
# Reflections are now in group order (sorted on key)
#
# Form R(int) and check for absences and centrics
#
# Done.




def generate_key(hkl,symmops):
    """ To be redone in C/fortran etc 
    Expect identity to be in symmops
    All symmetry matrices, not including inversion
    """
    
    res = np.zeros( (3, hkl.shape[1], len(symmops)*2) , np.int16 )
    j = 0
    for s in [-1,1]:
        for op in symmops:
            res[:,j] = np.round(s*np.dot( op, hkl )).astype(np.int16)

            

            









