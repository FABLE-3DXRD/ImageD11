
from __future__ import print_function, division
import collections, json

"""
Load and save a data analysis project

   Images part - all the frames in one long table
               - which frames are adjacent
   We use a namedtuple for each scan point and record all the things
   that we care about as a tuple: 
          (filename, Omega, dty, z, Load, etc) <- can be hashed
   The edges tell us which frames want/need merging.

   Processing parts:
   foreach node:
       Corrected images (dark, flat, ... )
       Background image(s)
       Threshold estimation
       Label images
       List of peaks

"""


RotationNode = collections.namedtuple( "RotationNode",
                                       "filename Omega" )

DiffTomoNode = collections.namedtuple( "DiffTomoNode",
                                       "filename Omega dty" )

class Scan( object ):
    """ Scan data as a graph
    nodes = list of datapoints 
    edges = which nodes are adjacent to this node
    """
    def __init__(self,
                 nodes = None,    # a List[]
                 edges = None,    # a List[[int,]]
                 ):      
        """
        nodes = measurement points in the scan
        edges = edges[i] = connections to point i
        hidden self.node_index that gives you position in scan
        """
        if nodes is None:
            self.nodes = []
        else:
            self.nodes = nodes
        if edges is None:
            self.edges = []
        else:
            self.edges = edges
        # Given a node, node_index finds it in the list
        self.node_index = {}
        for i, node in enumerate(self.nodes):
            self.node_index[node] = i
        if __debug__:
            for i, e in enumerate(self.edges): # i is implicit
                for j in e:
                    assert j>=0 and i<len(nodes)

    def addNode( self, node, neighbors = [] ):
        """
        Adding another node
         node = node to be added
         neighbors = images to link up to
                     usually a list of negative integers
                     e.g. [ -1,] = previous in 1D
                          [ -1, -180] = previous in omega, dty
                     Can be empty
                     Can be for interlaced, etc
        Always stored as positive
        """
        adr  = len(self.nodes)          # adr = address in array
        self.nodes.append( node )
        self.node_index[node] = adr
        nbs = []                        # nbs = neighbors
        err = 0                         # bad refs given
        for inode in neighbors:
            if inode > 0 and inode < adr:
                neighbor = inode
            elif inode < 0 and inode >= -adr:
                neighbor = inode + adr
            else:
                err +=1 # Error
                continue
            self.edges[ neighbor ].append( adr )
            nbs.append( neighbor )
        self.edges.append( nbs )
        assert len(self.edges) == len(self.nodes)
        return err

    def neighbors( self, node ):
        """
        node = scan point (Image, Omega, etc)
        node = integer
        """
        if isinstance( node, int ):
            return self.edges[node]
        if isinstance( node, tuple): 
            i = self.node_index[ node ]
            return self.edges[i]
        raise TypeError("node not understood, want tuple or int row index")

    def __len__(self):
        return len(self.nodes)

    def __getitem__(self, i):
        """ We index to the node, not the edges ? """
        return self.nodes[i]

    def todict( self ):
        """
        Convert to a dictionary representation for serialisation
        """
        data = { "titles": self.nodes[0]._fields,
                 "rows" : [[field for field in node ]
                                 for node in self.nodes],
                 "edges" : self.edges } 
        return data 

    def fromdict( self, dct ):
        """
        Convert from a dictionary representation for serialisation
        """
        self.edges = dct['edges']
        # named tuple
        tup = collections.namedtuple( "Node", dct['titles'] )
        for i, item in enumerate( dct['rows'] ):
            t = tup( *item ) 
            self.nodes.append( t )
            self.node_index[t] = i
    

def add_motor( inputscan, motorname, motorpos ):
    """
    Adds another motor to a scan
       e.g. dty position when not in header
    """
    fields = inputscan[0]._fields + (motorname,)
    tup = collections.namedtuple( "Node", fields )
    try:
        if len( inputscan ) == len( motorpos ):
            nodes = [ tup(*(vals + (mot,))) for vals, mot in
                      zip(inputscan, motorpos)]
        else:
            raise Exception("Motorpos has wrong length")
    except TypeError: # no len
        mot = float( motorpos )
        nodes = [ tup( *(vals + (mot,))) for vals in inputscan ]
    edges = [e for e in inputscan.edges]
    return Scan( nodes, edges )


def fablescandata_from_spt( spt, headeritems=["Omega",],
                            nodeclass = RotationNode ):
    """
    Read an ImageD11 peaksearch spt file to return a Scan object
    """
    myscan = Scan()
    # Unpack the [motor|counter]_mne/_pos pairs
    motor_mne = counter_mne = []
    with open(spt, "r") as sptfile:
        for line in sptfile.readlines():
            if line[0] != "#": # skip blank and peaks
                continue
            #             0123456
            if line.find("# File ")==0:
                filename = line[7:-1] # trim newline
                header = {}
                continue
            if line.find("# motor_mne") == 0:
                motor_mne = line.split("=")[1].split()
                continue
            if line.find("# motor_pos") == 0:
                motor_pos = [float(v) for v in line.split("=")[1].split()]
                for k,v in zip( motor_mne, motor_pos ):
                    header[ k ] = v
                continue
            if line.find("# counter_mne") == 0:
                counter_mne = line.split("=")[1].split()
                continue
            if line.find("# counter_pos") == 0:
                counter_pos = [float(v) for v in line.split("=")[1].split()]
                for k,v in zip( counter_mne, counter_pos ):
                    header[ k ] = v
                continue
            if line.find("# npks") == 0:
                args = { key:float(header[key]) for key in headeritems }
                mynode = nodeclass( filename=filename, **args )
                myscan.addNode( mynode, (-1,) )
                continue
            for hitem in headeritems:
                if line.find(hitem) > 0:
                    val = line[1:].split( "=" )[-1].strip()
                    header[ hitem ] = val
            # end of frame marker
    return myscan


def sinogram_from_spt_list( sptlist ):
    scans = [fablescandata_from_spt( sptfile ) for sptfile in sptlist ]


def mergeScans( scan1, scan2 ):
    """
    Merge two scans together
    """
    n1 = len(scan1)
    n2 = len(scan2)
    # "+" joins lists:
    nodes = scan1.nodes + scan2.nodes
    edges = scan1.edges + [ (i + n1, j + n1) for i,j in nodes ]
    return Scan( nodes, edges )
    

class Project( object ):
    """
    Holds all the information about the project

       scans / experiment data
       processing steps and results

    Mostly links to external files ?
    """
    def __init__(self):
        self.Scans = []
        self.processing = []
        
    def load(self):
        pass
    
    def save(self):
        pass

