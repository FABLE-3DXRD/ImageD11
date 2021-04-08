
from __future__ import print_function
import yaml, pprint


my_yml = """

Detectors:
    Frelon21 : 
        name : "Frelon21"
        type : "2Ddetector"
        dark : "/data/id11/.../dark1s.edf"
        flat : "/data/id11/.../flat.edf"
        spline : "/data/id11/.../frelon21_mar16.spline"
    Frelon4M :
        name : "Frelon4M"
        type  : "2Ddetector"
        dark : "/data/id11/.../dark1s.edf"
        flat : "/data/id11/.../flat.edf"
        spline : "/data/id11/.../Frelon4M.spline"
    pico0 :
        name : "pico0"
        type : "0Ddetector"
        # e.g. photons to pA current converter
        gain : 1.e7
    mca :
        name : "mca"
        type : "1Ddetector"
        nchannels : 4096
        inifile : SDD_53keV.ini
        gain :  #  4096/53.0

Positioners:
    FF_Detector_Mount :
        - {type : translation, name : ffdtx1, axis : [1.0, 0.0, 0.0], pos : 20.0}
        - {name : ffdtz1,  type: translation, axis : [0.0, 0.0, 1.0]}
        - {name : ffdtilt, type: rotation, axis : [0.0, 1.0, 0.0]} 
    3DXRD_Huber_Tower :
        - {name : diffty, type: translation, axis : [0.0, 1.0, 0.0] }
        - {name : difftz, type: translation, axis : [0.0, 0.0, 1.0] }
        - {name : diffry, type: rotation, axis : [0.0, 1.0, 0.0] }
        - {name : diffrz, type: rotation, axis : [0.0, 0.0, 1.0] }
        - {name : samry,  type: rotation, axis : [0.0, 1.0, 0.0] }
        - {name : samrx,  type: rotation, axis : [1.0, 0.0, 0.0] }
        - {name : samtx,  type: translation, axis : [1.0, 0.0, 0.0] }
        - {name : samty,  type: translation, axis : [0.0, 1.0, 0.0] }
        - {name : samtz,  type: translation, axis : [0.0, 0.0, 1.0]  }
    Fable_diffractometer :
        - { name : wedge, type : rotation, axis : [0.0, 1.0, 0.0] }
        - { name : omega, type : rotation, axis : [0.0, 0.0, 1.0] }
        - { name : t_x, type : translation, axis : [1.0, 0.0, 0.0] }
        - { name : t_y, type : translation, axis : [0.0, 1.0, 0.0] }
        - { name : t_z, type : translation, axis : [0.0, 0.0, 1.0] }
    Fable_detector :
        - { name : distance, type : translation, axis : [1.0, 0.0, 0.0] }
        - { name : tilt_x, type : rotation, axis : [1.0, 0.0, 0.0] }
        - { name : tilt_y, type : rotation, axis : [0.0, 1.0, 0.0] }
        - { name : tilt_z, type : rotation, axis : [0.0, 0.0, 1.0] }
        - { name : Oij, type : positioner,
            mat4 :  [[1,   0,   0, 0],
                     [0, o22, o21, 0],
                     [0, o12, o11, 0],
                     [0,   0,   0, 1]] }
        - { name : z_size, type : scale, axis : [0, 0, 1] }
        - { name : y_size, type : scale, axis : [0, 1, 0] }
        - { name : z_center, type : translation, axis : [0.0, 0.0,-1.0] }
        - { name : y_center, type : translation, axis : [0.0,-1.0, 0.0] }

Fe3O4: 
    - name : magnetite
    - space_group : Fd-3m
    - composition : { 'Fe':3, 'O':4 }
    - unit_cell : [ 8.39, 8.39, 8.39, 90., 90., 90. ]

Silicon:
    - name : silicon
    - space_group : Fd-3m
    - composition : { 'Si':1 }
    - unit_cell : [ 5.43094,5.43094, 5.43094,90,90,90 ] 

Experiment:
  Radiation :
    # energy : 0.124                  # Units match unit cell 
    wavelength : 0.124
    direction : [ 1.0, 0.0, 0.0 ] 
    bandpass : 2.0e-3
    divergence : 1.0e-6             # Units?
  SampleMount :                # from 0,0,0 to sample position
    3DXRD_Huber_Tower
  Detectors :                  # from 0,0,0 to detector position
    - { positioner : FF_Detector_Mount, detector : Frelon21 }
  Scans :
    scan1 :
      - Motor  : diffrz
        Step  : 0.25
        Start : 0.0
      - measurement :
        - monitor :
            name : pico6
            normvalue : 1e7  # to normalise to
            dataURL : "place:///to_find_the_data.dat"
        - Frelon21 : 
            dark : dark1s.edf
            darkoffset : 12
            binning : [1, 1]
            flips : [No, No]
            dataURL :
              namefmt : "{stem:s}{pass:d}_{number:04d}.edf.gz"
              folder : /data/id11/inhouse3/blah/toto
              stem : toto17_
              first : 0
              last  : 899
              interlaced : Yes
              iflip : No

    scan2 :
      - Motor  : diffry
        Step  : 0.05
        Start : 0.0
      - measurement :
        - monitor :
            name : pico6
            normvalue : 1e7  # to normalise to
            dataURL : "place:///to_find_the_data.dat"
        - Frelon21 : 
            dark : dark1s.edf
            darkoffset : 12
            binning : [1, 1]
            flips : [No, No]
            dataURL :
              namefmt : "{stem:s}{number:04d}.edf"
              folder : /data/id11/inhouse3/blah/toto2
              stem : toto17_
              first : 0
              last  : 20


    scan3 :
      - Motorpos : { "diffty": 0.1 }
      - Motor  : diffrz
        Step  : 0.05
        Start : 0.0
      - measurement :
        - monitor :
            name : pico6
            normvalue : 1e7  # to normalise to
            dataURL : "place:///to_find_the_data.dat"
        - Frelon21 : 
            dark : dark1s.edf
            darkoffset : 12
            binning : [1, 1]
            flips : [No, No]
            dataURL :
              namefmt : "{stem:s}{number:04d}.edf"
              folder : /data/id11/inhouse3/blah/toto3
              stem : toto17_
              first : 0
              last  : 20

"""

ImageD11_project_schema = """
 TODO : figure out how to format and validate this stuff properly
"""



y = yaml.load(my_yml)
pprint.pprint(y)




class ImageD11Project( object ):

    def __init__(self, project_dict = None, filename=None ):
        """ Takes filename or existing project dictionary """
        self.project_dict = { }
        if filename is not None:
            self.load( filename )
        if project_dict is not None:
            self.project_dict.update( project_dict )

    def load(self, filename):
        """ Updates the existing project from the file (e.g. adds and overwrites) """
        with open(filename, "r") as stream:
            pd = yaml.load( stream )
        self.project_dict.update( pd )
        
    def save(self, filename):
        """ Writes the current project to a file """
        with open(filename, "w") as stream:
            stream.write( yaml.dump( self.project_dict, default_flow_style=False, encoding='utf-8' ) )

    def validate(self):
        """ Check that we have what we need in the dictionary 
        TODO - try to use schema from nexus or cif or something ? """
        pass




