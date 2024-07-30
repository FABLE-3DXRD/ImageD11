
from __future__ import print_function
import yaml, pprint


my_yml = """

# Detectors
Frelon21 : 
    name : "Frelon21"
    type : "2Ddetector"
    dark : "/data/id11/.../dark1s.edf"
    flat : "/data/id11/.../flat.edf"
    spline : "/data/id11/.../frelon21_mar16.spline"

Frelon4M :
    name : "Frelon4M"
    type : "2Ddetector"
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

# Experiment Stations
3dxrd: 
    - name : 3dxrd
    - wavelength : 0.1234
    - Positioners:
        - omega : 0.0
        - wedge : 0.0
        - samtz : 0.0
    - Detectors:
        - name : Frelon21 
          detx : 120.0  
        - name : Frelon4M
          ffdtz1 : 220.0

Fe3O4: 
    - name : magnetite
    - space_group : Fd-3m
    - composition : { 'Fe':3, 'O':4 }
    - unit_cell : [ 8.39, 8.39, 8.39, 90., 90., 90. ]

Silicon:
    - name : silicon
    - space_group : Fd-3m
    - composition : { 'Si':1 }
    - unit_cell : [ 5.43094,5.430945.43094,90,90,90 ] 

scan1:
    - instrument:
      - 3dxrd
    - sample:
      - Silicon
      - Fe3O4
    - measurement:
      - Frelon21: ["test0001.edf", "test0002.edf", "test0003.edf" ]
      - omega:    [ 0.0 ,          1.0,            2.0]
      - epoch:    [ 0.0 ,          5.0,            10.0]
      - srcur:    [ 190.0 ,      189.9,           189.8]

"""

y = yaml.safe_load(my_yml)
pprint.pprint(y)
