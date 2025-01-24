
from __future__ import print_function


# ImageD11 Software for beamline ID11
# Copyright (C) 2008  Jon Wright
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  0211-1307  USA

"""
A script which goes through a peak search output to recover the bounding
boxes found and write them in separate edf files
Puts in some additional header keys
"""


from ImageD11.columnfile import columnfile
from fabio.openimage import openimage
from fabio import edfimage
import numpy as np
import sys, bisect

# First arg is a filtered peaks file containing the spots for this grain
flt = columnfile(sys.argv[1])
sptfile = open(sys.argv[2],"r")
if len(sys.argv)>3:
    bkg = openimage(sys.argv[3])
    bkg.data = bkg.data.astype(np.float32)
else:
    bkg = None

for col in ["sc","fc","omega",  # Slow, fast and omega center positions
            "Min_o","Max_o",    # Bounding box
            "Min_s","Max_s",
            "Min_f","Max_f" ] : 
    # check we have the columns we need
    assert hasattr(flt, col),col

print("Looking for",flt.nrows,"spots")

# Now get the list of omega angles for the images
omegas = {} ; filenames = {}; spots = {}
for line in sptfile:
    if line.find("# File")==0:
        filename = line.split(" ")[2].rstrip()
    if line.find("# Omega ")==0:
        om = float(line.split("=")[1])
        # print om, filename
        filenames[om] = filename
        omegas[filename] = om    
        spots[om] = []

oms = list(filenames.keys())
oms.sort()
print("Found",len(oms),"distinct omega values")


#
#
#
# Quantix
PIXELSIZE = 4.3 

# add a 5 pixel border on bounding box
BORDER = 5

# Dict of header keys
# If type == "String": just copy directly
# If type == [list] : columnfile labels l[0] put through function l[1]

mappings = { "HeaderID":"EH:000001:000000:000000",
            "Image":"1",
            "title": [ ("spot3d_id",), lambda x: "Spot ID %d"%(x) ],
            # Dim_1 = 85 ;
            # Dim_2 = 29 ;
            # DataType = FloatValue ;
            # ByteOrder = LowByteFirst ;
            # Size = 9860 ;
            "SpotID": [("spot3d_id",), lambda x: "%d"%(x) ],
            "COM_y" : [("fc",) ,       lambda x: PIXELSIZE*x ],
            "COM_z" : [("sc",) ,       lambda x: PIXELSIZE*x ], 
            "COM_omega" : [("omega",) ,       lambda x: x ],
            "spot_integint" : [("sum_intensity",),    lambda x: x ],
            "spot_pixelnum" : [("Number_of_pixels",), lambda x: x ], # Not bounding box pixels
            "ave_intensity" : [("avg_intensity",),    lambda x: x ],
#          col_colbegin = 408 ;
#            col_colend = 493 ;
#          row_rowbegin = 341 ;
#            row_rowend = 370 ;
#         y_colrowbegin = 3.692000E-04 ; # 284 pixel
#           y_colrowend = 2.587000E-04 ; # 199 pixel
#         z_colrowbegin = 4.707300E-04 ; # 362.1 pixel
#           z_colrowend = 5.084300E-04 ; # 391.1 pixel
#          DetSOdist_mm = 4.02300E+00 ;
#       OrigDetHoriz_px = 6.920000E+02 ;
#        OrigDetVert_px = -2.110000E+01 ;
"omega_min_deg" : [("Min_o",) ,  lambda x: x ],
"omega_max_deg" : [("Max_o",) ,  lambda x: x ],
            "IsTouchingAnotherSpot" :"0",
            "IsAtBound_y":"1",
            "IsAtBound_z":"1",
            # 0+0=0, 1+0=1, 1+1=2->1
            "IsAtBound_omega" : [ ("onfirst","onlast",), lambda a,b :  1-int(a or b) ],
#            "Convexity" :"9.850187E-01",
#            "Convexity_trials" :"267",
#            "Convexity_cutlevel" :"1.147189E+00",
            "IsSmoothed":"0",
            "PixSizeHoriz_um": "%f"%(PIXELSIZE),
            "PixSizeVert_um": "%f"%(PIXELSIZE),
            }


def writespot(filename, roi, box, hdr, flt, ispt):
    """
    Filename  - target filename
    roi       - the roi in the original image
    box       - the pixels
    hdr       - original file header
    flt       - peaksearch output in columnfile object
    spot3d_id - which peak this is
    """
    # Make the header
    print("writing",filename)
    myheader = {"ByteOrder":hdr["ByteOrder"]}
    ks = []
    for k in edfimage.MINIMUM_KEYS:
        ks.append(k)
        if k not in mappings:
            try:
                myheader[k] = edfimage.DEFAULT_VALUES[k]
            except:
                pass
    for k in flt.titles:
        myheader[k] = getattr(flt, k)[ispt]
        ks.append(k)
    ks.append("OriginalFile")
    myheader["OriginalFile"] = hdr['filename']
    myheader["DataType"] = "FloatValue"
    for k, v in  [ ("ROI_slow_start", roi[0].start ),
                   ("ROI_slow_stop" , roi[0].stop  ),
                   ("ROI_fast_start", roi[1].start ),
                   ("ROI_fast_stop" , roi[1].stop  ) ]:
        ks.append(k)
        myheader[k] = int(v)
    km = list(mappings.keys())
    km.sort()
    for k in km:
        if k not in ks:
            ks.append(k)    
        if k in km:
            v = mappings[k]
            if type(v) == type("String"):
                myheader[k] = v
                continue
            if type(v) == type(["List"]): # Funky
                try:
                    args = [getattr(flt, a)[ispt] for a in list(v[0])]
                    myheader[k] = v[1](*args)
                except:
                    print(v[0])
                    raise
    #b = box.astype(np.float32)
    #print np.maximum.reduce(b), np.minimum.reduce(b)
    o = edfimage.edfimage( box.astype(np.float32) , myheader )
    o.header_keys = ks
    o.header_keys.append("filename")
    o.write(filename, force_type=None)

    
    


# Now we want to convert the loop over spots into a loop over images
# Find out which spots are in each image

for i in range(flt.nrows): # looping over spots
    first =  bisect.bisect_left( oms, flt.Min_o[i] ) 
    last  = bisect.bisect_right( oms, flt.Max_o[i] ) 
    #print flt.Min_o[i], oms[first_om], flt.Max_o[i], oms[last_om]
    for o in oms[first:last]:
        spots[o].append(i)

boxes = {}
border = 5

integ_int = {}

for o in oms:
    print(o, filenames[o], spots[o])
    data_obj = openimage(filenames[o])
    if bkg is not None:
        data_obj.data = data_obj.data.astype(np.float32) - bkg.data.astype(np.float32)
    print(np.maximum.reduce(np.ravel(data_obj.data)), end=' ')
    print(np.minimum.reduce(np.ravel(data_obj.data)))
    for i in spots[o]:
        roi = ( slice( max(0, flt.Min_s[i] - border) ,
                       min( flt.Max_s[i] + border , data_obj.dim2 ) ) ,
                slice( max(0, flt.Min_f[i] - border) , 
                       min( flt.Max_f[i] + border ,  data_obj.dim1 ) ) )

        if o == flt.Min_o[i]:
            boxes[i] = data_obj.data[roi]
        else:
            boxes[i] += data_obj.data[roi]

        if o == flt.Max_o[i]:
            # This is the last image it is on
            writespot( "%d.edf"%(i), roi, boxes[i], data_obj.header,
                    flt, i )
            
            integ_int[i] = np.sum(boxes[i])

            del boxes[i]

        if False:
            from matplotlib.pylab import imshow, show
            imshow(boxes[i])
            show()
            sys.exit()

integ_int = [ integ_int[i] for i in range(flt.nrows) ]
flt.addcolumn( integ_int, "box_intensity")
flt.writefile(sys.argv[1]+".new")
