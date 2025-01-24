
from __future__ import print_function

from ImageD11.columnfile import columnfile
from fabio.openimage import openimage
import sys

# peak number 16

from matplotlib.pylab import figure,show,imshow,show, plot
DATA="../data_XRD-CT/IRIS4_rect5_1_"

def integrate_peak(number, colfile, bbs=10):
    i = list(colfile.spot3d_id.astype(int)).index( number )
    om = colfile.omega[i]
    print("omega",i, om)
    bb = [colfile.Min_s[i],
          colfile.Max_s[i],
          colfile.Min_f[i],
          colfile.Max_f[i],
          colfile.Min_o[i],
          colfile.Max_o[i]]

    first_image = int(om) * 147
    last_image = first_image + 146
    print(bb)

    s = colfile.s_raw[i]
    f = colfile.f_raw[i]
    
    intensity = []
    for j in range(first_image, last_image+1):
        print(j, end=' ')
        sys.stdout.flush()
        im = openimage( DATA + "%04d.edf"%( j ) ).data
        intensity.append(im[ s-bbs:s+bbs , f-bbs:f+bbs ].mean())

        if 0 and j == 54:

            imshow(im,vmax=1000)
            figure()
            imshow(im[ s-20:s+20 , f-20:f+20],vmax=1000)
            show()

            plot(intensity)
            show()
    return om, intensity

c  = columnfile("jon.flt")
number = 16


allpeaks = []
oms = []
for line in open("grain1.pks").readlines():
    try:
        num = int(line)
        om, intensity = integrate_peak( num, c )
        allpeaks.append( intensity )
        oms.append( om )
    except:
        pass


import numpy

ar = numpy.array( intensity )
numpy.save( "grain1.sino", ar)
om = numpy.array( oms )
numpy.save("grain1.omega",om )
