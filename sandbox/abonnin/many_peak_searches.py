
from __future__ import print_function


for i in range(180):
    f=147*i
    l=147*i+146
    HOME = '/users/abonnin/lien_ma1063/diffraction/IRIS4_rect5_1/'
    stem = HOME + 'data_XRD-CT/IRIS4_rect5_1_'
    dark = HOME + 'dark_0000.edf'
    spline = HOME + 'frelonID22_dist06102010.spline'
    outfile = HOME + "omegapeaksearch/%d.spt"%(i)
    cmd = 'peaksearch.py -n %s -f %d -l %d -d %s -s %s -T 0 -S 0.2 -t 500 -o %s'%(
        stem, f, l, dark, spline, outfile )

    print(cmd)
    
