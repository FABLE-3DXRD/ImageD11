#!/usr/bin/env python

from __future__ import print_function

help = """
Reads in a series of flt files to merge together different threshold levels

Takes all peaks from the highest level.
For lower levels, it adds new peaks which are found if they are not 
already present
If they correspond to a previous peak, it chooses between the old and the
new according to the centre of mass. If the lower thresholded peak's centre
of mass does not agree with the higher within pixel_tol, then the lower is
rejected

Peaks are identified as unique via their co-ordinates of their strongest
pixel (f,s,omega)
"""


from ImageD11 import  transformer
from ImageD11.columnfile import newcolumnfile, columnfile
import numpy
import sys, time

try:
    pars = sys.argv[1]
    stem = sys.argv[2]
    outf = sys.argv[3]
    dpix = float(sys.argv[4])
    thres = [int(v) for v in sys.argv[5:]]
except:
    print("Usage: pars stem outputfile pixel_tol thresholds1  thresholds2 ...")
    print(help)
    sys.exit()

assert outf[-4:] == ".flt", """output file should end in ".flt" """

thres.sort()
thres = thres[::-1]

print("Using parameters",pars)
print("Merging files", end=' ')
for v in thres:
    print("%s_t%d.flt"%(stem,v), end=' ')
print()
print("Into output file %s"%(outf))

#if raw_input("OK? [y/n]") not in ["Y","y"]:
#    sys.exit()

allpks = open(outf,"w")

allpeaks = {}
always_ignore = {}

goodthres = []

for v in thres:

    mytransformer = transformer.transformer()
    mytransformer.loadfileparameters( pars )
    
    flt = "%s_t%d.flt"%(stem,v)
    print(flt, end=' ')
    try:
        tc = columnfile( flt )
        if tc.nrows == 0:
            print("Skipped",tc," no peaks")
            continue
    except:
        print("Skipped",v," Exception reading",flt)
        continue
    goodthres.append( v )
    mytransformer.loadfiltered( flt )
    mytransformer.compute_tth_eta( )
    mytransformer.addcellpeaks( )

    print("npeaks", mytransformer.colfile.nrows, end=' ')
        
    # mytransformer.write_colfile(flt2)

    f = mytransformer.colfile.titles.index('sc')
    s = mytransformer.colfile.titles.index('fc')
    titles = mytransformer.colfile.titles
    nignore = 0
    nnew = 0
    nold = 0
    for i in range(mytransformer.colfile.nrows):
        # Position of max intensity
        key = ( int(mytransformer.colfile.IMax_o[i]*100) ,
                int(mytransformer.colfile.IMax_s[i]) ,
                int(mytransformer.colfile.IMax_f[i]) )
        
        if key in always_ignore:
            nignore = nignore + 1
            continue
        
        if key in allpeaks:
            if v is goodthres[0]:
                print(key)
                print("duplicate")
                #raise
            # This peak is already found
            # Should we replace it, or trash the lower threshold ??
            # previous is allpeaks[key]
            # current is mytransformer.colfile.bigarray[:,i]
            df = allpeaks[key][f]-mytransformer.colfile.bigarray[f,i]
            ds = allpeaks[key][s]-mytransformer.colfile.bigarray[s,i]
            dist2 = df*df + ds*ds
            if dist2 > dpix*dpix:
                # ignore the weaker peak
                # print "Ignoring weaker",
                nignore = nignore + 1
                always_ignore[key] = 1
            else:
                # Replace the stronger peak with the weaker peak
                allpeaks[key] = mytransformer.colfile.bigarray[:,i].copy()
                nold = nold + 1
        else:
            nnew = nnew + 1
            allpeaks[key] = mytransformer.colfile.bigarray[:,i].copy()
            
    print("total peaks",len(list(allpeaks.keys())), "ignored", nignore, "new",nnew, "replacements",nold)
    assert nignore + nold + nnew == mytransformer.colfile.nrows

                                                        

keys = list(allpeaks.keys())

keys.sort()

c = newcolumnfile( titles )

assert len(titles) == len(allpeaks[keys[0]])

bigarray = [allpeaks[k] for k in keys]

c.bigarray = numpy.array(bigarray).T.copy()
print(c.bigarray.shape)
c.nrows = len(keys)
c.set_attributes()
c.setcolumn(numpy.array(list(range(len(keys)))),"spot3d_id")
c.writefile( outf )

mytransformer = transformer.transformer()
mytransformer.loadfileparameters( pars )

mytransformer.loadfiltered( outf )
mytransformer.compute_tth_eta( )
mytransformer.addcellpeaks( )
mytransformer.computegv( )

mytransformer.write_colfile( outf )
mytransformer.savegv( outf.replace(".flt",".gve" ) )
