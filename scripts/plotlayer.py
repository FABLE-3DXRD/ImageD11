#!/usr/bin/env python

from __future__ import print_function



from matplotlib.pylab import *
from matplotlib.patches import Ellipse
from ImageD11.grain import *
allaxes = [ ((1,0,0), (0,1,0)) ,
            ((0,0,1), (1,0,0)) ,
            ((0,1,0), (0,0,1)) ]

def make_ellipses( grains , scalet = 1.0, scaler=1.0 ):
    el = []
    vol = []
    
    for g in grains:
        try:
            vol.append( float(g.intensity_info.split()[14]) )
        except:
            vol.append( int(g.npks) )
    # Sort by volume to plot biggest first (at back)
    tmp = list(zip( vol, grains))
    tmp.sort()
    from ImageD11.indexing import ubitoU
    for vol, g in tmp[::-1]:
        size = pow(vol, 1.0/3.0)
        el.append( (Ellipse( xy = g.translation[0:2],
                            width = size*scalet,
                            height = size*scalet,
                             ),
                    Ellipse( xy = (g.ubi[0,1],g.ubi[0,2]),
                             width = size*scaler,
                             height = size*scaler,
                             ), g)
                   )
        
        # print g.translation[0:2]
    return el

def graincolor( g ):
    from ImageD11.indexing import calc_drlv2, ubitoB
    hkl = array( [ ( 1,0,0 ),
                   ( 0,1,0 ),
                   ( 0,0,1 ),
                   ( 1,1,0 ),
                   ( 0,1,1 ),
                   ( 1,0,1 ),
                    ] )
    from numpy.linalg import inv
    
    B = inv(ubitoB( g.ubi ))

    gcalc = dot( B.T, hkl.T )
    d = calc_drlv2( g.ubi, gcalc.T )
    red = sum(d*d)

    v = 1/sqrt(2)
    
    u = array( [[   v,   0,  v ],
                [   0,   1,  0 ],
                [  -v,   0,  v ]] )
    
    gcalc = dot( dot(u,B).T , hkl.T )
    d = calc_drlv2( g.ubi, gcalc.T )
    green = sum(d*d)
    
    u2 = array([ [ 1.0, 0.0, 0 ],
                 [ 0. , v,  -v ],
                 [ 0.0, v,   v ] ])
    
    gcalc =  dot( dot(u2 ,B).T , hkl.T )
    d = calc_drlv2( g.ubi, gcalc.T )
    blue = sum(d*d)
    #print "Grain",red,green,blue,
    # Normalise the length of the vector to be 1
    modc = sqrt(red*red+blue*blue+green*green)
    #print "Grain",red,green,blue
    return red/modc, green/modc, blue/modc
    


    
def plubi(uf, first = False , fig=None):
    
    gl = read_grain_file(uf )
    pl = []
    for g in gl:
        if g.translation[2] > -15:
           pl.append(g)
           
    xyz = array(   [ g.translation for g in pl ] )
    ubis =  array( [ g.ubi for g in pl ] )
    # print ubis.shape
    el = make_ellipses( pl , scalet = 2.0, scaler=0.005 )

    ax = fig.add_subplot( 111, aspect = 'equal' )
    
    for et, er, g in el:
        ax.add_artist(et)
        et.set_clip_box(ax.bbox)
        gc = graincolor( g )
        et.set_facecolor( gc )
        et.set_alpha( 0.75 )

    min_x  = xyz[:,0].min()
    min_y  = xyz[:,1].min()
    max_x  = xyz[:,0].max()
    max_y  = xyz[:,1].max()
    dx = max( max_x - min_x, 1)
    dy = max( max_y - min_y, 1)
    ax.set_xlim( min_x - dx/10. , max_x + dx/10.)
    ax.set_ylim( min_y - dy/10. , max_y + dy/10.)

    title(uf)



    
    
if __name__=="__main__":

    import sys, os, glob, time

    

    
    
    fig = figure()
    
    plubi(  sys.argv[1], True, fig)
    if len(sys.argv[2])>2:
        savefig( sys.argv[2] )

