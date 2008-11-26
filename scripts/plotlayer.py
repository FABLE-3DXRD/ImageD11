


from matplotlib.pylab import *
from matplotlib.patches import Ellipse
from ImageD11.grain import *
    
def make_ellipses( grains , scalet = 1.0, scaler=1.0 ):
    el = []
    vol = []
    
    for g in grains:
        try:
            vol.append( float(g.intensity_info.split()[14]) )
        except:
            vol.append( int(g.npks) )
    # Sort by volume to plot biggest first (at back)
    tmp = zip( vol, grains)
    tmp.sort()
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

    gcalc = dot( B, hkl.T )
    d = calc_drlv2( g.ubi, gcalc.T )
#    print sum(d*d)
    red = 1 - sum(d*d)

    v = 1/sqrt(2)
    
    u = array( [[ v, v, 0 ],
                [-v, v, 0 ],
                [ 0.0, 0.0, 1 ]] )
    
    gcalc = dot( dot(u,B) , hkl.T )
    d = calc_drlv2( g.ubi, gcalc.T )
#    print sum(d*d)
    green = 1 - sum(d*d)
    
    u2 = array([ [ 1.0, 0.0, 0 ],
                 [ 0. , v,-v ],
                 [ 0.0, v, v ] ])
    
    gcalc =  dot( dot(u2 ,B) , hkl.T )
    d = calc_drlv2( g.ubi, gcalc.T )
 #   print sum(d*d)
    blue = 1 - sum(d*d)

    return clip( (red, green, blue), 0, 1 )
    


    
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

        
    
    ax.set_xlim( -400, 400)
    ax.set_ylim( -400, 400)

    title(uf)



    
    
if __name__=="__main__":

    import sys, os, glob, time

    

    
    
    fig = figure()
    
    plubi(  sys.argv[1], True, fig)
    show()
