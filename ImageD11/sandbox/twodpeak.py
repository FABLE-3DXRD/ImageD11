
import numpy as np

# TODO : derivatives, longer tails


def twodpeak( pars, x, y ):
    """
    pars : (m0, m1x, m1y, m2xx, m2xy, m2yy)
    x : x-coordinates, probably a 2D array (slow or tth)
    y : y-coordinates, probably a 2D array (fast or eta)

    Gaussian would be:
    dx = x - m1x
    dy = y - m1y
    peak  = m0 * exp( -0.5 (dx*dx/m2xx + dx*dy/m2xy + dy*dy/m2yy ) ) 
    """
    assert x.shape == y.shape
    (m0, m1x, m1y, m2xx, m2xy, m2yy) = pars
    dx = x - m1x
    dy = y - m1y
    # [dx, dy].T . [[ m2xx, m2xy],[m2xy, m2yy]] . [dx, dy]
    E = np.linalg.inv( (( m2xx, m2xy),(m2xy, m2yy) ) )
    det_E = m2xx*m2yy - m2xy*m2xy
    peak = np.exp(-0.5*( E[0,0]*dx*dx + E[0,1]*2*dy*dx + E[1,1]*dy*dy ) )
    norm = np.sqrt( 4 * np.pi * np.pi * det_E )
    return m0 * peak / norm

def testpeak( m0   = 20., 
              m1x  = 40.,
              m1y  = 71.5,
              m2xx = 8.,
              m2xy = 1.,
              m2yy = 3. ):
    x, y = np.mgrid[0:int(m1x*2), 0:int(m1y*2)]
    pars = ( m0, m1x, m1y, m2xx, m2xy, m2yy )
    return twodpeak( pars, x , y )

def check_normalised():
    m0=100
    pk = testpeak(m0=m0)
    pks = pk.sum()
    if not np.isclose( pks, m0 ):
        print(pks,m0)
        import pylab as pl
        pl.imshow( pk )
        pl.show()
    else:
        print("Normalisation looks OK")

    
# Compare to ImageD11 peaksearch moments / image moments
def check_moments():
    from ImageD11.cImageD11 import connectedpixels, blobproperties, \
        blob_moments, s_1, s_I, s_fI, s_sI, s_ssI, s_sfI, s_ffI,    \
        s_raw, f_raw, m_ss, m_ff, m_sf
    from ImageD11 import cImageD11
    p = {'m0'   : 10.  ,
         'm1x'  : 80.  ,
         'm1y'  : 77.5 ,
         'm2xx' : 5.,
         'm2xy' : 2.,
         'm2yy' : 3., }
    pk = testpeak(**p).astype(np.float32)
    blob = np.zeros( pk.shape, np.int32 )
    npk = connectedpixels( pk, blob, 0., 0 )
    res = blobproperties( pk, blob, npk, 0 )
    err = blob_moments( res )
    mss = res[0,s_ssI]/res[0,s_I] - res[0,s_raw]*res[0,s_raw]
    mff = res[0,s_ffI]/res[0,s_I] - res[0,f_raw]*res[0,f_raw]
    msf = res[0,s_sfI]/res[0,s_I] - res[0,f_raw]*res[0,s_raw]
    assert np.isclose( res[0,s_raw], p['m1x'] )
    assert np.isclose( res[0,f_raw], p['m1y'] )
    assert np.isclose( mss, p['m2xx'] )
    assert np.isclose( mff, p['m2yy'] )
    assert np.isclose( msf, p['m2xy'] )
    print("Moments look OK")
    if 0:
        for adr, name in  zip(
                (s_1,s_I,s_fI,s_sI,s_ssI,s_sfI,s_ffI,s_raw,f_raw,m_ss,m_ff,m_sf),
                "s_1,s_I,s_fI,s_sI,s_ssI,s_sfI,s_ffI,s_raw,f_raw,m_ss,m_ff,m_sf".split(",")):
            print(name, res[ 0, adr ])
        print( res[0,s_raw], p['m1x'] )
        print( res[0,f_raw], p['m1y'] )
        print( [[mss, msf],[msf,mff]] )
        import pylab as pl
        pl.imshow( np.log(pk), origin='lower' )
        pl.show()

        
def test():

    check_moments()
    check_normalised()
    
if __name__=="__main__":
    test()
        
    


