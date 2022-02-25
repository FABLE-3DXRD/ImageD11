from __future__ import print_function

import numpy as n
def compute_grain_origins(omega, wedge = 0.0, chi = 0.0,
                          t_x = 0.0, t_y = 0.0, t_z = 0.0):
    """
    # print "Using translations t_x %f t_y %f t_z %f"%(t_x,t_y,t_z)
    # Compute positions of grains
    # expecting tx, ty, tz for each diffraction spot
    #
    # g =  R . W . k
    #  g - is g-vector w.r.t crystal
    #  k is scattering vector in lab
    #  so we want displacement in lab from displacement in sample
    #  shift =  W-1  R-1 crystal_translation
    #
    # R = ( cos(omega) , sin(omega), 0 )
    #     (-sin(omega) , cos(omega), 0 )
    #     (         0  ,         0 , 1 )
    #
    # W = ( cos(wedge) ,  0  ,  sin(wedge) )
    #     (         0  ,  1  ,          0  )
    #     (-sin(wedge) ,  0  ,  cos(wedge) )
    #
    # C = (         1  ,          0  ,       0     ) ??? Use eta0 instead
    #     (         0  ,   cos(chi)  ,  sin(chi)   )  ??? Use eta0 instead
    #     (         0  ,  -sin(chi)  ,  cos(chi)   )  ??? Use eta0 instead
    """
    w=n.radians(wedge)
    WI = n.array( [ [ n.cos(w),         0, -n.sin(w)],
                    [      0,           1,         0],
                    [ n.sin(w),         0,  n.cos(w)] ] , n.float)
    c=n.radians(chi)
    print ("WI",WI)
    CI = n.array( [ [      1,            0,         0],
                    [      0,     n.cos(c), -n.sin(c)],
                    [      0,     n.sin(c),  n.cos(c)] ] , n.float)
    print( "CI",CI)
    t   = n.zeros((3,omega.shape[0]),n.float) # crystal translations
    # Rotations in reverse order compared to making g-vector
    # also reverse directions. this is trans at all zero to
    # current setting. gv is scattering vector to all zero
    om_r = n.radians(omega)
    print( "om_r",om_r)
    # This is the real rotation (right handed, g back to k)
    t[0,:] = n.cos(om_r)*t_x - n.sin(om_r)*t_y
    t[1,:] = n.sin(om_r)*t_x + n.cos(om_r)*t_y
    t[2,:] =                                  t_z
    print ("om.t",t)
    if chi != 0.0:
        c = n.cos(n.radians(chi))
        s = n.sin(n.radians(chi))
        u = n.zeros(t.shape,n.float)
        u[0,:]= t[0,:]  
        u[1,:]=        c * t[1,:]    + -s * t[2,:]
        u[2,:]=        s * t[1,:]    +  c * t[2,:]
        t = u
    print ("chi.om.t",t)
    if wedge != 0.0:
        c = n.cos(n.radians(wedge))
        s = n.sin(n.radians(wedge))
        u = n.zeros(t.shape,n.float)
        u[0,:]= c * t[0,:]           + -s * t[2,:]
        u[1,:]=            t[1,:]
        u[2,:]= s * t[0,:]           +  c * t[2,:]
        t = u
    print( "wedge.chi.om.t",t)
    return t

import ctypes, numpy as np

dll = ctypes.CDLL("./diffraction.so")

REAL = ctypes.c_double
NPREAL = np.float64

dll.c_choosegrains.argtypes = [
#void choosegrains( vector XL[], real omega[], 
   ctypes.POINTER(REAL),
   ctypes.POINTER(REAL),
#        real wedge, real chi, real wvln,
   REAL, REAL, REAL,
#        rmatrix UB[], rmatrix UBI[], vector T[], 
   ctypes.POINTER(REAL),
   ctypes.POINTER(REAL),
   ctypes.POINTER(REAL),
#        int ngrains, int npks,
   ctypes.c_int, ctypes.c_int,
#        int labels[], vector hkls[]){ 
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(REAL)
        ]

def flatfloat(a, shape=None):
    if shape is None:
        s = a.shape
    else:
        s=shape
    flat = np.array( np.array(a).ravel(), NPREAL )
    flat.shape = s
    return flat 


if __name__=="__main__":
    import sys, numpy as np

    from ImageD11.grain import read_grain_file
    from ImageD11.columnfile import columnfile
    from ImageD11.parameters import parameters
    from ImageD11.transform import compute_xyz_lab

    colfile = sys.argv[1]
    grainfile = sys.argv[2]
    parfile = sys.argv[3]

    pars = parameters()
    pars.loadparameters(parfile)
    
    wedge = np.radians(float(pars.get("wedge")))
    chi   = np.radians(float(pars.get("chi")))
    wvln  = float(pars.get("wavelength"))
    
    c = columnfile( colfile )
    try:
        sf = c.sc,c.fc
    except:
        sf = c.xc,c.yc
#    c.nrows = 2
#    c.bigarray = c.bigarray[:,:2]
    print (c.bigarray.shape)
    c.set_attributes()

    XL = flatfloat( compute_xyz_lab( sf , **pars.parameters ).T )
    om = flatfloat( np.radians(c.omega ))


#    XL[0,:]=[41,42,43]

    grains = read_grain_file( grainfile )
    UBIs = flatfloat( [ g.ubi for g in grains ], shape=(len(grains),3,3) )
    UBs =  flatfloat( [ np.linalg.inv(g.ubi) for g in grains ], shape=(len(grains),3,3) )
    Ts  =  flatfloat( [ g.translation for g in grains ], shape=(len(grains),3))


    labels = np.zeros( c.nrows, np.int32 )
    hkls = np.zeros( ( c.nrows, 3 ), NPREAL )

    print( XL.shape)
    print( om.shape)
    print( UBs.shape)
    print( UBIs.shape)
    print( Ts.shape)
    print( labels.shape)
    print( hkls.shape)

    dll.c_choosegrains( 
            XL.ctypes.data_as( ctypes.POINTER( REAL) ),
            om.ctypes.data_as( ctypes.POINTER( REAL) ),
            wedge, chi, wvln,
            UBs.ctypes.data_as( ctypes.POINTER( REAL) ),
            UBIs.ctypes.data_as( ctypes.POINTER( REAL) ),
            Ts.ctypes.data_as( ctypes.POINTER( REAL) ),
            len(grains), c.nrows,
            labels.ctypes.data_as( ctypes.POINTER( ctypes.c_int) ),
            hkls.ctypes.data_as( ctypes.POINTER( REAL) ) 
            )

    pars.set('t_x',Ts[-1,0])
    pars.set('t_y',Ts[-1,1])
    pars.set('t_z',Ts[-1,2])

    
    import ImageD11.transform
    peaks_xyz = ImageD11.transform.compute_xyz_lab( sf, **pars.parameters )

    
    print( peaks_xyz[:,-1])

    print( "translation",pars.get('t_x'),pars.get('t_y'),pars.get('t_z'))
    origin = compute_grain_origins(c.omega, wedge = pars.get('wedge'),
                                   chi=pars.get('chi'), t_x = pars.get('t_x'),
                                   t_y = pars.get('t_y'), t_z = pars.get('t_z'))
    print (origin[:,-1], origin.shape)
    print( "I got to the end of the script")

    tth,eta = ImageD11.transform.compute_tth_eta( sf, omega=c.omega, **pars.parameters)
    kvecs = ImageD11.transform.compute_k_vectors(tth, eta, pars.get('wavelength'))
    print( kvecs.T)
    gvecs = ImageD11.transform.compute_g_from_k( kvecs, c.omega, wedge=pars.get('wedge'), 
                                                 chi=pars.get('chi'))
    print( gvecs.T)
    for i,gr in enumerate(grains):
        oldhkls = np.dot( gr.ubi, gvecs )
        print( "grain %d"%(i))
        print( oldhkls.T)
        print( hkls)


#C:\Users\wright\Programming\fable\ImageD11\test\makemap>c:\mingw64\bin\gcc ..\..\src\diffraction.c -O3 -march=native -ffast-math -static-libgcc -share d -o diffraction.dll

# C:\Users\wright\Programming\fable\ImageD11\test\makemap>python ..\..\src\cdiffraction.py test.flt map.ubi test.prm | more


