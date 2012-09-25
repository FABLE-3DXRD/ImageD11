

import ctypes

dll = ctypes.CDLL("diffraction.dll")

dll.choosegrains.argtypes = [
#void choosegrains( vector XL[], real omega[], 
   ctypes.POINTER(ctypes.c_float),
   ctypes.POINTER(ctypes.c_float),
#        real wedge, real chi, real wvln,
   ctypes.c_float, ctypes.c_float, ctypes.c_float,
#        rmatrix UB[], rmatrix UBI[], vector T[], 
   ctypes.POINTER(ctypes.c_float),
   ctypes.POINTER(ctypes.c_float),
   ctypes.POINTER(ctypes.c_float),
#        int ngrains, int npks,
   ctypes.c_int, ctypes.c_int,
#        int labels[], vector hkls[]){ 
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_float)
        ]

def flatfloat(a, shape=None):
    if shape is None:
        s = a.shape
    else:
        s=shape
    flat = np.array( np.array(a).ravel(), np.float32 )
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
    
    wedge = 0.01#float(pars.get("wedge"))
    chi   = 0.02#float(pars.get("chi"))
    wvln  = float(pars.get("wavelength"))
    
    c = columnfile( colfile )

    XL = flatfloat( compute_xyz_lab( [ c.xc, c.yc] , **pars.parameters ).T )
    om = flatfloat( c.omega )


#    XL[0,:]=[41,42,43]

    grains = read_grain_file( grainfile )
    UBIs = flatfloat( [ g.ubi for g in grains ], shape=(len(grains),3,3) )
    UBs =  flatfloat( [ np.linalg.inv(g.ubi) for g in grains ], shape=(len(grains),3,3) )
    Ts  =  flatfloat( [ g.translation for g in grains ], shape=(len(grains),3))

    labels = np.zeros( c.nrows, np.int )
    hkls = np.zeros( ( c.nrows, 3 ), np.float32 )

    dll.choosegrains( 
            XL.ctypes.data_as( ctypes.POINTER( ctypes.c_float) ),
            om.ctypes.data_as( ctypes.POINTER( ctypes.c_float) ),
            wedge, chi, wvln,
            UBs.ctypes.data_as( ctypes.POINTER( ctypes.c_float) ),
            UBIs.ctypes.data_as( ctypes.POINTER( ctypes.c_float) ),
            Ts.ctypes.data_as( ctypes.POINTER( ctypes.c_float) ),
            len(grains), c.nrows,
            labels.ctypes.data_as( ctypes.POINTER( ctypes.c_int) ),
            hkls.ctypes.data_as( ctypes.POINTER( ctypes.c_float) ) 
            )
    for i in range(20):
        print labels[i],hkls[i]


#C:\Users\wright\Programming\fable\ImageD11\test\makemap>c:\mingw64\bin\gcc ..\..\src\diffraction.c -O3 -march=native -ffast-math -static-libgcc -share d -o diffraction.dll

# C:\Users\wright\Programming\fable\ImageD11\test\makemap>python ..\..\src\cdiffraction.py test.flt map.ubi test.prm | more


