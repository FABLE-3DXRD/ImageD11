
import h5py
import numpy as np
from ImageD11 import sparseframe, cImageD11

def moments( bfrm ):
    I = bfrm[:, cImageD11.s2D_I]
    s = bfrm[:, cImageD11.s2D_sI] / I
    f = bfrm[:, cImageD11.s2D_fI] / I
    n = bfrm[:, cImageD11.s2D_1]
    a = I / n
    return s, f, a, n

def newpks(hf, scans=None, npixels=1, monitor=None, monval=1e3, ):
    """ do a peak search in a sparse frame series """
    titles = 's_raw f_raw avg_intensity Number_of_pixels sum_intensity omega dty'.split()
    c = {}
    for name in titles:
        c[name] = []
    # using the python file open overcomes some threading crap
    with h5py.File(open(hf,"rb"),'r') as hin:
        # scan numbers
        if scans is None:
            scans = list(hin['/'])
        for scan in scans:
            gin = hin[scan]
            shape = gin.attrs['shape0'], gin.attrs['shape1']
            # import pdb; pdb.set_trace()
            omega = gin['measurement/rot'][:]
            difty = gin['instrument/positioners/dty'][()]
            row = gin['row'][()]
            col = gin['col'][()]
            sig = gin['intensity'][()]
            if monitor is not None:
                mon = monval/gin[monitor][()]
            ipt = np.cumsum( gin['nnz'][:] )
            iprev = 0
            for k,nnz in enumerate( gin['nnz'][()] ) :
                inext = iprev + nnz
                if nnz == 0:
                    continue
                f = sparseframe.sparse_frame( row[iprev:inext],
                                              col[iprev:inext], shape )
                f.set_pixels("intensity", sig[iprev:inext] )
                sparseframe.sparse_connected_pixels(f, threshold=0.1)
                pks = sparseframe.sparse_moments(f, "intensity", "connectedpixels")
                # sparseframe.sparse_localmax(f)
                # pks = sparseframe.sparse_moments(f, "intensity", "localmax")
                s,f,a,n = moments( pks )
                m = n > npixels
                c['s_raw'].append( s[m] )
                c['f_raw'].append( f[m] )
                if monitor is not None:
                    c['avg_intensity'].append( a[m]*mon[k] )
                    c['sum_intensity'].append( a[m]*n[m]*mon[k] )
                else:
                    c['avg_intensity'].append( a[m] )
                    c['sum_intensity'].append( a[m]*n[m] )                
                c['Number_of_pixels'].append( n[m] )
                npk = m.sum()
                c['omega' ].append( np.full( npk , omega[k] ) )
                c['dty'].append( np.full( npk , difty ) )
                iprev = inext
    for t in titles:
        c[t] = np.concatenate( c[t] )
    return c

def pfun( *args ):
#    print(args[0])
    s, h5name, npixels = args[0]
    return (s, newpks( h5name, [s,], npixels = npixels))  
