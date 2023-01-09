
from __future__ import print_function

import numpy as np
from PIL import Image
REAL = float 

def normalise( im ):
    """ Set mean to 0, std to 1 """
    avg = np.ravel(im).mean()
    std = np.ravel(im).std()
    if std > 0:
        return ( im - avg ) / std
    else:
        return im - avg 

def cross_correlate(a2, b2, fftshape):
    """ Make 2D cross correlation image """
    # FFT:
    a2 = np.fft.rfft2(a2, s = fftshape)
    b2 = np.fft.rfft2(b2, s = fftshape)
    # Cross correlate
    c = np.fft.irfft2( a2 * b2.conj() ) 
    return c

def find_offset( c ):
    """ Co-ordinates of max in a 2D cross correlation """
    flatmax = c.argmax()
    # print flatmax, c.flat[flatmax]
    dim0 = int(flatmax/c.shape[1])
    dim1 = flatmax - c.shape[1] * dim0
    # print dim0,dim1,c[dim0,dim1]
    roi = c[ dim0-6:dim0+7, dim1-6:dim1+7 ] 
    troi = roi.ravel().sum()
    x  = np.dot( roi.sum(axis=1), np.arange(dim0-6,dim0+7) ) / troi
    y  = np.dot( roi.sum(axis=0), np.arange(dim1-6,dim1+7) ) / troi
    # Average the local area ?
    return x, y
    

def register_image( rname, cname ):
    """ We take the current and place it onto reference """
    ref = normalise(np.asarray(Image.open(rname)).sum(axis=2, dtype=REAL))
    cur = normalise(np.asarray(Image.open(cname)).sum(axis=2, dtype=REAL))
    fftshape = ( ref.shape[0] + cur.shape[0] + 4,
                 ref.shape[1] + cur.shape[1] + 4 )
    cor = cross_correlate( ref, cur, fftshape)
    x, y = find_offset( cor )
    
    if False:
        from matplotlib.pylab import imshow, show, figure, colorbar, plot
        figure(1)
        imshow(ref)
        figure(2)
        imshow(cur)
        figure(3)
        imshow(cor)
        colorbar()
        print(y,x) 
        plot( [y], [x], "+", mec='g', lw=2, ms=10 )
        show()

    print(x, y)
    return cur, ref, x, y

def display_registered( current, reference, xf, yf ):
    merged = reference.copy()
    x, y = int(xf+0.5) , int( yf+0.5)
    print("cur",current.shape, "ref",reference.shape,x,y)
    print(x,x+current.shape[0] , y,y+current.shape[1]) 
    merged[ x:x+current.shape[0] , y:y+current.shape[1] ] = current
    from matplotlib.pylab import imshow, show, title
    imshow(merged)
    title("x=%f y=%f"%(xf,yf))
    show()
    


def test(filename):
    rgb = np.asarray(Image.open(filename))
    vals = 1.0*rgb[:,:,0] + rgb[:,:,1] + rgb[:,:,2]
    print(vals.shape)
    cen = vals.shape[0]/2, vals.shape[1]/2
    obj = vals[ 190:263 , 460:523 ]
    fftshape = ( vals.shape[0] + obj.shape[0] + 3, vals.shape[1] + obj.shape[1] + 3 )
    c = cross_correlate( 
            normalise(vals),
            normalise(obj),
            fftshape )
    c = c[:vals.shape[0],:vals.shape[1]]
    x, y  = find_offset(c)
    # ideal offset? position of the (0,0) pixel of obj in vals
    from matplotlib.pylab import imshow, show, figure, colorbar, plot
    figure(1)
    imshow(obj)
    figure(2)
    imshow(vals)
    figure(3)
    imshow(c)
    colorbar()
    print(y,x) 
    plot( [y], [x], "+", mec='g', lw=2, ms=10 )
    show()
    return c
    

if __name__=="__main__":
    import sys
    #    test(sys.argv[1])
    c,r,x,y = register_image( sys.argv[1], sys.argv[2] )
    display_registered( c,r, x, y )
