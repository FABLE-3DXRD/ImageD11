import math, Image, numpy

# Returns a byte-scaled image
def bytescale(data, cmin=None, cmax=None, high=255, low=0):
    if data.dtype is numpy.uint8:
        return data
    high = high - low
    if cmin is None:
        cmin = numpy.amin(numpy.ravel(data))
    if cmax is None:
        cmax = numpy.amax(numpy.ravel(data))
    scale = high *1.0 / (cmax-cmin or 1)
    bytedata = ((data*1.0-cmin)*scale + 0.4999).astype(numpy.uint8)
    return bytedata + numpy.cast[numpy.uint8](low)

def fromimage(im, flatten=0):
    """Takes a PIL image and returns a copy of the image in a numpy container.
    If the image is RGB returns a 3-dimensional array:  arr[:,:,n] is each channel

    Optional arguments:

    - flatten (0): if true, the image is flattened by calling convert('F') on
    the image object before extracting the numerical data.  This flattens the
    color layers into a single grayscale layer.  Note that the supplied image
    object is NOT modified.
    """
    assert Image.isImageType(im), "Not a PIL image."
    if flatten:
        im = im.convert('F')
    mode = im.mode
    adjust = 0
    if mode == '1':
        im = im.convert(mode='L')
        mode = 'L'
        adjust = 1
    str = im.tostring()
    type = numpy.uint8
    if mode == 'F':
        type = numpy.float32
    elif mode == 'I':
        type = numpy.uint32
    elif mode == 'I;16':
        type = numpy.uint16
    arr = numpy.fromstring(str,type)
    shape = list(im.size)
    shape.reverse()
    if mode == 'P':
        arr.shape = shape
        if im.palette.rawmode != 'RGB':
            print "Warning: Image has invalid palette."
            return arr
        pal = numpy.fromstring(im.palette.data,type)
        N = len(pal)
        pal.shape = (int(N/3.0),3)
        return arr, pal
    if mode in ['RGB','YCbCr']:
        shape += [3]
    elif mode in ['CMYK','RGBA']:
        shape += [4]
    arr.shape = shape
    if adjust:
        arr = (arr != 0)
    return arr

_errstr = "Mode is unknown or incompatible with input array shape."
def toimage(arr,high=255,low=0,cmin=None,cmax=None,pal=None,
            mode=None,channel_axis=None):
    """Takes a numpy array and returns a PIL image.  The mode of the
    PIL image depends on the array shape, the pal keyword, and the mode
    keyword.

    For 2-D arrays, if pal is a valid (N,3) byte-array giving the RGB values
    (from 0 to 255) then mode='P', otherwise mode='L', unless mode is given
    as 'F' or 'I' in which case a float and/or integer array is made

    For 3-D arrays, the channel_axis argument tells which dimension of the
      array holds the channel data.
    For 3-D arrays if one of the dimensions is 3, the mode is 'RGB'
      by default or 'YCbCr' if selected.
    if the

    The numpy array must be either 2 dimensional or 3 dimensional.
    """
    data = numpy.asarray(arr)
    if numpy.iscomplexobj(data):
        raise ValueError, "Cannot convert a complex-valued array."
    shape = list(data.shape)
    valid = len(shape)==2 or ((len(shape)==3) and \
                              ((3 in shape) or (4 in shape)))
    assert valid, "Not a suitable array shape for any mode."
    if len(shape) == 2:
        shape = (shape[1],shape[0]) # columns show up first
        if mode == 'F':
            data32 = data.astype(numpy.float32)
            image = Image.fromstring(mode,shape,data32.tostring())
            return image
        if mode in [None, 'L', 'P']:
            bytedata = bytescale(data,high=high,low=low,cmin=cmin,cmax=cmax)
            image = Image.fromstring('L',shape,bytedata.tostring())
            if pal is not None:
                image.putpalette(numpy.asarray(pal,dtype=numpy.uint8).tostring())
                # Becomes a mode='P' automagically.
            elif mode == 'P':  # default gray-scale
                pal = arange(0,256,1,dtype=numpy.uint8)[:,newaxis] * \
                      ones((3,),dtype=numpy.uint8)[newaxis,:]
                image.putpalette(numpy.asarray(pal,dtype=numpy.uint8).tostring())
            return image
        if mode == '1':  # high input gives threshold for 1
            bytedata = (data > high)
            image = Image.fromstring('1',shape,bytedata.tostring())
            return image
        if cmin is None:
            cmin = amin(ravel(data))
        if cmax is None:
            cmax = amax(ravel(data))
        data = (data*1.0 - cmin)*(high-low)/(cmax-cmin) + low
        if mode == 'I':
            data32 = data.astype(numpy.uint32)
            image = Image.fromstring(mode,shape,data32.tostring())
        else:
            raise ValueError, _errstr
        return image

    # if here then 3-d array with a 3 or a 4 in the shape length.
    # Check for 3 in datacube shape --- 'RGB' or 'YCbCr'
    if channel_axis is None:
        if (3 in shape):
            ca = numpy.flatnonzero(numpy.asarray(shape) == 3)[0]
        else:
            ca = numpy.flatnonzero(numpy.asarray(shape) == 4)
            if len(ca):
                ca = ca[0]
            else:
                raise ValueError, "Could not find channel dimension."
    else:
        ca = channel_axis

    numch = shape[ca]
    if numch not in [3,4]:
        raise ValueError, "Channel axis dimension is not valid."

    bytedata = bytescale(data,high=high,low=low,cmin=cmin,cmax=cmax)
    if ca == 2:
        strdata = bytedata.tostring()
        shape = (shape[1],shape[0])
    elif ca == 1:
        strdata = transpose(bytedata,(0,2,1)).tostring()
        shape = (shape[2],shape[0])
    elif ca == 0:
        strdata = transpose(bytedata,(1,2,0)).tostring()
        shape = (shape[2],shape[1])
    if mode is None:
        if numch == 3: mode = 'RGB'
        else: mode = 'RGBA'


    if mode not in ['RGB','RGBA','YCbCr','CMYK']:
        raise ValueError, _errstr

    if mode in ['RGB', 'YCbCr']:
        assert numch == 3, "Invalid array shape for mode."
    if mode in ['RGBA', 'CMYK']:
        assert numch == 4, "Invalid array shape for mode."

    # Here we know data and mode is coorect
    image = Image.fromstring(mode, shape, strdata)
    return image


def imrotate(arr,angle,interp='bilinear'):
    """Rotate an image counter-clockwise by angle degrees.

    Interpolation methods can be:
        'nearest' :  for nearest neighbor
        'bilinear' : for bilinear
        'cubic' or 'bicubic' : for bicubic
    """
    arr = numpy.asarray(arr)
    func = {'nearest':0,'bilinear':2,'bicubic':3,'cubic':3}
    im = toimage(arr)
    im = im.rotate(angle,resample=func[interp])
    return fromimage(im)


def radon(arr,theta=None):
    if theta is None:
        theta = numpy.mgrid[0:180]
    s = numpy.zeros((arr.shape[1],len(theta)), numpy.float)
    k = 0
    for th in theta:
        im = imrotate(arr,-th)
        s[:,k] = numpy.sum(im, axis=1)
        k += 1
    return s

def irad( ar ):
    na = ar.shape[1]
    th = numpy.arange(na)*math.pi/180.0
    n = ar.shape[0]
    inv = numpy.zeros( (n, n) , ar.dtype).ravel()
    r = numpy.arange(n) - n/2.0
    x = numpy.outer( r, numpy.ones( len(r) ))
    y = x.T
    ai = numpy.arctan2( x, y )
    angle_image = (ai*180.0/math.pi).round().astype(int)
    angle_image = numpy.where( angle_image < 0, angle_image+180, angle_image)
    angle_image = numpy.where( angle_image >179, angle_image-180, angle_image)
    radius_image = numpy.sqrt(x*x+y*y).ravel().round().astype(int)
    
    inv_goto = numpy.zeros( (n,n), ar.dtype).ravel()
    mon      = numpy.zeros( (n,n), numpy.int).ravel()
    
    rg = numpy.fft.fftshift( ar, axes=(0,))
    for i in range(na):
        data = ar[:,i]
        np = len(data)
        # Come from algorithm - compute indices in filter
        inds = numpy.compress( angle_image.ravel() == i, 
                                numpy.arange( n*n))
        for j in inds:
            k = radius_image[j] 
            if  k >= 0 and k < n:
                inv[j] = inv[j] + data[radius_image[j]]

        # Go To algorithm - compute destination in image
        
        xi = numpy.clip( r*numpy.sin( th[i] ) + n/2., 0, n).round()
        yi = numpy.clip( r*numpy.cos( th[i] ) + n/2., 0, n).round()
        inds = (xi*n + yi).astype(int)
        tmp = numpy.zeros( n*n, ar.dtype)
        tnp = numpy.zeros( n*n, numpy.int)
        data = rg[:,i]
        tmp[inds[:-1]] = data
        tnp[inds[:-1]] = 1
        inv_goto = inv_goto+tmp
        mon      = mon + tnp
    inv.shape= (n,n)
    mon.shape= (n,n)
    inv_goto.shape=(n,n)
    comb = numpy.where( mon == 0, inv, inv_goto)
    comb = comb / numpy.where( mon ==0, 1, mon)
    return comb
    pl.figure()
    pl.imshow(  numpy.log(abs(inv)))
    pl.title("come_from")
    pl.figure()
    pl.title("goto")
    pl.imshow( numpy.log(abs( inv_goto)))
    pl.figure()
    pl.title("mon")
    pl.imshow( mon)
    pl.colorbar()
    pl.figure()
    norman = inv_goto/mon
    pl.title("mon")
    pl.imshow( numpy.log(abs(norman)))
    pl.colorbar()
    pl.show()

    return inv
    imshow(numpy.log((inv*inv.conj()).real))
    show()
    sys.exit()


def fourier_radial( sinogram , theta=None ):
    """
    sinogram is from the data
        dimensions [npixels, nangles]
    th is the list of angles in degrees of the projections
        assumed 1 degree steps otherwise
    returns the radon transform

    Currently is not interpolating in fourier space
    """
    if theta is not None:
        assert len(theta) == sinogram.shape[1]
        th = numpy.array(theta)*numpy.pi/180.0
    else:
        th = numpy.array( [ i*math.pi/180.0 for i in range(sinogram.shape[1]) ] )
    # sinogram has len(512)
    print "Sinogram shape", sinogram.shape
    ar = numpy.fft.rfft(sinogram, axis=0)
    print "FT sinogram shape", ar.shape
    # This has len(257)
    # Vertical and horizontal pixels ... to use irfft, rfft
    nv = (ar.shape[0]-1)*2
    nh = ar.shape[0]
    fim = numpy.zeros( ( nv, nh), ar.dtype )
    nim = numpy.zeros( ( nv, nh), numpy.int )
    print "Target ft image shape",fim.shape
    # Use a GOTO algorithm to begin with
    nlost = {}
    for fproj, theta in zip(ar.T, th):
        nlost[theta]=0
        # fproj is the fourier transform of a projection
        # theta is the projection angle
        cth = math.cos(theta)
        sth = math.sin(theta)
        for r in range(len(fproj)):
            i = int(round( r * cth ))
            j = int(round( r * sth ))
            f = fproj[r]
            if j < 0:
                i = -i
                j = -j
                f = f.conj()
            # Assume the numpy array indexing is OK for negatives
            # this works for th 0->180
            if i >= -nv/2 and i <= nv/2 and j >= 0 and j < nh:
                fim[i,j] += f
                nim[i,j] += 1
            else:
                nlost[theta] += 1
                print theta*180./math.pi, i, j, nv/2, nh
    for theta in th:
        if nlost[theta] > 0:
            print theta*180./math.pi, nlost[theta],"   ",
    print
    return fim/nim.clip(1,len(th)*nv*nh)




#if __name__=="__main__":
if 1:
    im = Image.open( "phantom3.gif" )
    im = im.transform( (512,512), Image.AFFINE, ( 1,0,-128, 0,1,-128) )
    a = numpy.asarray(im).astype(numpy.float32)
    theta = range(0,110,5)
    r = radon( a , theta = theta)
    #    r[:,angle]
    print r.shape
    r = numpy.concatenate( (  r[-r.shape[0]/2:,:], r[:r.shape[0]/2,:], ), axis=0)
    # shift to put the centre of the sinogram at the origin(!)
    f = fourier_radial( r, theta )
    import pylab as pl
    pl.figure(1)
    pl.clf()
    pl.imshow(a)
    pl.figure(2)
    pl.clf()
    pl.imshow(r)
    pl.figure(3)
    pl.clf()
    pl.imshow( abs(f))
    ans = numpy.fft.irfft2( f )
    pl.figure(4)
    pl.clf()
    pl.imshow( numpy.fft.fftshift(ans )  )
    pl.show()



