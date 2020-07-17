
from __future__ import print_function

import fabio, pyFAI, sys, numpy as np, scipy.ndimage, os

class imagecleaner(object):
    def __init__(self, poni, npchi, npfilt, nptth,
                 maskfile=None, dims=(2048,2048) ):
        self.npchi = npchi
        self.npfilt = npfilt
        self.nptth = nptth
        
        self.ai = pyFAI.load( poni )
        self.ai.setChiDiscAtZero() # Does not help, wanted pi/2
        
        if maskfile is not None:
            self.immask = 1-fabio.open(maskfile).data
        else:
            self.immask = np.ones( dims, np.float32)
        tmps, x, y = self.ai.integrate2d(
            self.immask, self.nptth, self.npchi, unit='2th_rad',
            azimuth_range=(0,360),
            correctSolidAngle=False,
            polarization_factor=None)
        self.mask = (tmps > 0.8).T

        tthstep = (x[-1]-x[0])/len(x)
        self.tthbin = np.round(
            (self.ai.twoThetaArray() - x[0]) / tthstep ).astype(int)
        chid = np.degrees( self.ai.chiArray() )
        chistep = 360.0/(len(y)-1)
        self.chibin = np.round( (chid ) / chistep ).astype(int)
        

    def clean(self, img):
        import time
        start = time.time()
        surf, x, y = self.ai.integrate2d( img, nptth, npchi, unit='2th_rad',
                                          azimuth_range=(0,360),
                                          correctSolidAngle=False,
                                          polarization_factor=None)
        print("pyFAI",time.time()-start, end=' ')
        start = time.time()
        for i in range(self.mask.shape[0]):
            sel = np.compress( self.mask[i] , surf[:,i] )
            if len(sel)>2:
                val = np.median( sel )
                surf[~self.mask[i],i]=val
            
        print("fill",time.time()-start, end=' ')
        start = time.time()
#        powder = scipy.ndimage.minimum_filter1d( surf, npfilt, axis=0 )
#        powder = scipy.ndimage.median_filter( surf, (npfilt,1) )
        powder = scipy.ndimage.percentile_filter( surf.T.copy(), percentile=10,
                                                  size=(1,npfilt))
        print("percentile",time.time()-start, end=' ')
        start= time.time()
        powderimage = powder[ self.tthbin, self.chibin ]
        print("splat",time.time()-start)
        return powderimage
        


if __name__=="__main__":


    poni = sys.argv[1]
    maskfile = sys.argv[2]
    dark = sys.argv[3]

    imgfiles = sys.argv[4:] 
    inputfolder = os.path.split(sys.argv[4])[0]
    outputfolder = inputfolder+"_edna"
    
    import glob, os
    
    assert os.path.isdir( inputfolder )
    assert os.path.isdir( outputfolder )
    
    npchi  = 3600 # int(sys.argv[3])
    npfilt =  300 # int(sys.argv[4])
    nptth  = 2000
    
    dark = fabio.open(dark).data
    dims = dark.shape

    o = imagecleaner( poni, npchi, npfilt, nptth,
                      maskfile=maskfile, dims=dims )
    
    mask = o.immask
    showplot= True
    
    for imgfile in imgfiles:# 
        img = fabio.open( imgfile ).data.astype(np.float32) - dark
        powd = o.clean( img ) 
        outname = os.path.join( outputfolder,
                                os.path.split( imgfile )[-1].replace("edf","pow" ) )
        tth, powI = o.ai.integrate1d( powd, 2000, mask= 1-o.immask, unit='2th_deg',
                                      correctSolidAngle=True,
                                      filename=outname )
        pks = img - powd
        outname = os.path.join( outputfolder,
                                os.path.split( imgfile )[-1].replace("edf","cor" ) )
        print("Treated",imgfile, outname)
        fabio.edfimage.edfimage( pks.clip(0,65535).astype(np.uint16) ).write( outname )
        os.system("gzip --fast %s"%(outname))
        if showplot:
            import pylab as pl
            pl.figure()
            pl.subplot(231)
            pl.imshow( np.log(img.astype(np.float32)*mask),vmin=np.log(10),
                       vmax=np.log(1000))
            pl.title("Total")
            pl.colorbar()

            pl.subplot(232)
            pl.imshow( np.log((img-powd+10)*mask),vmin=np.log(10),vmax=np.log(1000))
            pl.title( "Peaks only")
            pl.colorbar()

            pl.subplot(233)
            pl.imshow( np.log(powd*mask+10 ),vmin=np.log(10) )
            pl.title("Powder Background")
            pl.colorbar()

    
            pl.subplot(212)
            tth, I = o.ai.integrate1d( img, 2000, mask= 1-o.immask, unit='2th_deg',
                                      correctSolidAngle=True,
                                      filename=outname )
            pl.plot(tth, I, "-", label='all')
            pl.plot(tth, powI,"-",label='powder')
            pl.ylim( 10,1000 )
            pl.xlim( tth[0], tth[-1])
            #pl.semilogy()
            pl.legend()
            pl.show()

            showplot=False

