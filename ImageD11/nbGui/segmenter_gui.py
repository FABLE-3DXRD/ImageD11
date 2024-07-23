


import os
import numpy as np
import matplotlib.pyplot as plt
from ipywidgets import interact, interactive, widgets, fixed, Layout
import h5py
import ImageD11.sinograms.lima_segmenter
import ImageD11.sparseframe

class SegmenterGui:
    
    """ UI for a jupyter notebook to set the segmentation parameters
    From @jadball notebook, @jonwright refactored to put in a python file
    """
    
    def __init__(self, dset, cut=1, pixels_in_spot=3, howmany=100000, counter="_roi1"):
        self.dset = dset
        self.options = { "maskfile" : dset.maskfile ,
                        "cut" : cut,
                        "pixels_in_spot" : pixels_in_spot,
                        "howmany": howmany }
        self.scan, self.idx = self.chooseframe(counter)
        cut_slider = widgets.IntSlider(value=cut, min=0, max=200, step=1, description='Cut:')
        pixels_in_spot_slider = widgets.IntSlider(value=pixels_in_spot, min=0, max=20, step=1, description='Pixels in Spot:')
        howmany_slider = widgets.IntSlider(value=np.log10(howmany), min=1, max=15, step=1, description='log(howmany):')
        self.display()
        self.widget = widgets.interactive(self.update_image, cut=cut_slider, pixels_in_spot=pixels_in_spot_slider, howmany=howmany_slider)
        display( self.widget )
        
    def chooseframe(self, counter):
        ds = self.dset
        # Locate a busy image to look at
        with h5py.File(ds.masterfile,'r') as hin:
            scan = ds.scans[len(ds.scans)//2]
            if scan.find("::"): # 1.1::[10000:12000]  etc
                lo, hi = [int(v) for v in scan[:-1].split("[")[1].split(":")]
                scan = scan.split("::")[0]
                roi1 = hin[scan]['measurement'][f'{ds.detector}{counter}'][lo:hi]
                idx = np.argmax(roi1) + lo
            else: # "1.1"
                roi1 = hin[scan]['measurement'][f'{ds.detector}{counter}'][:]
                idx = np.argmax(roi1) + lo
        print("Using frame", idx, "from scan", scan)
        return scan, idx
        
    def segment_frame(self):
        """
        ds = ImageD11.sinograms.dataset object
        options = dict to pass to ImageD11.sinograms.lima_segmenter.SegmenterOptions
        image_file_num = which scan0123/eiger_0000.h5 to look at
        frame_num = which frame in the file
        """
        ds = self.dset
        opts = ImageD11.sinograms.lima_segmenter.OPTIONS = ImageD11.sinograms.lima_segmenter.SegmenterOptions(**self.options)
        opts.setup()
        with h5py.File( ds.masterfile, 'r' ) as hin:
            frms = hin[self.scan]['measurement'][ds.detector]
            for i, spf in enumerate( ImageD11.sinograms.lima_segmenter.reader( frms, opts.mask, opts.cut, start = self.idx ) ): 
                if spf is None:
                    print("Warning: no pixels found",scan,i,self.idx)
                ref = frms[i+self.idx]
                break    
        if spf is None:
            spi = np.zeros_like(ref)
            npeaks = 0
        else:
            spi = spf.to_dense('intensity')
            npeaks = ImageD11.sparseframe.sparse_localmax( spf )
        if opts.mask is not None:
            ref = ref * opts.mask
            spi = spi * opts.mask    
        return ref, spi, npeaks

    def display(self):
        # Display the image initially
        self.fig, self.axs = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(16, 9))
        raw_image, segmented_image, nblobs = self.segment_frame()
        self.im1 = self.axs[0].imshow(raw_image, cmap="viridis", norm='log', vmin=1, vmax=1000, interpolation="nearest")
        self.im2 = self.axs[1].imshow(segmented_image, cmap="viridis", norm='log', vmin=0.5, vmax=1000, interpolation="nearest")
        self.axs[0].set_title("Raw image")
        self.axs[1].set_title("Segmented image")
        self.fig.suptitle(f"{nblobs} peaks found\n cut={self.options['cut']}, pixels_in_spot={self.options['pixels_in_spot']}, howmany={self.options['howmany']}")
        plt.show()


    def update_image(self, cut, pixels_in_spot, howmany):
        howmany_exp = 10**howmany
        self.options["cut"] = cut
        self.options["pixels_in_spot"] = pixels_in_spot
        self.options["howmany"] = howmany_exp
        raw_image, segmented_image, nblobs = self.segment_frame()
        self.im1.set_data(raw_image)
        self.im2.set_data(segmented_image)
        self.fig.suptitle(f"{nblobs} peaks found\n cut={cut}, pixels_in_spot={pixels_in_spot}, howmany={howmany_exp}")
        plt.draw()


    def getopts(self):
        opts = { name: self.options[name] for name in ('cut','pixels_in_spot','howmany') }
        print("options = ",repr(opts))
        return opts

