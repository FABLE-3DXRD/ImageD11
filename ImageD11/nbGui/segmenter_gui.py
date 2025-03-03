import numpy as np
import matplotlib.pyplot as plt
from ipywidgets import interact, interactive, widgets, fixed, Layout
from IPython.display import display
import h5py
import ImageD11.sinograms.lima_segmenter
import ImageD11.sparseframe
from ImageD11.nbGui.segmenter_gui_utils import * 

class SegmenterGui:
    """ UI for a jupyter notebook to set the segmentation parameters
    From @jadball notebook, @jonwright refactored to put in a python file
    """

    def __init__(self, dset, counter="_roi1", scan=None, frame=None, cut=1, pixels_in_spot=3, howmany=100000, ):
        self.dset = dset
        self.fig = None
        self.options = {"maskfile": dset.maskfile,
                        "cut": cut,
                        "pixels_in_spot": pixels_in_spot,
                        "howmany": howmany}
        self.scan = scan
        self.idx = frame
        self.chooseframe(counter)
        cut_slider = widgets.IntSlider(value=cut, min=0, max=200, step=1, description='Cut:')
        pixels_in_spot_slider = widgets.IntSlider(value=pixels_in_spot, min=0, max=20, step=1,
                                                  description='Pixels in Spot:')
        howmany_slider = widgets.IntSlider(value=np.log10(howmany), min=1, max=15, step=1, description='log(howmany):')
        self.widget = widgets.interactive(self.update_image, cut=cut_slider, pixels_in_spot=pixels_in_spot_slider,
                                          howmany=howmany_slider)
        display(self.widget)

    def chooseframe(self, counter):
        # Use the base class method without fetching the raw image
        self.scan, self.idx, _ = chooseframe(self.dset, self.scan, self.idx, counter, fetch_raw_image=False)

    def segment_frame(self):
        """
        ds = ImageD11.sinograms.dataset object
        options = dict to pass to ImageD11.sinograms.lima_segmenter.SegmenterOptions
        image_file_num = which scan0123/eiger_0000.h5 to look at
        frame_num = which frame in the file
        """
        ds = self.dset
        opts = ImageD11.sinograms.lima_segmenter.OPTIONS = ImageD11.sinograms.lima_segmenter.SegmenterOptions \
            (**self.options)
        opts.setup()
        with h5py.File(ds.masterfile, 'r') as hin:
            frms = hin[self.scan]['measurement'][ds.detector]
            for i, spf in enumerate(
                    ImageD11.sinograms.lima_segmenter.reader(frms, opts.mask, opts.cut, start=self.idx)):
                if spf is None:
                    print("Warning: no pixels found", self.scan, i, self.idx)
                ref = frms[i + self.idx]
                break
        if spf is None:
            spi = np.zeros_like(ref)
            npeaks = 0
        else:
            spi = spf.to_dense('intensity')
            npeaks = ImageD11.sparseframe.sparse_localmax(spf)
        if opts.mask is not None:
            ref = ref * opts.mask
            spi = spi * opts.mask
        return ref, spi, npeaks

    def display(self):
        # Display the image initially
        self.fig, self.axs = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(10, 5), constrained_layout=True)
        self.im1 = self.axs[0].imshow(self.raw_image, cmap="viridis", norm='log', vmin=0.5, vmax=1000,
                                      interpolation="nearest")
        self.im2 = self.axs[1].imshow(self.segmented_image, cmap="viridis", norm='log', vmin=0.5, vmax=1000,
                                      interpolation="nearest")
        self.axs[0].set_title("Raw image")
        self.axs[1].set_title("Segmented image")
        self.fig.suptitle("%d peaks found with cut=%d, pixels_in_spot=%d, howmany=%d " %
                          (self.nblobs, self.options['cut'], self.options['pixels_in_spot'], self.options['howmany']))
        # self.fig.show()

    def update_image(self, cut, pixels_in_spot, howmany):
        howmany_exp = 10 ** howmany
        self.options["cut"] = cut
        self.options["pixels_in_spot"] = pixels_in_spot
        self.options["howmany"] = howmany_exp
        self.raw_image, self.segmented_image, self.nblobs = self.segment_frame()
        if self.fig is None:
            self.display()
        self.im1.set_data(self.raw_image)
        self.im2.set_data(self.segmented_image)
        self.fig.suptitle("%d peaks found with cut=%d, pixels_in_spot=%d, howmany=%d " %
                          (self.nblobs, cut, pixels_in_spot, howmany))
        self.fig.canvas.draw()

    def getopts(self):
        opts = {name: self.options[name] for name in ('cut', 'pixels_in_spot', 'howmany')}
        print("options = ", repr(opts))
        return opts


class FrelonSegmenterGui:
    """ UI for a jupyter notebook to set the segmentation parameters
    From @jadball notebook, @jadball refactored to put in a python file
    """

    def __init__(self, dset, worker_func, process_func, counter="_roi1", scan=None, frame=None, **options):
        self.dset = dset
        self.worker_func = worker_func
        self.process_func = process_func
        self.fig = None
        self.options = options
        self.scan = scan
        self.idx = frame
        self.chooseframe(counter)
        thresh_slider = widgets.IntSlider(value=options["threshold"], min=0, max=10000, step=1,
                                          description='Threshold:')
        smsig_slider = widgets.FloatSlider(value=options["smoothsigma"], min=0.0, max=3.0, step=0.05,
                                           description='Smoothsigma:')
        bgc_slider = widgets.FloatSlider(value=options["bgc"], min=0.0, max=1.0, step=0.05, description='bgc:')
        minpx_slider = widgets.IntSlider(value=options["minpx"], min=0, max=50, step=1, description='minpx:')
        mofft_slider = widgets.IntSlider(value=options["m_offset_thresh"], min=0, max=200, step=1,
                                         description='m_offset_thresh:')
        mratt_slider = widgets.IntSlider(value=options["m_ratio_thresh"], min=0, max=2000, step=1,
                                         description='m_ratio_thresh:')
        self.widget = widgets.interactive(self.update_image,
                                          threshold=thresh_slider,
                                          smoothsigma=smsig_slider,
                                          bgc=bgc_slider,
                                          minpx=minpx_slider,
                                          m_offset_thresh=mofft_slider,
                                          m_ratio_thresh=mratt_slider)
        display(self.widget)

    def chooseframe(self, counter):
        self.scan, self.idx, self.raw_image = chooseframe(self.dset, self.scan, self.idx, counter, fetch_raw_image=True)
        self.raw_image = self.raw_image.astype("uint16")

    def segment_frame(self):
        image_worker = self.worker_func(**self.options)
        goodpeaks = image_worker.peaksearch(img=self.raw_image, omega=0)
        fc, sc = goodpeaks[:, 23:25].T  # 23 and 24 are the columns for fc and sc from blob properties
        return image_worker, fc, sc, len(fc)

    def display(self):
        self.fig, self.axs = plt.subplots(1, 3, sharex=True, sharey=True, figsize=(10, 3.5), layout="constrained")
        self.im1 = self.axs[0].imshow(self.raw_image, norm='log', vmin=100, vmax=1000)
        self.axs[0].set_title("Original image")
        self.im2 = self.axs[1].imshow(self.smoothed_image, cmap="viridis", norm='log', vmin=0.5, vmax=1000,
                                      interpolation="nearest")
        self.axs[1].set_title("Background corrected")
        self.im3 = self.axs[2].imshow(self.smoothed_image, cmap="viridis", norm='log', vmin=0.5, vmax=1000,
                                      interpolation="nearest")
        self.axs[2].set_title(str(self.npeaks) + " peaks")
        self.sc1, = self.axs[2].plot(self.fc, self.sc, marker='+', c="r", ls="")
        self.axs[2].set_aspect(1)

    def update_image(self, threshold, smoothsigma, bgc, minpx, m_offset_thresh, m_ratio_thresh):
        self.options["threshold"] = threshold
        self.options["smoothsigma"] = smoothsigma
        self.options["bgc"] = bgc
        self.options["minpx"] = minpx
        self.options["m_offset_thresh"] = m_offset_thresh
        self.options["m_ratio_thresh"] = m_ratio_thresh

        image_worker, self.fc, self.sc, self.npeaks = self.segment_frame()
        self.smoothed_image = image_worker.smoothed
        if self.fig is None:
            self.display()

        self.im1.set_data(self.raw_image)
        self.im2.set_data(self.smoothed_image)
        self.im3.set_data(self.smoothed_image)
        self.sc1.set_data(self.fc, self.sc)
        self.axs[2].set_title(str(self.npeaks) + " peaks")
        self.fig.canvas.draw()

    def getopts(self):
        opts = {name: self.options[name] for name in
                ('bgfile', 'maskfile', 'flatfile', 'darkfile', 'threshold', 'smoothsigma', 'bgc', 'minpx', 'm_offset_thresh', 'm_ratio_thresh') if name in self.options.keys()}
        print("options = ", repr(opts))
        return opts
