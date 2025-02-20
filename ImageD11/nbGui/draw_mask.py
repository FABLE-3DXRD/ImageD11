import numpy as np

from matplotlib import pyplot as plt
from matplotlib.widgets import PolygonSelector
from matplotlib.path import Path

from skimage.filters import threshold_otsu
from skimage.morphology import convex_hull_image


def plot_mask_result(image, mask):
    fig, axs = plt.subplots(1, 3, sharex=True, sharey=True, layout='constrained', figsize=(10, 4))
    axs[0].imshow(image, origin='lower')
    axs[1].imshow(mask, origin='lower')
    axs[2].imshow(np.where(mask, image, 0), origin='lower')
    axs[0].set_title('Reconstruction')
    axs[1].set_title('Mask')
    axs[2].set_title('Masked reconstruction')
    fig.supxlabel("<-- Sample Y axis")
    fig.supylabel("Sample x axis (Beam) >")
    plt.show()


class InteractiveMask:
    def __init__(self, image):
        self.image = image
        self.polygon_points = None  # Store polygon points
        self.fig, self.ax = plt.subplots(layout='constrained')
        self.ax.imshow(image, origin='lower')
        self.selector = PolygonSelector(self.ax, self.on_select)
        
        # Add instruction text
        self.ax.set_title("Draw mask and press Enter to finish")
        
        # Connect event to close plot on Enter key
        self.fig.canvas.mpl_connect("key_press_event", self.on_key)
        self.fig.supxlabel("<-- Sample Y axis")
        self.fig.supylabel("Sample x axis (Beam) >")
        
        plt.show(block=False)  # Non-blocking for Jupyter

    def on_select(self, verts):
        """Callback function for polygon selection."""
        self.polygon_points = np.array(verts, dtype=np.int32)

    def on_key(self, event):
        """Closes the figure when Enter is pressed."""
        if event.key == "enter":
            plt.close(self.fig)  # Close the figure

    def get_mask(self, doplot=True):
        """Waits for user to draw polygon and then returns the mask."""
        # Generate mask after plot is closed
        mask = np.zeros_like(self.image, dtype=bool)

        if self.polygon_points is not None and len(self.polygon_points) > 0:
            yy, xx = np.meshgrid(np.arange(self.image.shape[0]), np.arange(self.image.shape[1]), indexing="ij")
            polygon_path = Path(self.polygon_points)
            mask = polygon_path.contains_points(np.vstack((xx.ravel(), yy.ravel())).T).reshape(self.image.shape).astype(bool)
        
        if doplot:
            plot_mask_result(self.image, mask)
        
        return mask

    
def threshold_mask(image, manual_threshold=None, doplot=True):
    if manual_threshold is None:
        thresh = threshold_otsu(image)
    else:
        thresh = manual_threshold

    binary = image > thresh
    chull = convex_hull_image(binary)
    whole_sample_mask = chull
    
    if doplot:
        plot_mask_result(image, whole_sample_mask)
    
    return whole_sample_mask