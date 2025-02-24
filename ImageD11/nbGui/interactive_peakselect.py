import numpy as np
from matplotlib import pyplot as plt

import ImageD11.sinograms.geometry


class FriedelPairSelector:
    def __init__(self, cring, cpk, y0):
        self.cring = cring
        self.cpk = cpk
        self.y0 = y0
        
        self.p1 = None  # Index of first peak
        self.p2 = None  # Index of second peak
        self.xy = None  # Fitted position
        self.peak_ycalc = None  # Computed dty values
        
        self._setup_plot()
        
    def _setup_plot(self):
        plt.ion()  # Enable interactive mode
        self.fig, self.ax = plt.subplots(1, 2, figsize=(10, 7), sharex=True, sharey=True, layout='constrained')
        
        # Left plot (selectable peaks)
        self.ax[0].scatter(
            self.cpk.omega, self.cpk.dty,
            c=np.log(self.cpk.sum_intensity) * np.sign(self.cpk.yl),
            s=4, cmap='RdBu', picker=True, pickradius=5
        )
        self.selected, = self.ax[0].plot([], [], "o", mfc='none')
        
        # Right plot (reference peaks)
        self.ax[1].scatter(
            self.cring.omega, self.cring.dty,
            c=np.log(self.cring.sum_intensity) * np.sign(self.cring.yl),
            s=4, cmap='RdBu'
        )
        
        # Fit lines
        self.om = np.linspace(self.cring.omega.min(), self.cring.omega.max(), 90)
        self.fitlinel, = self.ax[0].plot(self.om, np.zeros_like(self.om), '-')
        self.fitliner, = self.ax[1].plot(self.om, np.zeros_like(self.om), '-')
        
        # Titles
        self.ax[0].set(title="Peaks from the ring you selected")
        self.ax[1].set(title="All peaks")
        self.fig.supxlabel(r"$\omega~(\degree)$")
        self.fig.supylabel("dty")
        self.fig.suptitle("Pick a peak on the left plot")
        
        # Connect interactive event
        self.fig.canvas.mpl_connect('pick_event', self._onpick)
        plt.show()
        
    def _find_pair(self, idx):
        """Find the Friedel pair for a selected peak."""
        dg2 = (self.cpk.gx[idx] + self.cpk.gx) ** 2 + \
              (self.cpk.gy[idx] + self.cpk.gy) ** 2 + \
              (self.cpk.gz[idx] + self.cpk.gz) ** 2
        return np.argmin(dg2)
        
    def _fit_pair(self, i, j):
        """Fit sine wave to the selected peak pair."""
        self.p1, self.p2 = i, j
        
        ci, cj = np.cos(np.radians([self.cpk.omega[i], self.cpk.omega[j]]))
        si, sj = np.sin(np.radians([self.cpk.omega[i], self.cpk.omega[j]]))
        yi, yj = self.cpk.dty[i] - self.y0, self.cpk.dty[j] - self.y0
        
        R = np.array([[-si, -ci], [-sj, -cj]])
        self.xy = np.linalg.inv(R).dot([yi, yj])
        sx, sy = self.xy
        
        ycalc = ImageD11.sinograms.geometry.x_y_y0_omega_to_dty(self.om, sx, sy, self.y0)
        self.fitlinel.set_ydata(ycalc)
        self.fitliner.set_ydata(ycalc)
        self.peak_ycalc = ImageD11.sinograms.geometry.x_y_y0_omega_to_dty(self.cring.omega, sx, sy, self.y0)
        
    def _onpick(self, event):
        """Handle peak selection by user."""
        ind = event.ind[-1]
        pair = self._find_pair(ind)
        self._fit_pair(ind, pair)
        
        self.selected.set_xdata([self.cpk.omega[ind], self.cpk.omega[pair]])
        self.selected.set_ydata([self.cpk.dty[ind], self.cpk.dty[pair]])
        self.ax[0].set(title='Selected Peaks: ' + str(ind) + ' ' + str(pair))
        self.fig.canvas.draw_idle()
        
    def get_selection(self):
        """Return the selected indices and fitted values."""
        print("Selected points are",self.p1,self.p2)
        print("Omega, dty")
        print(self.cpk.omega[self.p1],self.cpk.dty[self.p1])
        print(self.cpk.omega[self.p2],self.cpk.dty[self.p2])
        return self.p1, self.p2, self.xy, self.peak_ycalc
