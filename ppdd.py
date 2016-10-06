import sys
import os
import numpy as np
from scipy import signal
import cunwrap
import abel
from matplotlib import pyplot as plt, patches

import fittools

class PPDD(object):
    """
    Python Plasma Density Diagnostics(PPDD) is the main class to read input data, perform abel transfrom and output transfrom result.
    A run should consists of methods in the following order:
    readfile
    find_peaks
    filt_move
    find_symmetry_axis
    abel
    """
    def __init__(self, xmin = 0, xmax = 800, ymin = 400, ymax = 600, xband = 0.01, yband = 0.1, symin = 50, symax = 150, method = 'hansenlaw', **kwargs):
        self.guess = fittools.Guess(**kwargs);
        self.xmin = xmin        #crop the region [ymin:ymax, xmin:xmax] from raw input data
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.xband = xband      #half passbands for filter in x and y direction
        self.yband = yband
        self.symin = symin      #limits the symmetry axis finding to [symin, symax]
        self.symax = symax
        self.method = method
        self.peak_fitted = False

        self.abel_methods = {
            "hansenlaw": self.abel_hansenlaw,
            "onion_bordas": self.abel_onion_bordas,
            "basex": self.abel_basex
        }

    def readfile(self, filename):
        """
        Read input file filename
        """
        self.rawdata = np.loadtxt(filename, dtype=int)
        self.peak_fitted = False

    def find_peaks(self):
        """
        Find the three peaks in the frequency spectrum. This procedure includes cropping the appropriate region.
        """
        if not self.peak_fitted :  #if already fitted peaks, skip to speed up
            #loaded new file or region has changed
            self.xy2d = self.rawdata[self.ymin:self.ymax, self.xmin:self.xmax]
            #create the shifted amplitude spectrum to fit
            XYf2d = np.fft.fftn(self.xy2d)
            self.XYf2d_shifted = np.abs(np.fft.fftshift(XYf2d))                #shift frequency of (0,0) to the center

            #try hot start
            try :
                self.find_peaks_hot_start()
            except RuntimeError :
                #hot start failed, fall back to cold start
                self.find_peaks_cold_start()    #if this fail, a RuntimeError will be raised
            self.peak_fitted = True

    def find_peaks_hot_start(self):
        """
        Find the three peaks in the frequency spectrum. 
        """
        self.fx, self.fy, newguess = fittools.find_peaks(self.XYf2d_shifted, self.guess) 
        self.guess = newguess

    def find_peaks_cold_start(self):
        """
        A cold start without relying on provided initial guess to fit the three peaks. This might be slow.
        """
        length_x = self.XYf2d_shifted.shape[1]
        length_y = self.XYf2d_shifted.shape[0]
        dXf = 1/length_x
        dYf = 1/length_y

        y_x0 = self.XYf2d_shifted[:, XYf2d_shifted.shape[1]//2]         #the x center line, gg if the main peak isn't there!
        popt_y = fittools.fit_gaussian(y_x0, -np.inf, np.inf)           #fit the 1D center line to get a good starting point for sigam
        #get the y sigam
        sigma_y0 = popt_y[2]

        x_sum = np.sum(self.XYf2d_shifted, axis = 0)                    #the projection on x axis     
        peaks = signal.find_peaks_cwt(x_sum, np.arange(1, 0.1*length_x))        #TODO further polish wavelet coefficients 
        acceptance = 0.1                                                #the left and right peaks should be symmetric, this is the maximum accepted difference ratio 
        center_index = fittool.find_nearest(peaks, length_x//2)
        left_peak = center_index - 1
        right_peak = center_index + 1
        while 1 :
            if ((right_peak >= peaks.shape[0]) or (left_peak < 0)) :
                break       #no more peaks in the left or right

            left_dist = peaks[center_index] - peaks[left_peak]
            right_dist = peaks[right_peak] - peaks[center_index] 
            if left_dist > right_dist*(1+acceptance) :
                right_peak += 1
                continue
            if right_dist > left_dist*(1+acceptance) :
                left_peak -= 1
                continue

            #peaks are in acceptance
            dist = (left_dist + right_dist) / 2
            mean_peak = (x_sum[peaks[left_peak]] + x_sum[peaks[right_peak]]) / 2
            try :
                popt_x = fittools.find_peaks_1d(x_sum, x_sum[peaks[center_index]], peaks[center_index], sigma_y0, mean_peak, peaks[center_index] + dist, sigma_y0, 0)
                break           #fit successful
            except RuntimeError :
                #1D fit failed, change to other peaks and try again
                left_peak -= 1
                right_peak += 1
                continue

        if not popt_x :
            #tough luck, find_peaks_cwt doen't return a pair of left and right peaks in accepetance
            #deprerate try
            self.guess.sigma_x0 = sigma_y0
            self.guess.sigma_y0 = sigma_y0
            self.guess.sigma_x1 = sigma_y0
            self.guess.sigma_y1 = sigma_y0
            self.find_peaks_hot_start()                                 #if fail, a RuntimeError will be raised
            return

        #We have popt_x and popt_y now
        #popt_x: a0, x0, sigma_x0, a1, x1, sigma_x1, offset
        coldguess = fittools.Guess(popt_x[3]/popt_x[0], popt_x[2], sigma_y0, popt_x[5], sigma_y0, popt_x[6]/popt[0], fx = (popt_x[4]-popt_x[1])*dXf, fy = 0)
        self.guess = coldguess
        self.find_peaks_hot_start()                                     #this shouldn't fail. Were it to fail, a RuntimeError will be raised

    def filt_move(self):
        """
        Filt the image(2d-array) xy2d and move the second peak to center. fx and fy is generated by find_peaks. xband and yband are the passband halfwidth of filter on x and y direction respectively. Return the phase spectrum.
        """
        length_x = self.xy2d.shape[1]
        length_y = self.xy2d.shape[0]
        x_filter_length = (length_x//3+1)//2*2-1      #must be odd
        y_filter_length = (length_y//3+1)//2*2-1      #must be odd
        #Filter on x direction
        b = signal.firwin(x_filter_length, cutoff=[self.fx*2-self.xband, self.fx*2+self.xband], window=('kaiser',8), pass_zero=False)
        a = np.zeros([x_filter_length])
        a[0] = 1
        xy2df = signal.filtfilt(b, a, self.xy2d)
        #Filter on y direction
        b = signal.firwin(y_filter_length, cutoff=self.fy+self.yband, window=('kaiser',8), pass_zero=True)
        a = np.zeros([y_filter_length])
        a[0] = 1
        xy2df = signal.filtfilt(b, a, xy2df, axis = 0)

        #Remove negative frequencies
        XYf2df = np.fft.fftn(xy2df)
        XYf2df[:,length_x//2:]=0
        #Shift second peak to center
        xy2df0 = np.fft.ifftn(XYf2df)
        phase = np.angle(xy2df0)
        shifter_x = np.arange(length_x)
        phase += shifter_x*(-2*np.pi*self.fx)
        shifter_y = np.arange(length_y)[:,np.newaxis]
        phase += shifter_y*(-2*np.pi*self.fy)

        #Unwrap
        phase = (phase+np.pi) % (2*np.pi) - np.pi
        self.phase = cunwrap.unwrap(phase)

    def find_symmetry_axis(self):
        self.ycenter = fittools.find_symmetry_axis(self.phase, self.symin, self.symax)

    def abel(self):
        IM = fittools.half_image(self.phase.transpose(), self.ycenter)
        self.abel_methods[self.method](IM)

    def abel_hansenlaw(self, IM):
        self.AIM = abel.hansenlaw.hansenlaw_transform(IM, direction = 'inverse').transpose()

    def abel_onion_bordas(self, IM):
        self.AIM = abel.onion_bordas.onion_bordas_transform(IM, direction = 'inverse').transpose()

    def abel_basex(self, IM):
        basex_path = os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])), 'basex')
        os.makedirs(basex_path, exist_ok=True)
        self.AIM = abel.basex.basex_transform(IM, basis_dir = basex_path, direction='inverse').transpose()

    def plot_raw(self, ax, region = None):
        """
        region: tuple of (xmin, xmax, ymin, ymax). region not equal None will add a rectagular to display selected region.
        """
        ax.set_title('Raw data')
        ax.pcolormesh(self.rawdata)
        ax.set_xlim(0, self.rawdata.shape[1])
        ax.set_ylim(0, self.rawdata.shape[0])
        if region :
            rect = patches.Rectangle((region[0], region[2]), region[1]-region[0], region[3]-region[2], linewidth=2, edgecolor='r', facecolor='none')
            ax.add_patch(rect)
        

    def plot_amplitude(self, ax, vmax=1e6, bands = None):
        """
        bands: tuple of (xband, yband). bands not equal None will add a rectagular to display passbands.
        """
        ax.set_title('Amplitude spectrum')
        xfreq = np.fft.fftshift(np.fft.fftfreq(self.XYf2d_shifted.shape[1]))
        yfreq = np.fft.fftshift(np.fft.fftfreq(self.XYf2d_shifted.shape[0]))
        ax.pcolormesh(xfreq, yfreq, self.XYf2d_shifted, vmax=vmax)
        ax.set_xlim(-0.2,0.2)
        ax.set_ylim(-0.2,0.2)
        if bands :
            rect = patches.Rectangle((self.fx-bands[0], -(np.abs(self.fy)+bands[1])), 2*bands[0], 2*(bands[1]+np.abs(self.fy)), linewidth=2, edgecolor='r', facecolor='none')
            ax.add_patch(rect)

    def plot_phase(self, ax, cax, limits = None, symmetry = None):
        """
        limits: tuple of (symin, symax). limits not equal None will add two lineouts of the limitation on symmetry axis finding.
        symmetry: y index of symmetry axis. symmetry not equal None will add the symmetry axis.
        """
        ax.set_title('Phase spectrum')
        im = ax.pcolormesh(self.phase)
        plt.colorbar(im, cax)
        if limits :
            ax.hlines(limits[0], 0, self.phase.shape[1], linewidth=2, colors='r')
            ax.hlines(limits[1], 0, self.phase.shape[1], linewidth=2, colors='r')
        if symmetry :
            ax.hlines(symmetry, 0, self.phase.shape[1], linewidth=3, colors='black')

    def plot_density(self, ax, cax, vmin=0, vmax=0.1):
        """
        Plot the result density of abel transform.
        """
        ax.set_title('Relative Refractivity')
        im = ax.pcolormesh(self.AIM, vmin=vmin, vmax=vmax)
        plt.colorbar(im, cax)
        ax.set_xlim(0, self.AIM.shape[1])
        ax.set_ylim(0, self.AIM.shape[0])

 

#def main():
#    #create a PPDD object
#    pypdd = PPDD()
#    failed_reads = []
#    failed_peaks = []
#    failed_symmetries = []
#
#    #Read Data
#    for filename in sys.argv[1:]:
#        try :
#            pypdd.readfile(filename)
#        except :
#            failed_reads.append(filename)
#            continue
#        #Fit three peaks to find the secondary peak
#        try :
#            pypdd.find_peaks()
#        except RuntimeError :
#            failed_peaks.append(filename)
#            continue
#        #Filter
#        pypdd.filt_move()
#        #Find the center of phase spectrum
#        try :
#            pypdd.find_symmetry_axis()
#        except RuntimeError :       #currently not possible because find_symmetry_axis always give a center in [ymin, ymax]
#            failed_symmetries.append(filename)
#            continue
#        #Abel transform
#        try :
#            pypdd.abel()
#        except ValueError :     #given invalid symmetry axis
#            failed_symmetries.append(filename)
#            continue
#
#        #Plot
#        plt.close("all")
#        f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(20,10))
#
#        ax1.set_title('Raw data')
#        im1 = ax1.pcolormesh(pypdd.rawdata)
#        rect1 = patches.Rectangle((pypdd.xmin, pypdd.ymin), pypdd.xmax-pypdd.xmin, pypdd.ymax-pypdd.ymin, linewidth=2, edgecolor='r', facecolor='none')
#        ax1.set_xlim(0, pypdd.rawdata.shape[1])
#        ax1.set_ylim(0, 800)
#        ax1.add_patch(rect1)
#
#        ax2.set_title('Phase spectrum')
#        im2 = ax2.pcolormesh(pypdd.phase)
#        ax2.hlines(pypdd.ycenter, 0, pypdd.phase.shape[1], linewidth=3, colors='black')
#        divider2 = make_axes_locatable(ax2)
#        cax2 = divider2.append_axes("right", size="5%", pad=0.05)
#        plt.colorbar(im2, cax2)
#
#        ax3.set_title('Amplitude spectrum')
#        XYf2d_shifted = pypdd.XYf2d_shifted
#        im3 = ax3.pcolormesh(np.fft.fftshift(np.fft.fftfreq(XYf2d_shifted.shape[1])), np.fft.fftshift(np.fft.fftfreq(XYf2d_shifted.shape[0])), XYf2d_shifted,vmax=1e6)
#        ax3.set_xlim(-0.2,0.2)
#        ax3.set_ylim(-0.2,0.2)
#        rect3 = patches.Rectangle((pypdd.fx-pypdd.xband,-np.abs(pypdd.fy)-pypdd.yband), 2*pypdd.xband, 2*(pypdd.yband+np.abs(pypdd.fy)), linewidth=2, edgecolor='r', facecolor='none')
#        ax3.add_patch(rect3)
#
#        ax4.set_title('Relative Refractivity')
#        im4 = ax4.pcolormesh(pypdd.AIM, vmax=0.1, vmin=0)
#        divider4 = make_axes_locatable(ax4)
#        cax4 = divider4.append_axes("right", size="5%", pad=0.05)
#        plt.colorbar(im4, cax4)
#
#        outputpath = os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])), 'output')
#        os.makedirs(outputpath, exist_ok=True)
#        plt.savefig(os.path.join(outputpath, os.path.basename(filename).rsplit('.', 1)[0]+'.png'), bbox_inches='tight')
#        plt.close()
#
#    if failed_reads:
#        print("Failed to read these input files:", file=sys.stderr)
#        for i in failed_reads:
#            print(i, file=sys.stderr)
#    if failed_peaks:
#        print("Failed to find the secondary peak in these input files:", file=sys.stderr)
#        for i in failed_peaks:
#            print(i, file=sys.stderr)
#    if failed_symmetries:
#        print("Failed to find the symmetry axis of phase spectrum in these input files:", file=sys.stderr)
#        for i in failed_symmetries:
#            print(i, file=sys.stderr)
#
#    return 0
#
#if __name__ == "__main__":
#    main()
