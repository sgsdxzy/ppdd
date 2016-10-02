from __future__ import print_function

import sys
import os
import numpy as np
import abel
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.axes_grid1 import make_axes_locatable

import filters
import fittools

def main():
    #create a Guess object
    guess = fittools.Guess(0.15, 0.5, 0.4, 0.8, 0.5, 3e-5, 0.06, 0.0)
    failed_peaks = []
    failed_symmetries = []

    #Read Data
    for filename in sys.argv[1:]:
        rawdata = np.loadtxt(filename, dtype=int)
        xy2d = rawdata[400:600,0:800] #TODO region
        #Fit three peaks to find the secondary peak
        try :
            fx, fy = fittools.find_peaks(xy2d, guess)
        except RuntimeError :
            failed_peaks.append(filename)
            continue

        #Filter
        phase = filters.filt_move(xy2d, fx, fy, xband=0.01, yband=0.1) #TODO passband

        #Find the center of phase spectrum
        try :
            ymin = 50
            ymax = 150
            ycenter = fittools.find_symmetry_axis(phase, ymin, ymax) #TODO ymin and ymax
        except RuntimeError :       #currently not possible because find_symmetry_axis always give a center in [ymin, ymax]
            failed_symmetries.append(filename)
            continue

        IM = abel.tools.center.center_image(phase.transpose(), center=(phase.shape[1]/2-1, ycenter), odd_size=True, crop='valid_region')
        Q0, _, _, Q3 = abel.tools.symmetry.get_image_quadrants(IM, symmetry_axis=0, symmetrize_method='average')
        Q = np.concatenate((Q0,Q3[::-1]),axis=0)[:,:] #TODO region

        #Abel transform         #TODO method selection
        try :
            AIM = abel.hansenlaw.hansenlaw_transform(Q, direction="inverse")
            #AIM = abel.onion_bordas.onion_bordas_transform(Q, direction="inverse")
            #AIM = abel.basex.basex_transform(Q, nbf='auto', basis_dir='/home/clapa/.cache/abel', dr=1.0, verbose=True, direction='inverse')
        except ValueError :     #given invalid symmetry axis
            failed_symmetries.append(filename)
            continue

        #Plot
        plt.close("all")
        f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(20,10))

        ax1.set_title('Raw data')
        im1 = ax1.pcolormesh(rawdata)
        rect1 = patches.Rectangle((0, 400), 800, 200,linewidth=2, edgecolor='r', facecolor='none')
        ax1.set_xlim(0,rawdata.shape[1])
        ax1.set_ylim(0,rawdata.shape[0])
        ax1.add_patch(rect1)

        ax2.set_title('Phase spectrum')
        im2 = ax2.pcolormesh(phase)
        ax2.hlines(ycenter, 0, phase.shape[1], linewidth=3, colors='black')
        divider2 = make_axes_locatable(ax2)
        cax2 = divider2.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(im2, cax2)

        ax3.set_title('Amplitude spectrum')
        XYf2d = np.fft.fftn(xy2d)
        XYf2d_shifted = np.abs(np.fft.fftshift(XYf2d))
        im3 = ax3.pcolormesh(np.fft.fftshift(np.fft.fftfreq(XYf2d_shifted.shape[1])), np.fft.fftshift(np.fft.fftfreq(XYf2d_shifted.shape[0])), XYf2d_shifted,vmax=1e6)
        ax3.set_xlim(-0.2,0.2)
        ax3.set_ylim(-0.2,0.2)
        rect3 = patches.Rectangle((fx-0.01,-np.abs(fy)-0.1), 0.02, 0.2+2*np.abs(fy), linewidth=2, edgecolor='r', facecolor='none')
        ax3.add_patch(rect3)

        ax4.set_title('Relative Refractivity')
        im4 = ax4.pcolormesh(AIM.transpose(), vmax=0.1, vmin=0)
        divider4 = make_axes_locatable(ax4)
        cax4 = divider4.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(im4, cax4)

        plt.savefig(os.path.basename(filename).rsplit('.', 1)[0]+'.png', bbox_inches='tight')
        plt.close()


    if failed_peaks:
        print("Failed to find the secondary peak in these input files:", file=sys.stderr)
        for i in failed_peaks:
            print(i, file=sys.stderr)
    if failed_symmetries:
        print("Failed to find the symmetry axis of phase spectrum in these input files:", file=sys.stderr)
        for i in failed_symmetries:
            print(i, file=sys.stderr)

    return 0

if __name__ == "__main__":
    main()
