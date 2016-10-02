from __future__ import print_function

import sys
import os
import numpy as np
import scipy
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import abel

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
        xy2d=rawdata[400:600,300:800] #TODO region
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
            ycenter, _ = abel.tools.center.find_center(phase, center='gaussian')
        except RuntimeError :
            failed_symmetries.append(filename)
            continue

        IM = abel.tools.center.center_image(phase.transpose(), center=(phase.shape[1]/2-1, ycenter), odd_size=True, crop='valid_region')
        Q0, _, _, Q3 = abel.tools.symmetry.get_image_quadrants(IM, symmetry_axis=0, symmetrize_method='average')
        Q = np.concatenate((Q0,Q3[::-1]),axis=0)[:,:] #TODO region

        #Abel transform         #TODO method selection
        try :
            AIM = abel.hansenlaw.hansenlaw_transform(Q, direction="inverse")
        except ValueError :
            failed_symmetries.append(filename)
            continue
        #AIM = abel.onion_bordas.onion_bordas_transform(Q, direction="inverse")
        #AIM = abel.basex.basex_transform(Q, nbf='auto', basis_dir='/home/clapa/.cache/abel', dr=1.0, verbose=True, direction='inverse')

        #Plot
        plt.close("all")
        fig = plt.figure(figsize=(20,8))
        plt.title('Relative Refractivity')
        plt.pcolormesh(AIM.transpose(), vmax=0.08, vmin=0)
        plt.colorbar()
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
