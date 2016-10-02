import numpy as np
from scipy import signal
import unwrap

def filt_move(xy2d, fx, fy, xband=0.01, yband=0.1):
    """
    Filt the image(2d-array) xy2d and move the second peak to center. fx and fy is generated by find_peaks. xband and yband are the passband halfwidth of filter on x and y direction respectively. Return the phase spectrum.
    """
    length_x = xy2d.shape[1]
    length_y = xy2d.shape[0]
    x_filter_length = (length_x/3+1)/2*2-1      #must be odd
    y_filter_length = (length_y/3+1)/2*2-1      #must be odd
    #Filter on x direction
    b = signal.firwin(x_filter_length, cutoff=[fx*2-xband, fx*2+xband], window=('kaiser',8), pass_zero=False)
    a = np.zeros([x_filter_length])
    a[0] = 1
    xy2df = signal.filtfilt(b, a, xy2d)
    #Filter on y direction
    b = signal.firwin(y_filter_length, cutoff=fy+yband, window=('kaiser',8), pass_zero=True)
    a = np.zeros([y_filter_length])
    a[0] = 1
    xy2df = signal.filtfilt(b, a, xy2df, axis = 0)

    #Remove negative frequencies
    XYf2df = np.fft.fftn(xy2df)
    XYf2df[:,length_x/2:]=0
    #Shift second peak to center
    xy2df0 = np.fft.ifftn(XYf2df)
    phase = np.angle(xy2df0)
    shifter_x = np.arange(length_x)
    phase += shifter_x*(-2*np.pi*fx)
    shifter_y = np.arange(length_y)[:,np.newaxis]
    phase += shifter_y*(-2*np.pi*fy)

    #Unwrap
    phase = (phase+np.pi) % (2*np.pi) - np.pi
    phase = unwrap.unwrap(phase)

    return phase
