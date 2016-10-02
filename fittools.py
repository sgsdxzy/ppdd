import numpy as np
import scipy

class Guess(object):
    """
    Container of guesses for fitting, used on initial fit guesses and learning.
    """
    def __init__(self, peak_ratio, sigma_x0, sigma_y0, sigma_x1, sigma_y1, offset_ratio, fx, fy):
        self.peak_ratio = peak_ratio
        self.sigma_x0 = sigma_x0
        self.sigma_y0 = sigma_y0
        self.sigma_x1 = sigma_x1
        self.sigma_y1 = sigma_y1
        self.offset_ratio = offset_ratio
        self.fx = fx
        self.fy = fy

def three_peaks((x, y), a0, x0, y0, sigma_x0, sigma_y0, a1, x1, y1, sigma_x1, sigma_y1, offset):
    """
    The fitting function of three peaks.
    """
    x0 = float(x0)
    y0 = float(y0)
    x1 = float(x1)
    y1 = float(y1)
    peak0 = a0*np.exp((-(x-x0)**2)/(2*sigma_x0**2) + (-(y-y0)**2)/(2*sigma_y0**2))
    peak1 = a1*np.exp((-(x-x1)**2)/(2*sigma_x1**2) + (-(y-y1)**2)/(2*sigma_y1**2))
    peakm1 = a1*np.exp((-(x+x1-2*x0)**2)/(2*sigma_x1**2) + (-(y+y1-2*y0)**2)/(2*sigma_y1**2))
    g = peak0 + peak1 + peakm1 + offset
    return g.ravel()
            
def find_peaks(xy2d, guess):
    """
    Fit the three peaks in the frequency spectrum of image(2d-array) xy2d. 
    Return the phase shift of the secondary peak in x and y direction.
    """    
    length_x = xy2d.shape[1]
    length_y = xy2d.shape[0]
    dXf = float(1.0/length_x)
    dYf = float(1.0/length_y)
    XYf2d = np.fft.fftn(xy2d)
    XYf2d_shifted = np.abs(np.fft.fftshift(XYf2d))                #shift frequency of (0,0) to the center
                 
    a0 = np.max(XYf2d_shifted)                                               #compose initial fit condition from guess
    x0 = length_x/2
    y0 = length_y/2
    a1 = guess.peak_ratio*a0
    x1 = x0 + guess.fx/dXf
    y1 = y0 + guess.fy/dYf
    offset = guess.offset_ratio*a0 
    initial_guess = (a0, x0, y0, guess.sigma_x0, guess.sigma_y0, a1, x1, y1, guess.sigma_x1, guess.sigma_y1, offset)     
    x, y = np.meshgrid(np.arange(length_x), np.arange(length_y))       
    popt,_ = scipy.optimize.curve_fit(three_peaks, (x, y), XYf2d_shifted.ravel(), p0=initial_guess)
    
    fx = (popt[6]-popt[1])*dXf
    fy = (popt[7]-popt[2])*dYf
    guess.peak_ratio = popt[5]/popt[0]                            #update guess
    guess.sigma_x0 = popt[3]
    guess.sigma_y0 = popt[4]
    guess.sigma_x1 = popt[8]
    guess.sigma_y1 = popt[9]
    guess.offset_ratio = popt[10]/popt[0]
    guess.fx = fx
    guess.fy = fy
    
    #xband1 = 0.09#100*popt[3]*dXf/0.5                            #not used
    #xband2 = 0.16#(popt[6]-popt[1]+30*popt[8])*dXf/0.5
    #yband = 0.12#80*popt[9]*dYf/0.5
    
    return fx, fy
