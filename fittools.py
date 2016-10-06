import numpy as np
from scipy.optimize import curve_fit, brentq
from scipy.interpolate import interp1d

class Guess(object):
    """
    Container of guesses for fitting, used on initial fit guesses and learning.
    """
    def __init__(self, peak_ratio = 0.1, sigma_x0 = 1, sigma_y0 = 1, sigma_x1 = 1, sigma_y1 = 1, offset_ratio = 3e-5, fx = 0.06, fy = 0):
        self.peak_ratio = peak_ratio
        self.sigma_x0 = sigma_x0
        self.sigma_y0 = sigma_y0
        self.sigma_x1 = sigma_x1
        self.sigma_y1 = sigma_y1
        self.offset_ratio = offset_ratio
        self.fx = fx
        self.fy = fy

def find_nearest(array, value):
    """
    Find the index of nearest element in array to value.
    """
    idx = (np.abs(array-value)).argmin()
    return idx

def gaussian(x, a, mu, sigma, c):
    """
    Gaussian function

    :math:`f(x)=a e^{-(x - \mu)^2 / (2 \\sigma^2)} + c`

    ref: https://en.wikipedia.org/wiki/Gaussian_function

    Parameters
    ----------
    x : 1D np.array
        coordinate

    a : float
        the height of the curve's peak

    mu : float
        the position of the center of the peak

    sigma : float
        the standard deviation, sometimes called the Gaussian RMS width

    c : float
        non-zero background

    Returns
    -------
    out : 1D np.array
        the Gaussian profile
    """
    return a * np.exp(-((x - mu) ** 2) / 2 / sigma ** 2) + c

def guss_gaussian(x):
    """
    Find a set of better starting parameters for Gaussian function fitting

    Parameters
    ----------
    x : 1D np.array
        1D profile of your data

    Returns
    -------
    out : tuple of float
        estimated value of (a, mu, sigma, c)
    """
    c_guess = (x[0] + x[-1]) / 2
    a_guess = x.max() - c_guess
    mu_guess = x.argmax()
    x_inter = interp1d(np.arange(len(x)), x)

    def _(i):
        return x_inter(i) - a_guess / 2 - c_guess

    try:
        sigma_l_guess = brentq(_, 0, mu_guess)
    except:
        sigma_l_guess = len(x) / 4
    try:
        sigma_r_guess = brentq(_, mu_guess, len(x) - 1)
    except:
        sigma_r_guess = 3 * len(x) / 4
    return a_guess, mu_guess, (sigma_r_guess -
                               sigma_l_guess) / 2.35482, c_guess

def fit_gaussian(x, xmin, xmax):
    """
    Fit a Gaussian function to x and return its parameters, with mu in [xmin, xmax]

    Parameters
    ----------
    x : 1D np.array
        1D profile of your data

    Returns
    -------
    out : tuple of float
        (a, mu, sigma, c)
    """
    p, q = curve_fit(gaussian, np.arange(x.size), x, p0=guss_gaussian(x), bounds=([-np.inf, xmin, -np.inf, -np.inf], [np.inf, xmax, np.inf, np.inf]))
    return p

def find_center_by_gaussian_fit(IM, ymin, ymax):
    """
    Find image center by fitting the summation along x and y axis of the data to two 1D Gaussian function
    """
    y = np.sum(IM, axis=1)
    return fit_gaussian(y, ymin, ymax)[1]

def find_center_by_convolution(IM, ymin, ymax):
    """ Center the image by convolution of two projections along each axis.
        code from the ``linbasex`` juptyer notebook
    Parameter
    -------
    IM: numpy 2D array
        image data
    Returns
    -------
        y-center
    """
    # projection along axis=0 of image (rows)
    QL_raw0 = IM.sum(axis=1)

    # autocorrelate projections
    conv_0 = np.convolve(QL_raw0, QL_raw0, mode='full')

    #Take the first max, should there be several equal maxima.
    # 10May16 - axes swapped - check this
    return np.argmax(conv_0[ymin*2:ymax*2])/2 + ymin

def find_symmetry_axis(phase, ymin, ymax):
    """
    Find symmetry axis of phase spectrum in range [ymin, ymax]. It will try different methods in the following order:
    find_center_by_gaussian_fit
    find_center_by_convolution
    If none of the methods could find a valid symmetry axis, a RuntimeError will be raised.

    Return the y index of the symmetry axis.
    """
    try :
        center = find_center_by_gaussian_fit(phase, ymin, ymax)
        return center
    except (RuntimeError, ValueError) :
        #find_center_by_gaussian_fit failed, just pass to use next method
        pass
    
    #find_center_by_convolution always succeeds
    center = find_center_by_convolution(phase, ymin, ymax)
    return center

def three_peaks_1d(x, a0, x0, sigma_x0, a1, x1, sigma_x1, offset):
    """
    The 1D fitting function for fitting three peaks in projection on x axis.
    """
    peak0 = gaussian(x, a0, x0, sigma_x0, 0)
    peak1 = gaussian(x, a1, x1, sigma_x1, 0)
    peakm1 = gaussian(x, a1, 2*x0-x1, sigma_x1, 0)
    return peak0 + peak1 + peakm1 + offset

def find_peaks_1d(x, a0, x0, sigma_x0, a1, x1, sigma_x1, offset):
    length_x = x.shape[0]
    popt,_ = curve_fit(three_peaks_1d, np.arange(length_x), x, p0 = (a0, x0, sigma_x0, a1, x1, sigma_x1, offset), 
            bounds = ([-np.inf, 0, 0, -np.inf, length_x//2, 0, -np.inf], [np.inf, length_x, np.inf, np.inf, length_x, max(0.01*length_x, 5), np.inf]))
            #needs to limit sigma to avoid unsense results
    return popt

def three_peaks(xy_tuple, a0, x0, y0, sigma_x0, sigma_y0, a1, x1, y1, sigma_x1, sigma_y1, offset):
    """
    The fitting function of three peaks.
    """
    (x, y) = xy_tuple
    peak0 = a0*np.exp((-(x-x0)**2)/(2*sigma_x0**2) + (-(y-y0)**2)/(2*sigma_y0**2))
    peak1 = a1*np.exp((-(x-x1)**2)/(2*sigma_x1**2) + (-(y-y1)**2)/(2*sigma_y1**2))
    peakm1 = a1*np.exp((-(x+x1-2*x0)**2)/(2*sigma_x1**2) + (-(y+y1-2*y0)**2)/(2*sigma_y1**2))
    g = peak0 + peak1 + peakm1 + offset
    return g.ravel()
            
def find_peaks(XYf2d_shifted, guess):
    """
    Fit the three peaks in the shifted 2d amplitude spectrum XYf2d_shifted. 
    Return the phase shift of the secondary peak in x and y direction.
    """    
    length_x = XYf2d_shifted.shape[1]
    length_y = XYf2d_shifted.shape[0]
    dXf = 1/length_x
    dYf = 1/length_y
                 
    a0 = np.max(XYf2d_shifted)                                    #compose initial fit condition from guess
    x0 = length_x//2
    y0 = length_y//2
    a1 = guess.peak_ratio*a0
    x1 = x0 + guess.fx/dXf
    y1 = y0 + guess.fy/dYf
    offset = guess.offset_ratio*a0 
    initial_guess = (a0, x0, y0, guess.sigma_x0, guess.sigma_y0, a1, x1, y1, guess.sigma_x1, guess.sigma_y1, offset)     
    x, y = np.meshgrid(np.arange(length_x), np.arange(length_y))       
    popt,_ = curve_fit(three_peaks, (x, y), XYf2d_shifted.ravel(), p0=initial_guess, 
            bounds = ([-np.inf, 0, 0, 0, 0, -np.inf, length_x//2, 0, 0, 0, -np.inf], 
                [np.inf, length_x, length_y, np.inf, np.inf, np.inf, length_x, length_y, max(0.01*length_x, 5), max(0.01*length_y, 5), np.inf]))   
            #needs to limit sigma to avoid unsense results
    
    fx = (popt[6]-popt[1])*dXf
    fy = (popt[7]-popt[2])*dYf

    newguess = Guess()
    newguess.peak_ratio = popt[5]/popt[0]                            #update guess
    newguess.sigma_x0 = popt[3]
    newguess.sigma_y0 = popt[4]
    newguess.sigma_x1 = popt[8]
    newguess.sigma_y1 = popt[9]
    newguess.offset_ratio = popt[10]/popt[0]
    newguess.fx = fx
    newguess.fy = fy
    
    #xband1 = 0.09#100*popt[3]*dXf/0.5                            #not used
    #xband2 = 0.16#(popt[6]-popt[1]+30*popt[8])*dXf/0.5
    #yband = 0.12#80*popt[9]*dYf/0.5
    
    return fx, fy, newguess

def half_image(IM, xcenter):
    """
    Generate half of image IM by the image center in the x direction. This function is used to prepare for abel transfrom.
    """
    xcenter = int(np.rint(xcenter))
    new_width = min(IM.shape[1] - xcenter - 1, xcenter)
    left = IM[:, xcenter-new_width:xcenter+1][:, ::-1]
    right = IM[:, xcenter:xcenter+new_width+1]
    return (left + right) / 2
