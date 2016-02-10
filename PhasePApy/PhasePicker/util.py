import numpy as np

def rms(x, axis=None):
    """ Function to calculate the root mean square value of an array.
    """
    return np.sqrt(np.mean(x**2, axis=axis))   
    
def rolling_window(a, window):
    """ Efficient rolling statistics with NumPy: This is applied to Picker._statistics() to calculate statistics
        and Summary.threshold() to calcuate threshold to trigger event
        Reference from:
        http://www.rigtorp.se/2011/01/01/rolling-statistics-numpy.html
    """
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides) 