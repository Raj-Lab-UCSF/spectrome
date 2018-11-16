import numpy as np

def mag2db(y):
    """Convert magnitude response to decibels.

    Args:
        y (numpy array): Power spectrum, raw magnitude response.

    Returns:
        dby (numpy array): Power spectrum in dB

    """
    dby = 20*np.log10(y)
    return dby
