'''Collection of functions to use MCMC methods on MEG data with forward network model'''
from forward import network_transfer as nt
import numpy as np

def ln_uniform_prior(param, range):
    """ln_uniform_prior: returns log of uniform prior for a parameter,
    up to a constant -- so returns zero if parameter 'in range' and -inf
    otherwise.

    Args:
        param (type): the parameter.
        range (type): range over which a non-zero prior probability is defined
        Should be a 2 element list, tuple, or array.

    Returns:
        type: ln p for the prior probability of the parameter on question.

    """
    if len(range)%2 != 0:
        raise TypeError('range given must be of length multiple of two.')

    if range[0]<param<range[1]:
        return 0.0 #ln prior up to a constant that doesn't matter
    return -np.inf

def model_error(data, **kwargs):
    #calculate error between model and data


def ln_gauss_likelihood():
    #assume our errors are drawn independently from Gaussian and calculate likelihood.
    var = []
    return -0.5*np.sum((model_error()**2)*invar - np.log(invar))

def lnprobs(theta, w, y, yerr):
    """Calculate the log of posterior probability (with no marginalization) for a single
    parameter set..

    Args:
        theta (type): the parameters of the model that will be called.
        w (type): frequency array associated with data points.
        y (type): real data points.
        yerr (type): variance of the distribution from which errors are drawn, as array.

    Returns:
        type: for single theta set, the log(likelihood)xlog(prior).

    """

    return
