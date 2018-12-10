'''Collection of functions to use MCMC methods on MEG data with forward network model'''
from forward import network_transfer as nt
import numpy as np
from inverse import cost
from forward import runforward as rf


def ln_likelihood(theta, brain, fvec, FMEGdata):
    """ln_likelihood. The ln likelihood function for MCMC. Takes a brain, parameter dict,
    vector of frequency and data dict (in same order to compare to model), and returns
    negative log likelihood.

    Args:
        theta (dict): Parameters of ntf model as a dictionary (need to change to arrange for MCMC).
        brain (brain class): An instance of the brain class on which to run the ntf model.
        fvec (array): List of frequencies to calculate for.
        FMEGdata (dict): Dictionary of the experimental data.

    Returns:
        float: Negative log likelihood for data given passed set of parameters.

    """
    #First need to convert parameter array to dictionary required by the run_forward func.

    #first, run the ntf model for this setting of parameters
    freq_model = rf.run_forward(brain, theta, fvec)

    #now calculate the pearson r based error
    errors, list_errors = cost.pearson_cost(FMEGdata, freq_model)

    list_errors = np.asarray(list_errors)

    #calculate the total negative lnlikelihood based on the errors of the pearson cost function
    return -0.5*(np.sum(list_errors**2))

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
