'''Collection of functions to use MCMC methods on MEG data with forward network model'''
import emcee #the affine invariant MCMC sampler from Foreman-Mackey, D., Hogg, D. W., Lang, D., Goodman, J., Feb. 2012. emcee: The MCMC Hammer.
import forward.network_transfer as fwd
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
        type: ln p for the prior probability p of param.

    """
    if len(range)%2 != 0:
        raise TypeError('range given must be of length multiple of two.')

    if range[0]<param<range[1]:
        return 0.0
    return -np.inf

def model_error(data, **kwargs):
    """This calculates the error between the data and the model.

    Args:
        data (type): Description of parameter `data`.
        params (type): Description of parameter `params`.

    Returns:
        type: Description of returned object.

    """


def ln_gauss_likelihood():
    """Returns log of likelihood of all data, where the assumption is made that errors in the
    data are independently drawn from Gaussian distributions.

    Returns:
        type: Description of returned object.

    """
    var = []
    return -0.5*np.sum((model_error()**2)*invar - np.log(invar))

def lnprobs():
    """Calculate the log of posterior probability (with no marginalization) for a single
    parameter set.

    Returns:
        type: Description of returned object.

    """
