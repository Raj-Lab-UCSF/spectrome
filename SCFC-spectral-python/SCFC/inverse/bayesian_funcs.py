'''Collection of functions to use MCMC methods on MEG data with forward network model'''
from forward import network_transfer as nt
import numpy as np
from inverse import cost
from forward import runforward as rf


def ln_likelihood_pearson(theta, brain, fvec, FMEGdata):
    """ln_likelihood. The ln likelihood function for MCMC. Takes a brain, parameter dict,
    vector of frequency and data dict (in same order to compare to model), and returns
    negative log likelihood.

    Args:
        theta (array): Parameters of ntf model as an array (has to be converted to a dictionary).
        brain (brain class): An instance of the brain class on which to run the ntf model.
        fvec (array): List of frequencies to calculate for.
        FMEGdata (dict): Dictionary of the experimental data.

    Returns:
        float: Negative log likelihood for data given passed set of parameters.

    """
    #First need to convert parameter array to dictionary required by the run_forward func.
    if type(theta) != dict:
        parameters = {'tau_e':theta[0],
                       'tau_i':theta[1],
                       'alpha':theta[2],
                       'speed':theta[3],
                       'gei':theta[4],
                       'gii':theta[5],
                       'tauC':theta[6]
                       }
    else:
        parameters = theta

    #first, run the ntf model for this setting of parameters
    freq_model = rf.run_forward(brain, parameters, fvec)

    #now calculate the pearson r based error
    errors, list_errors = cost.pearson_cost_oneminus(FMEGdata, freq_model)

    list_errors = np.asarray(list_errors)

    #calculate the total negative lnlikelihood based on the errors of the pearson cost function
    return -0.5*(np.sum(list_errors))

#def lnprobs(theta, w, y, yerr):
