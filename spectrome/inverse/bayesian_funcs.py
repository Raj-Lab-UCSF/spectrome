'''Collection of functions to use MCMC methods on MEG data with forward network model
-- uses cost (to find likelihood) and prior functions'''
import numpy as np
import math
from ..inverse import cost
from ..inverse import priors
from ..forward import runforward as rf
from ..inverse.cost import pearson_cost


# TO DO: make this aspect of the code a bit better.
#this is an unfortunately hard-coded dictionary of the parameters of prior distributions.
param_dists = {'tau_e':{'type':'gamma','alpha':2, 'loc':0.004492887509221829, 'scale':0.0032375996670863947},
               'tau_i':{'type':'gamma','alpha':2.006014687419703, 'loc':0.004662441067103153, 'scale':0.0025497764353055712},
               'alpha':{'type':'uniform','lower':0, 'upper':5},
               'speed':{'type':'uniform','lower':0, 'upper':25},
               'gei':{'type':'uniform','lower':0, 'upper':10},
               'gii':{'type':'uniform','lower':0, 'upper':10},
               'tauC':{'type':'gamma','alpha':2, 'loc':0.004211819821836749, 'scale':0.0029499360144432463}}

def ln_likelihood_pearson(theta, brain, fvec, FMEGdata):
    """Calculate negative log likelihood given parameters.

    Args:
        theta (array): Parameters of ntf model as an array.
        brain (Brain class): An instance of the brain class on which to run the ntf model.
        fvec (array): List of frequencies to calculate for.
        FMEGdata (dict): Dictionary of the experimental data.

    Returns:
        float: Negative log likelihood for data given passed set of parameters.

    """
    tau_e, tau_i, alpha, speed, gei, gii, tauC = theta
    #First need to convert parameter array to dictionary required by the run_forward func.
    if type(theta) != dict:
        parameters = {'tau_e':tau_e,
                       'tau_i':tau_i,
                       'alpha':alpha,
                       'speed':speed,
                       'gei':gei,
                       'gii':gii,
                       'tauC':tauC
                       }
    else:
        parameters = theta

    #first, run the ntf model for this setting of parameters
    freq_model = rf.run_forward(brain, parameters, fvec)

    #now calculate the pearson r based error
    errors, list_errors = pearson_cost(FMEGdata, freq_model)

    list_errors = np.asarray(list_errors)

    #calculate the total lnlikelihood based on the errors of the pearson cost function
    # change this if you think the likelihood function isn't quite right.
    # eg. might be: -0.5*(np.sum(list_errors))
    return -(np.sum(list_errors)) + 68*(1-np.log(math.e**2 - 1))

def lnprobs(theta, param_dists, brain, fvec, data):
    """Calculate unnormalised posterior probability distribution
    by summing ln(prior) and ln(p(data|parameters)).

    Args:
        theta (array): array of parameters in the standard order.
        param_dists (type): Description of parameter `param_dists`.
        brain (type): Description of parameter `brain`.
        fvec (type): Description of parameter `fvec`.
        data (type): Description of parameter `data`.

    Returns:
        type: Description of returned object.

    """
    #print(theta)
    lnprior = priors.lnpriors(theta, param_dists)
    if not np.isfinite(lnprior):
        return -np.inf
    return lnprior + ln_likelihood_pearson(theta, brain, fvec, data)
