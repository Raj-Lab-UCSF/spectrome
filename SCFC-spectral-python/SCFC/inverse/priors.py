from scipy import stats
import math
import numpy as np

#this is an unfortunately hard-coded dictionary of the parameters of prior distributions.
param_dists = {'tau_e':{'type':'gamma','alpha':2, 'loc':0.004492887509221829, 'scale':0.0032375996670863947},
    'tau_i':{'type':'gamma','alpha':2.006014687419703, 'loc':0.004662441067103153, 'scale':0.0025497764353055712},
              'alpha':{'type':'uniform','lower':0, 'upper':5},
              'speed':{'type':'uniform','lower':0, 'upper':25},
              'gei':{'type':'uniform','lower':0, 'upper':10},
              'gii':{'type':'uniform','lower':0, 'upper':10},
              'tauC':{'type':'gamma','alpha':2, 'loc':0.004211819821836749, 'scale':0.0029499360144432463}}


def lnprior_gamma(param, key, param_dists):
    """lnprior_gamma. Generates a -ln(prior) probability for a gamma distributed parameter.

    Args:
        param (float): value of parameter at which the -ln(prob) is required.
        key (string): name of the parameter.
        param_dists (dict): dictionary of the parameter distribution required parameters. .

    Returns:
        float: -ln(p) for the parameter, assuming the parameter prior distribution given.

    """
    tau = param
    x = (tau - param_dists[key]['loc'])/param_dists[key]['scale']
    p = stats.gamma.pdf(x, param_dists[key]['alpha'])
    p = p.astype(float)
    if p == 0:
        return -np.inf #use the same as for outside uniform distribution bounds
    return -math.log(p)

def lnprior_uniform(param, key, param_dists):
    """lnprior_uniform. Generates a -ln(prior) probability for a uniformly distributed parameter.

    Args:
        param (float): value of parameter at which the -ln(prob) is required.
        key (string): name of the parameter.
        param_dists (dict): dictionary of the parameter distribution required parameters.

    Returns:
        float: -ln(p) for the parameter, assuming the parameter prior distribution given.

    """
    if param_dists[key]['lower'] < param < param_dists[key]['upper']:
        return 0.0
    return -np.inf

def lnpriors(params, param_dists):
    i = 0
    neg_lnprior = 0
    priorlist = []
    for key in param_dists.keys():
        if param_dists[key]['type'] == 'gamma':
            #print(params[i], type(params[i]))
            neg_lnp = lnprior_gamma(params[i], key, param_dists)
        elif param_dists[key]['type'] == 'uniform':
            #print(params[i], type(params[i]))
            neg_lnp = lnprior_uniform(params[i], key, param_dists)
        priorlist.append(neg_lnp)
        neg_lnprior += neg_lnp
        i += 1
    return neg_lnprior, priorlist
