''' A script to run a full MCMC for a single set of patient data'''

#generic modules
import sys, os
import numpy as np
import math, h5py, emcee

from brain import Brain
from utils import path as pth
from multiprocessing import Pool

from inverse.cost import pearson_cost
from inverse.bayesian_funcs import ln_likelihood_pearson, lnprobs
from inverse.priors import lnprior_uniform, lnprior_gamma, lnpriors

# from scipy.stats import pearsonr, gamma, uniform

FMEGdownsample = pth.read_hdf5('data/freqMEGdata_8008_101.h5')

hcp_dir = pth.get_sibling_path('SCFC/data') # connectome information is in /data/ dir

fmin = 2 # 2Hz - 45Hz signal range
fmax = 45
fvec = np.linspace(fmin,fmax,40)
newbrain = Brain.Brain()
newbrain.add_connectome(hcp_dir) # Use default files in /data/
newbrain.reorder_connectome(newbrain.connectome, newbrain.distance_matrix)
newbrain.bi_symmetric_c()
newbrain.reduce_extreme_dir()

param_dists = {'tau_e':{'type':'gamma','alpha':2, 'loc':0.004492887509221829, 'scale':0.0032375996670863947},
              'tau_i':{'type':'gamma','alpha':2.006014687419703, 'loc':0.004662441067103153, 'scale':0.0025497764353055712},
              'alpha':{'type':'uniform','lower':0, 'upper':5},
              'speed':{'type':'uniform','lower':0, 'upper':25},
              'gei':{'type':'uniform','lower':0, 'upper':10},
              'gii':{'type':'uniform','lower':0, 'upper':10},
              'tauC':{'type':'gamma','alpha':2, 'loc':0.004211819821836749, 'scale':0.0029499360144432463}}

ndim, nwalkers = 7, 16

initial_parameters = [0.3,0.3,0.8,14.4,3.8,1.7,0.2]

pos = [np.transpose(initial_parameters) + 1e-4*np.random.randn(ndim) for i in range(nwalkers)] #starting position in parameter space.


# Set up the backend
# Don't forget to clear it in case the file already exists
filename = "saved_sampler.h5"
backend = emcee.backends.HDFBackend(filename)
backend.reset(nwalkers, ndim)

## PARALLEL
max_n = 1
burnout_time = 100
index = 0
autocorr = np.empty(max_n)
old_tau = np.inf

from emcee.autocorr import integrated_time

with Pool() as pool:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprobs,
                                    args=(param_dists, newbrain, fvec, FMEGdownsample),
                                    backend=backend,
                                    pool=pool)

    for sample in sampler.sample(pos, iterations=max_n, progress=True):
        if sampler.iteration > burnout_time and sampler.iteration % 10:
            chain = sampler.chain

            tau = np.mean([integrated_time(walker, c=1, tol=1, quiet=True) for walker in chain], axis=0)
            # tau = sampler.get_autocorr_time(tol=0)
            autocorr[index] = np.mean(tau)
            index += 1

            converged = np.all(tau * 10 < sampler.iteration)
            converged &= np.all(np.abs(old_tau - tau) / tau < 0.01)
            if converged:
                break
            old_tau = tau


sampler.chain.shape
# import matplotlib.pyplot as plt
# autocorr.shape
#
# n = np.arange(1, index+1)
# y = autocorr[:index]
# plt.plot(n, n / 10.0, "--k")
# plt.plot(n, y)
# plt.xlim(0, n.max())
# plt.ylim(0, y.max() + 0.1*(y.max() - y.min()))
# plt.xlabel("number of steps")
# plt.ylabel(r"mean $\hat{\tau}$");
#
# autocorr
# index
#
# for i in range(sampler.chain.shape[0]):
#     mpl.plot(np.arange(500), sampler.chain[i,-500:,0], linewidth = 0.1, color = 'black')
#
# sampler.chain[0,-50:,0]
#
# ax3 = mpl.gca()
# ax3.set_xlabel('step')
# ax3.set_ylabel('tau_e')
#
# sampler.chain[0,:50,0]
