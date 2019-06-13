""" A script to run a full MCMC for a single set of patient data"""

# generic modules
import sys, os

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import numpy as np
import math, h5py, emcee

from brain import Brain
from utils import path as pth
from multiprocessing import Pool

from inverse.cost import pearson_cost
from inverse.bayesian_funcs import ln_likelihood_pearson, lnprobs
from inverse.priors import lnprior_uniform, lnprior_gamma, lnpriors

# from scipy.stats import pearsonr, gamma, uniform

filename = "freqMEGdata_8008-101.h5"
FMEGdownsample = pth.read_hdf5("data/downsampled_patients/freqMEGdata_8008-101.h5")
FMEGdownsample = pth.read_hdf5("data/downsampled_patients/" + str(filename))

hcp_dir = pth.get_sibling_path("SCFC/data")  # connectome information is in /data/ dir

fmin = 2  # 2Hz - 45Hz signal range
fmax = 45
fvec = np.linspace(fmin, fmax, 40)
newbrain = Brain.Brain()
newbrain.add_connectome(hcp_dir)  # Use default files in /data/
newbrain.reorder_connectome(newbrain.connectome, newbrain.distance_matrix)
newbrain.bi_symmetric_c()
newbrain.reduce_extreme_dir()

param_dists = {
    "tau_e": {
        "type": "gamma",
        "alpha": 2,
        "loc": 0.004492887509221829,
        "scale": 0.0032375996670863947,
    },
    "tau_i": {
        "type": "gamma",
        "alpha": 2.006014687419703,
        "loc": 0.004662441067103153,
        "scale": 0.0025497764353055712,
    },
    "alpha": {"type": "uniform", "lower": 0, "upper": 5},
    "speed": {"type": "uniform", "lower": 0, "upper": 25},
    "gei": {"type": "uniform", "lower": 0, "upper": 10},
    "gii": {"type": "uniform", "lower": 0, "upper": 10},
    "tauC": {
        "type": "gamma",
        "alpha": 2,
        "loc": 0.004211819821836749,
        "scale": 0.0029499360144432463,
    },
}

ndim, nwalkers = 7, 16

initial_parameters = [0.3, 0.3, 0.8, 14.4, 3.8, 1.7, 0.2]

pos = [
    np.transpose(initial_parameters) + 1e-4 * np.random.randn(ndim)
    for i in range(nwalkers)
]  # starting position in parameter space.


# Set up the backend
# Don't forget to clear it in case the file already exists
output = "chains/" + str(filename)
backend = emcee.backends.HDFBackend(output)
backend.reset(nwalkers, ndim)

## PARALLEL
n_steps = 5000
with Pool(processes=9) as pool:
    sampler = emcee.EnsembleSampler(
        nwalkers,
        ndim,
        lnprobs,
        args=(param_dists, newbrain, fvec, FMEGdownsample),
        backend=backend,
        pool=pool,
    )
    #    sampler.run_mcmc(pos, n_steps, progress=True)
    for result in sampler.sample(pos, iterations=n_steps, store=True, progress=True):
        var = 1
