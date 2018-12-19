''' A script to run a full MCMC for a single set of patient data'''

#generic modules
import sys, os
sys.path.append("..")
import matplotlib.pyplot as mpl
%matplotlib inline
from scipy.io import loadmat
import numpy as np
import emcee
from multiprocessing import Pool

#our modules
from utils import path
from brain import Brain



param_start = [0.2,
 0.2,
 0.9670233044667013,
 14.479055812222878,
 3.9552727546221096,
 1.7950301129315296,
 0.2]
