"""
Use the basinhopping algorithm to find best alpha, speed, and frequency
that produces the best DICE scores for a given canonical network
"""

# number stuff imports
import h5py
import numpy as np
import pandas as pd
from scipy.optimize import basinhopping

import sys

# spectrome imports
from spectrome.brain import Brain
from spectrome.utils import functions, path
from spectrome.forward import eigenmode, get_complex_laplacian

# hcp template connectome directory
hcp_dir = "../data"

HCP_brain = Brain.Brain()
HCP_brain.add_connectome(hcp_dir)
HCP_brain.reorder_connectome(HCP_brain.connectome, HCP_brain.distance_matrix)
HCP_brain.bi_symmetric_c()
HCP_brain.reduce_extreme_dir()

# Load Pablo's Yeo 2017 canonical network maps
com_dk = np.load(
    "../data/com_dk.npy",
    allow_pickle = True
).item()
DK_df_normalized = pd.read_csv(
    "../data/DK_dictionary_normalized.csv"
).set_index("Unnamed: 0")

# binarize:
ub, lb = 1, 0 #define binary boundaries

DKfc_binarized = pd.DataFrame(
    [], index=DK_df_normalized.index, columns=DK_df_normalized.columns
)
for name in DK_df_normalized.index:
    u = np.mean(np.nan_to_num(DK_df_normalized.loc[name].values))
    s = np.std(np.nan_to_num(DK_df_normalized.loc[name].values))
    threshold = u - s * 0.1
    DKfc_binarized.loc[name] = np.where(
        DK_df_normalized.loc[name].values > threshold, ub, lb
    )

def laplacian_dice(x, Brain, FC_networks, network_name):
    # Define frequency range of interest
    fmin = 0
    fmax = 45
    fvec = np.linspace(fmin, fmax, 50)
    f2w = np.abs(fvec - x[0]).argmin()  # 8th index = alpha ~10hz
    w = 2 * np.pi * fvec[f2w]
    # Laplacian, Brain already prep-ed with connectomes outside of function:
    Brain.add_laplacian_eigenmodes(w=w, alpha = x[1], speed = x[2], num_ev = 86)
    # Dice:
    HCP_brain.binary_eigenmodes = np.where(HCP_brain.norm_eigenmodes > 0.6, 1, 0)
    hcp_dice = eigenmode.get_dice_df(HCP_brain.binary_eigenmodes, FC_networks)
    # Compute mean Dice for chosen network:
    ntw_dice = np.round(hcp_dice[network_name].values.astype(np.double),3)
    mean_dice = np.mean(ntw_dice)
    return mean_dice

opt_res = basinhopping(laplacian_dice, x0 = (2,0.5,10), minimizer_kwargs = {"args":(HCP_brain, DKfc_binarized, str(sys.argv[1]))}, niter=250,  disp=True)
file_name = str(sys.argv[1]) + "_BH_dice.h5"
file_path = os.path.join(hcp_dir, file_name)
path.save_hdf5(file_path, opt_res)