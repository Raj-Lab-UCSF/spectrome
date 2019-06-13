"""Looking at changes in Eigenmode-FC similarity as a function of alpha and frequency
*** This code needs debugging, currently producing same figures every iteration *** 
"""
# number stuff imports
import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

sys.path.append("../../")

# spectrome imports
from spectrome.brain import Brain
from spectrome.utils import functions, path
from spectrome.forward import eigenmode, get_complex_laplacian

# from brain import Brain
# from utils import functions, path
# from forward import eigenmode, get_complex_laplacian


# house keeping things: directory for hcp dataset
hcp_dir = "/home/axiezai/lab/brain-scfc/spectrome/spectrome/data"

# Define frequency range of interest
fmin = 2  # 2Hz - 45Hz signal range, filter for this with hbp
fmax = 45
fvec = np.linspace(fmin, fmax, 44)

# Define range of alpha/frequency parameter
alpha = np.linspace(0.5, 15, 5)
freq = np.linspace(2, 45, 5)


# load the canonical networks (Yeo 2017)
com_dk = np.load(
    "/home/axiezai/lab/brain-scfc/spectrome/spectrome/data/com_dk.npy"
).item()
DK_df_normalized = pd.read_csv(
    "/home/axiezai/lab/brain-scfc/spectrome/spectrome/data/DK_dictionary_normalized.csv"
).set_index("Unnamed: 0")
coords = np.array([com_dk[region] for region in DK_df_normalized.columns])

# we need to threshold so that the sparsity of all canonical networks are similar to each other
upperb, lowerb = 1, 0
Dkfc_binarized = pd.DataFrame(
    [], index=DK_df_normalized.index, columns=DK_df_normalized.columns
)

for names in DK_df_normalized.index:
    u = np.mean(np.nan_to_num(DK_df_normalized.loc[names].values))
    s = np.std(np.nan_to_num(DK_df_normalized.loc[names].values))
    threshold = u - s * 0.1
    Dkfc_binarized.loc[names] = np.where(
        DK_df_normalized.loc[names].values > threshold, upperb, lowerb
    )
    # print(
    #    "Number of non-zero elements in {} network is {}".format(
    #        names, np.count_nonzero(DKfc_binarized.loc[names].values)
    #    )
    # )


# Creating HCP brain object for spectrome and eigenmodes
HCP_brain = Brain.Brain()
HCP_brain.add_connectome(hcp_dir)  # loads HCP connectome
HCP_brain.reorder_connectome(HCP_brain.connectome, HCP_brain.distance_matrix)
HCP_brain.bi_symmetric_c()
HCP_brain.reduce_extreme_dir()


# First scan through alpha:
omega = 2 * np.pi * fvec[np.abs(fvec - 10).argmin()]  # current frequency of interest
eigthresh = 0.6  # see jupyter notebook for how this was obtained
networks_alpha = pd.DataFrame(
    np.zeros([len(alpha), Dkfc_binarized.shape[0]]), index=alpha, columns=Dkfc_binarized.index
)
for param_alpha in alpha:
    HCP_brain.add_laplacian_eigenmodes(
        # decompose laplacian and obtain eigenmodes
        w=omega,
        alpha=param_alpha,
        speed=10,
        num_ev=86,
    )

    # binarize eigen modes based on threshold
    HCP_brain.binary_eigenmodes = np.where(
        HCP_brain.norm_eigenmodes > eigthresh, upperb, lowerb
    )

    # compute dice dissimilarity values for all eigenmodes
    hcp_dice = eigenmode.get_dice_df(HCP_brain.binary_eigenmodes, Dkfc_binarized)

    # compute number of eigenmodes most similar to each network
    for eignum in hcp_dice.index:
        iter_em = pd.to_numeric(hcp_dice.loc[eignum, :])
        best_network = iter_em.idxmin(axis=0)  # find lowest dice dissimilarity
        networks_alpha.at[param_alpha, best_network] += 1

    # networks_count.rename(index={i:'{}'.format(alpha[i])})
    # print(networks_count)

networks_alpha.rename(columns={"Unnamed: 0": "Canonical Networks"}, inplace=True)
# plot results:
ax_alpha = networks_alpha.plot.bar(rot=0, title="alpha")


# scan through frequencies
networks_freq = pd.DataFrame(
    np.zeros([len(alpha), hcp_dice.shape[1]]), index=freq, columns=hcp_dice.columns
)

for f in freq:
    omega = 2 * np.pi * fvec[np.abs(fvec - f).argmin()]
    HCP_brain.add_laplacian_eigenmodes(
        # decompose laplacian and obtain eigenmodes
        w=omega,
        alpha=1,
        speed=10,
        num_ev=86,
    )

    # binarize eigen modes based on threshold
    HCP_brain.binary_eigenmodes = np.where(
        HCP_brain.norm_eigenmodes > eigthresh, upperb, lowerb
    )

    # compute dice dissimilarity values for all eigenmodes
    hcp_dice = eigenmode.get_dice_df(HCP_brain.binary_eigenmodes, Dkfc_binarized)
    # compute number of eigenmodes most similar to each network
    for eignum in hcp_dice.index:
        iter_em = pd.to_numeric(hcp_dice.loc[eignum, :])
        best_network = iter_em.idxmin(axis=0)  # find lowest dice dissimilarity
        networks_freq.at[f, best_network] += 1

# plot results:
ax_freq = networks_freq.plot.bar(rot=0, title="frequency")
