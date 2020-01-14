import os
import sys
sys.path.append("../../")

# spectrome imports:
from spectrome.forward import runforward
from spectrome.utils import functions, path, generate
from spectrome.brain import Brain

# Other modules:
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat
from scipy.stats import pearsonr
from tqdm import tqdm_notebook

# An external custom function for re-arranging the individual connectomes to match the Desikan-Killiany atlas brain regions indexing.
def individual_connectome_order():
    """Get individual connectome's brain region ordering (specific for DK86 atlas).

    Args:

    Returns:
        perm (type): Brain region orders for all regions
        empty (type): Brain regions with no MEG
        cort (type): Brain cortical regions.

    """
    cort_lh = np.array([0, 1, 2, 3, 4, 6, 7, 8, 10, 11, 12, 13, 14,
                             15, 17, 16, 18, 19, 20, 21, 22, 23, 24, 25,
                             26, 27, 28, 29, 30, 31, 5, 32, 33, 9])
    qsubcort_lh = np.array([0, 40, 36, 39, 38, 37, 35, 34, 0])
    qsubcort_rh = qsubcort_lh + 34 + 1
    cort = np.concatenate([cort_lh, 34 + cort_lh])
    cort_rh = cort_lh + 34 + 7
    perm = np.concatenate([cort_lh, cort_rh,
                                qsubcort_lh, qsubcort_rh])
    empty = np.array([68, 77, 76, 85])
    return perm, empty, cort

# define data directory
data_dir = path.get_data_path()

# cortical areas with MEG collected + source reconstructed
rois_with_MEG = np.arange(0,68)

# define frequency bins:
fmin = 2
fmax = 45
fvec = np.linspace(fmin, fmax, 40) # 40 frequency bins between 2Hz-45Hz

# filter coefficients for smoothing
lpf = np.array([1, 2, 5, 2, 1])
lpf = lpf/np.sum(lpf)

# set a RNG seed
np.random.seed(24)

## individual connectomes, this is a Nregion x Nregion x Nsubjects array:
ind_data = loadmat(data_dir + '/individual_subjects.mat')
ind_cdk = ind_data['A_all_subjs_final']
perm, empty, cort = individual_connectome_order()

## individual optimization results:
ind_optr = loadmat(data_dir + '/SCFC_opparam_individual.mat')
# extract the spectral correlation values:
nsubs = ind_optr['output']['feval'].shape[1]
ind_corr = np.zeros(nsubs)
for i in np.arange(0,nsubs):
    if ind_optr['output']['feval'][0,i].shape[1] is 1:
        ind_corr[i] = ind_optr['output']['feval'][0,i]

ind_corr = ind_corr[ind_corr != 0] #there are 3 subjects without connectomes/results
print(ind_corr.shape)
#nsubs = ind_corr.shape[0]

# extract parameters for use later:
ind_params = np.zeros([nsubs, 7])
for i, params in enumerate(np.squeeze(ind_optr['output']['param'])):
    if ind_optr['output']['param'][0,i].shape[1] == 1:
        ind_params[i] = np.squeeze(params)

ind_params = ind_params[~np.all(ind_params == 0, axis=1)] #remove 3 subjects again
print(ind_params.shape)

## individual MEG frequency spectrum:
ind_freq = loadmat(data_dir + '/freqMEGdata.mat')
ind_psd = np.zeros([np.squeeze(ind_freq['freqMEGdata']['psd'])[0].shape[0], len(fvec), nsubs])
empty_psd = np.zeros(3)
e = 0
for i, psd in enumerate(np.squeeze(ind_freq['freqMEGdata']['psd'])):
    if psd.shape[1] != 0:
        # smooth
        for q in np.arange(0,len(psd)):
            ind_psd[q,:,i] = np.convolve(psd[q,:], lpf, 'same')
        
        # de-mean:
        #ind_psd[:,:,i] = ind_psd[:,:,i] - np.mean(ind_psd[:,:,i], axis = 0)
    else:
        empty_psd[e] = i
        e += 1

ind_psd = np.delete(ind_psd, empty_psd, axis = 2)
print(ind_psd.shape) # Nbins x Nregions x Nsubjects
nsubs = ind_params.shape[0]

num_it = 25 # 25 iterations x 36 parameters = 900 simulations total
r95_corr = np.zeros([num_it, nsubs])

for r in tqdm_notebook(np.arange(0, num_it), desc = 'Random connectomes'):
    for s in tqdm_notebook(np.arange(0, nsubs), desc = 'Current subject'):
        C_ind = ind_cdk[:,:,s] # grab current subject's individual connectome
        F_ind = ind_psd[:,:,s] # grab current subject's MEG
        
        # permute to fix ordering:
        C_ind = C_ind[perm,:][:,perm]
        C_ind[empty,:] = 0
        C_ind[:,empty] = 0
        
        # create spectrome brain:
        brain = Brain.Brain()
        brain.add_connectome(data_dir) # grabs distance matrix
        brain.connectome = C_ind # re-assign connectome to individual connectome
        # re-ordering for DK atlas and normalizing the connectomes:
        brain.reorder_connectome(brain.connectome, brain.distance_matrix)
        brain.bi_symmetric_c()
        brain.reduce_extreme_dir()
        
        # Create random connectivity matrix:
        V = 86 # number of nodes
        E = int(np.floor((86 ** 2)*(1-0.95))) # sparsity on number of edges
        u = np.mean(np.triu(brain.reducedConnectome)) #mean
        v = np.std(np.triu(brain.reducedConnectome)) #variance
        Crand = generate.random_Cij_und(V, E)
        Crand = generate.add_weights(Crand, u, v)
        brain.reducedConnectome = Crand
        
        # simulate model spectra:
        brain.ntf_params['tau_e'] = ind_params[s,0]
        brain.ntf_params['tau_i'] = ind_params[s,1]
        brain.ntf_params['alpha'] = ind_params[s,2]
        brain.ntf_params['speed'] = ind_params[s,3]
        brain.ntf_params['gei'] = ind_params[s,4]
        brain.ntf_params['gii'] = ind_params[s,5]
        brain.ntf_params['tauC'] = ind_params[s,6]
        freq_mdl, freq_resp, _, _ = runforward.run_local_coupling_forward(brain, brain.ntf_params, fvec)
        freq_mdl = freq_mdl[rois_with_MEG,:]
        # smooth out spectra
        # smooth out spectra
        freq_out = np.zeros(freq_mdl.shape)
        for p in np.arange(0,len(freq_mdl)):
            freq_out[p,:] = np.convolve(np.abs(freq_mdl[p,:]), lpf, 'same')
            
        corrs = np.zeros(len(freq_out))
        for c in np.arange(0, len(freq_out)):
            corrs[c] = pearsonr(ind_psd[c,:,s], freq_out[c,:])[0]
        
        r95_corr[r,s] = np.mean(corrs)
        
file_name = "spectral_corr_random.h5"
file_path = os.path.join(data_dir, file_name)
path.save_hdf5(file_path, r95_corr)
print(np.mean(r95_corr))