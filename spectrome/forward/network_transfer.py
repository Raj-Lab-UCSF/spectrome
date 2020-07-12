"""Module for computing basic quantities from a spectral graph model: the forward model
Makes the calculation for a single frequency only. This variation separates the global 
coupling alpha (in laplacian) and the local coupling alpha = 1. """

import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft

def network_transfer_local_alpha(brain, parameters, w, use_smalleigs=True):
    """Network Transfer Function for spectral graph model.

    Args:
        brain (Brain): specific brain to calculate NTF
        parameters (dict): parameters for ntf. We shall keep this separate from Brain
        for now, as we want to change and update according to fitting.
        frequency (float): frequency at which to calculate NTF
        use_smalleigs (boolean): how many eigen modes to use, True = using only 2/3 (cortical), leaving out subcortical

    Returns:
        frequency_response (numpy asarray):
        ev (numpy asarray): Eigen values
        Vv (numpy asarray): Eigen vectors
        model_out (numpy asarray):  Each region's frequency response for
        the given frequency (w)
        FCmodel (numpy asarray): Functional connectivity - still in the works

    """
    C = brain.reducedConnectome
    D = brain.distance_matrix

    tau_e = parameters["tau_e"]
    tau_i = parameters["tau_i"]
    speed = parameters["speed"]
    gei = parameters[
        "gei"
    ]  # excitatory-inhibitory synaptic conductance as ratio of E-E syn
    gii = parameters[
        "gii"
    ]  # inhibitory-inhibitory synaptic conductance as ratio of E-E syn
    tauC = parameters["tauC"]  #  tauC = 0.5*tau_e
    alpha = parameters["alpha"]
    local_alpha  = 1
#     alpha = 3
    ## Values from Ashish code
#     tau_e = 0.012
#     tau_i = 0.008
#     alpha = 0.6
#     tauC = 0.5*tau_e
#     speed = 10
#     zero_thr = 0.01
#     alpha = 1.6
#     speed = 1
    # Not being used: Pin = 1 and tau_syn = 0.002
    # Defining some other parameters used:
    zero_thr = 0.05
    a = 0.5  # fraction of signal at a node that is recurrent excitatory

    # define sum of degrees for rows and columns for laplacian normalization
    rowdegree = np.transpose(np.sum(C, axis=1))
    coldegree = np.sum(C, axis=0)
    qind = rowdegree + coldegree < 0.2 * np.mean(rowdegree + coldegree)
    rowdegree[qind] = np.inf
    coldegree[qind] = np.inf

    nroi = C.shape[0]
    if use_smalleigs is True:
        K = np.round(2 / 3 * C.shape[0])  # 2/3
        K = K.astype(int)
    else:
        K = nroi

    Tau = 0.001 * D / speed
    Cc = C * np.exp(-1j * Tau * w)

    # Eigen Decomposition of Complex Laplacian Here
    #L1 = 0.8 * np.identity(nroi)  # 0.8I in matlab
    L1 = np.identity(nroi)
    L2 = np.divide(1, np.sqrt(np.multiply(rowdegree, coldegree)) + np.spacing(1))
    L = L1 - alpha * np.matmul(np.diag(L2), Cc)

    d, v = np.linalg.eig(L)  # decomposition with scipy.linalg.eig
    eig_ind = np.argsort(np.abs(d))  # sorting in ascending order and absolute value
    eig_vec = v[:, eig_ind]  # re-indexing eigen vectors according to sorted index
    eig_val = d[eig_ind]  # re-indexing eigen values with same sorted index

    eigenvalues = np.transpose(eig_val)
    eigenvectors = eig_vec[:, 0:K]

#     # Cortical model
    Fe = np.divide(1 / tau_e ** 2, (1j * w + 1 / tau_e) ** 2)
    Fi = np.divide(1 / tau_i ** 2, (1j * w + 1 / tau_i) ** 2)


    Hed = 1 / (1j * w + 1 / tau_e * Fe)
    Hid = 1 / (1j * w + gii / tau_i * Fi)

#     Fei = Fe * Fi / (1 + gei * Fe * Fi)
#     Heid = 1 / (1j * w + gei / tau_e * Fei)
    Heid = Hed * Hid / (1 + gei * Hed * Hid)
#     Htotal = a * Hed + (1 - a) / 2 * Hid + (1 - a) / 2 * Heid
    Htotal = Hed + Hid + Heid
    

    q1 = (1j * w + 1 / tauC * Fe * eigenvalues)
    # q1 = tauC*(1j*w + 1/tauC*He*ev)
    qthr = zero_thr * np.abs(q1[:]).max()
    magq1 = np.maximum(np.abs(q1), qthr)
    angq1 = np.angle(q1)
    q1 = np.multiply(magq1, np.exp(1j * angq1))
    frequency_response = np.divide(Htotal, q1)
    
    model_out = 0
    ## Changed model_out to take into account P(w) which is a vector of ones. Check if conjugate is correct!!!
    for k in range(1, K):
        model_out += (frequency_response[k]) * eigenvectors[:, k] * sum(np.conjugate(eigenvectors[:, k])) 

    den = np.sqrt(np.abs(model_out))
#     FCmodel = np.matmul(np.matmul(np.diag(1 / den), FCmodel), np.diag(1 / den))
    return model_out, frequency_response, eigenvalues, eigenvectors

