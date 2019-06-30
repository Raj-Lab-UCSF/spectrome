"""Module for computing basic quantities from a spectral graph model: the forward model
Makes the calculation for a single frequency only. This variation separates the global 
coupling alpha (in laplacian) and the local coupling alpha = 1. """

import numpy as np


def network_transfer_function(brain, parameters, w, use_smalleigs = True):
    """Network Transfer Function for spectral graph model.

    Args:
        brain (Brain): specific brain to calculate NTF
        parameters (dict): parameters for ntf. We shall keep this separate from Brain
        for now, as we want to change and update according to fitting.
        frequency (float): frequency at which to calculate NTF
        use_smalleigs (boolean): how many eigen modes to use, True = using only 2/3 (cortical), leaving out subcortical

    Returns:
        frequency_response (numpy asarray): frequency response of local oscillators
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
    gei = parameters["gei"]
    gii = parameters["gii"]
    tauC = parameters["tauC"]
    global_alpha = parameters["alpha"]
    local_alpha = 1

    # Not being used: Pin = 1 and tau_syn = 0.002
    # Defining some other parameters used:
    zero_thr = 0.05
    #use_smalleigs = True  # otherwise uses full eig()
    a = 0.5  # fraction of signal at a node that is recurrent excitatory
    #  gei = 4 # excitatory-inhibitory synaptic conductance as ratio of E-E syn
    #  gii = 1 # inhibitory-inhibitory synaptic conductance as ratio of E-E syn
    #  tauC = 0.5*tau_e

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
    L1 = np.identity(nroi)
    L2 = np.divide(1, np.sqrt(np.multiply(rowdegree, coldegree)) + np.spacing(1))
    L = L1 - global_alpha * np.matmul(np.diag(L2), Cc)

    d, v = np.linalg.eig(L)  # decomposition with scipy.linalg.eig
    eig_ind = np.argsort(np.abs(d))  # sorting in ascending order and absolute value
    eig_vec = v[:, eig_ind]  # re-indexing eigen vectors according to sorted index
    eig_val = d[eig_ind]  # re-indexing eigen values with same sorted index

    eigenvalues = np.transpose(eig_val)
    eigenvectors = eig_vec[:, 0:K]

    # Cortical model
    Fe = np.divide(1 / tau_e ** 2, (1j * w + 1 / tau_e) ** 2)
    Fi = np.divide(gii * 1 / tau_i ** 2, (1j * w + 1 / tau_i) ** 2)

    # Hed = 1/tau_e/(1j*w + 1/tau_e*He)
    Hed = local_alpha / tau_e / (1j * w + local_alpha / tau_e * Fe)
    # Hid = 1/tau_i/(1j*w + 1/tau_i*Hi)
    Hid = local_alpha / tau_i / (1j * w + local_alpha / tau_i * Fi)

    Heid = gei * Fe * Fi / (1 + gei * Fe * Fi)
    Htotal = a * Hed + (1 - a) / 2 * Hid + (1 - a) / 2 * Heid

    q1 = 1 / local_alpha * tauC * (1j * w + local_alpha / tauC * Fe * eigenvalues)
    # q1 = tauC*(1j*w + 1/tauC*He*ev)
    qthr = zero_thr * np.abs(q1[:]).max()
    magq1 = np.maximum(np.abs(q1), qthr)
    angq1 = np.angle(q1)
    q1 = np.multiply(magq1, np.exp(1j * angq1))
    frequency_response = np.divide(Htotal, q1)

    model_out = 0
    for k in range(1, K):
        model_out += frequency_response[k] * eigenvectors[:, k]

    FCmodel = np.matmul(
        np.matmul(eigenvectors[:, 1:K], np.diag(frequency_response[1:K] ** 2)), np.transpose(eigenvectors[:, 1:K])
    )

    den = np.sqrt(np.abs(model_out))
    FCmodel = np.matmul(np.matmul(np.diag(1 / den), FCmodel), np.diag(1 / den))
    return frequency_response, eigenvalues, eigenvectors, model_out, FCmodel


# def paramlist_todict(paramlist):


def network_transfer_local_alpha(brain, parameters, w, use_smalleigs = True):
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
    gei = parameters["gei"]
    gii = parameters["gii"]
    tauC = parameters["tauC"]
    alpha = parameters["alpha"]
    # local_alpha  = 1

    # Not being used: Pin = 1 and tau_syn = 0.002
    # Defining some other parameters used:
    zero_thr = 0.05
    #use_smalleigs = True  # otherwise uses full eig()
    a = 0.5  # fraction of signal at a node that is recurrent excitatory
    #  gei = 4 # excitatory-inhibitory synaptic conductance as ratio of E-E syn
    #  gii = 1 # inhibitory-inhibitory synaptic conductance as ratio of E-E syn
    #  tauC = 0.5*tau_e

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
    L1 = 0.8 * np.identity(nroi)  # 0.8I in matlab
    L2 = np.divide(1, np.sqrt(np.multiply(rowdegree, coldegree)) + np.spacing(1))
    L = L1 - alpha * np.matmul(np.diag(L2), Cc)

    d, v = np.linalg.eig(L)  # decomposition with scipy.linalg.eig
    eig_ind = np.argsort(np.abs(d))  # sorting in ascending order and absolute value
    eig_vec = v[:, eig_ind]  # re-indexing eigen vectors according to sorted index
    eig_val = d[eig_ind]  # re-indexing eigen values with same sorted index

    eigenvalues = np.transpose(eig_val)
    eigenvectors = eig_vec[:, 0:K]

    # Cortical model
    Fe = np.divide(1 / tau_e ** 2, (1j * w + 1 / tau_e) ** 2)
    Fi = np.divide(gii * 1 / tau_i ** 2, (1j * w + 1 / tau_i) ** 2)

    # Hed = 1/tau_e/(1j*w + 1/tau_e*He)
    Hed = alpha / tau_e / (1j * w + alpha / tau_e * Fe)
    # Hid = 1/tau_i/(1j*w + 1/tau_i*Hi)
    Hid = alpha / tau_i / (1j * w + alpha / tau_i * Fi)

    Heid = gei * Fe * Fi / (1 + gei * Fe * Fi)
    Htotal = a * Hed + (1 - a) / 2 * Hid + (1 - a) / 2 * Heid

    q1 = 1 / alpha * tauC * (1j * w + alpha / tauC * Fe * eigenvalues)
    # q1 = tauC*(1j*w + 1/tauC*He*ev)
    qthr = zero_thr * np.abs(q1[:]).max()
    magq1 = np.maximum(np.abs(q1), qthr)
    angq1 = np.angle(q1)
    q1 = np.multiply(magq1, np.exp(1j * angq1))
    frequency_response = np.divide(Htotal, q1)

    model_out = 0
    for k in range(1, K):
        model_out += frequency_response[k] * eigenvectors[:, k]

    FCmodel = np.matmul(
        np.matmul(eigenvectors[:, 1:K], np.diag(frequency_response[1:K] ** 2)), np.transpose(eigenvectors[:, 1:K])
    )

    den = np.sqrt(np.abs(model_out))
    FCmodel = np.matmul(np.matmul(np.diag(1 / den), FCmodel), np.diag(1 / den))
    return frequency_response, eigenvalues, eigenvectors, model_out, FCmodel

def network_transfer_HM(brain, parameters, w, use_smalleigs = True):
    """Network transfer function for spectral graph model, the local oscillator model is modified by HM.
    
    Args:
        brain (Brain): Brain class object with connectome and distance matrix
        parameters (dict): model parameters
        w (float): Frequency of interest
        use_smalleigs (boolean): how many eigen modes to use, True = using only 2/3 (cortical), leaving out subcortical
    
    Returns:
        [type]: [description]
    """
    
    # Housing keeping - defining connectomes, distance matrix, and model parameters.
    C = brain.reducedConnectome
    D = brain.distance_matrix
    tau_e = parameters["tau_e"]
    tau_i = parameters["tau_i"]
    speed = parameters["speed"]
    gei = parameters["gei"]
    gii = parameters["gii"]
    tauC = parameters["tauC"]
    global_alpha = parameters["alpha"]
    local_alpha = 1

    # Not being used: Pin = 1 and tau_syn = 0.002
    # Defining some other parameters used:
    zero_thr = 0.05
    #use_smalleigs = True  # otherwise uses full eig()
    numsmalleigs = np.round(2 / 3 * C.shape[0])  # 2/3
    a = 0.5  # fraction of signal at a node that is recurrent excitatory
    #  gei = 4 # excitatory-inhibitory synaptic conductance as ratio of E-E syn
    #  gii = 1 # inhibitory-inhibitory synaptic conductance as ratio of E-E syn
    #  tauC = 0.5*tau_e

    # define sum of degrees in rows and columns for laplacian normalization
    rowdegree = np.transpose(np.sum(C, axis=1))
    coldegree = np.sum(C, axis=0)
    qind = rowdegree + coldegree < 0.2 * np.mean(rowdegree + coldegree)
    rowdegree[qind] = np.inf
    coldegree[qind] = np.inf

    # Use all eigenmodes or 2/3 eigenmodes excluding the subcortical ones
    nroi = C.shape[0]
    if use_smalleigs is True:
        K = np.round(2 / 3 * C.shape[0])  # 2/3
        K = K.astype(int)
    else:
        K = nroi

    # Complex connectivity:
    Tau = 0.001 * D / speed # divide distance by speed, which is in meters per second, 0.001 converts D to meters
    Cc = C * np.exp(-1j * Tau * w)

    # Complex Laplacian:
    L1 = np.identity(nroi)
    L2 = np.divide(1, np.sqrt(np.multiply(rowdegree, coldegree)) + np.spacing(1))
    L = L1 - global_alpha * np.matmul(np.diag(L2), Cc)

    # eigen decomposition:
    d, v = np.linalg.eig(L)  # decomposition with scipy.linalg.eig
    eig_ind = np.argsort(np.abs(d))  # sorting in ascending order and absolute value
    eig_vec = v[:, eig_ind]  # re-indexing eigen vectors according to sorted index
    eig_val = d[eig_ind]  # re-indexing eigen values with same sorted index

    eigenvalues = np.transpose(eig_val)
    eigenvectors = eig_vec[:, 0:K] # K is either 2/3 or all eigenmodes

    # Cortical model:
    Fe = np.divide(1 / tau_e ** 2, (1j * w + 1 / tau_e) ** 2)
    Fi = np.divide(gii * 1 / tau_i ** 2, (1j * w + 1 / tau_i) ** 2)

    He = local_alpha / tau_e / (1j * w + local_alpha / tau_e * Fe)
    Hi = local_alpha / tau_i / (1j * w + local_alpha / tau_e * Fi)
    
    #denominator term for alternative model proposed by HM
    denom = 1 + (gei**2 / tau_e*tau_i) * Fe * Fi * He * Hi
    He_alt = np.divide(He,denom)
    Hi_alt = np.divide(Hi,denom)

    Hoffdiag_alt = np.divide(gei * ((-1 / tau_e) * Fe + (1/tau_i) * Fi) * He * Hi, denom)
    Htotal = He_alt + Hi_alt + Hoffdiag_alt

    q1 = 1 / local_alpha * tauC * (1j * w + local_alpha / tauC * Fe * eigenvalues)
    # q1 = tauC*(1j*w + 1/tauC*He*ev)
    qthr = zero_thr * np.abs(q1[:]).max()
    magq1 = np.maximum(np.abs(q1), qthr)
    angq1 = np.angle(q1)
    q1 = np.multiply(magq1, np.exp(1j * angq1))
    frequency_response = np.divide(Htotal, q1)

    model_out = 0
    for k in range(1, K):
        model_out += frequency_response[k] * eigenvectors[:, k]

    #FCmodel = np.matmul(
    #    np.matmul(eigenvectors[:, 1:K], np.diag(frequency_response[1:K] ** 2)), np.transpose(eigenvectors[:, 1:K])
    #)

    #den = np.sqrt(np.abs(model_out))
    #FCmodel = np.matmul(np.matmul(np.diag(1 / den), FCmodel), np.diag(1 / den))
    return frequency_response, eigenvalues, eigenvectors, model_out, Htotal