"""Module for computing basic quantities from a spectral graph model: the forward model
Makes the calculation for a single frequency only. This variation separates the global 
coupling alpha (in laplacian) and the local coupling alpha = 1. """

import numpy as np

def network_transfer_function(brain, parameters, w):
    """Network Transfer Function for spectral graph model.

    Args:
        brain (Brain): specific brain to calculate NTF
        parameters (dict): parameters for ntf. We shall keep this separate from Brain
        for now, as we want to change and update according to fitting.
        frequency (float): frequency at which to calculate NTF

    Returns:
        freqresp (numpy asarray):
        ev (numpy asarray): Eigen values
        Vv (numpy asarray): Eigen vectors
        freqresp_out (numpy asarray):  Each region's frequency response for
        the given frequency (w)
        FCmodel (numpy asarray): Functional connectivity - still in the works

    """
    C = brain.reducedConnectome
    D = brain.distance_matrix

    tau_e = parameters['tau_e']
    tau_i = parameters['tau_i']
    speed = parameters['speed']
    gei   = parameters['gei'  ]
    gii   = parameters['gii'  ]
    tauC  = parameters['tauC' ]
    global_alpha = parameters['alpha']
    local_alpha  = 1

    # Not being used: Pin = 1 and tau_syn = 0.002
    # Defining some other parameters used:
    zero_thr = 0.05
    use_smalleigs = True  # otherwise uses full eig()
    numsmalleigs = np.round(2/3*C.shape[0])  # 2/3
    a = 0.5  # fraction of signal at a node that is recurrent excitatory
    #  gei = 4 # excitatory-inhibitory synaptic conductance as ratio of E-E syn
    #  gii = 1 # inhibitory-inhibitory synaptic conductance as ratio of E-E syn
    #  tauC = 0.5*tau_e

    rowdegree = np.transpose(np.sum(C, axis=1))
    coldegree = np.sum(C, axis=0)

    qind = rowdegree + coldegree < 0.2*np.mean(rowdegree + coldegree)
    rowdegree[qind] = np.inf
    coldegree[qind] = np.inf

    nroi = C.shape[0]
    if use_smalleigs is True:
        K = numsmalleigs
        K = K.astype(int)
    else:
        K = nroi

    Tau = 0.001*D/speed
    Cc = C*np.exp(-1j*Tau*w)

    # Eigen Decomposition of Complex LaPlacian Here
    L1 = np.identity(nroi)
    L2 = np.divide(1, np.sqrt(np.multiply(rowdegree,coldegree))+np.spacing(1))
    L = L1 - global_alpha*np.matmul(np.diag(L2), Cc)

    d, v = np.linalg.eig(L) # decomposition with scipy.linalg.eig
    eig_ind = np.argsort(np.abs(d)) # sorting in ascending order and absolute value
    eig_vec = v[:, eig_ind] # re-indexing eigen vectors according to sorted index
    eig_val = d[eig_ind] # re-indexing eigen values with same sorted index

    ev = np.transpose(eig_val)
    Vv = eig_vec[:, 0:K]

    # Cortical model
    He = np.divide(1/tau_e**2, (1j*w+1/tau_e)**2)
    Hi = np.divide(gii*1/tau_i**2, (1j*w+1/tau_i)**2)

    #Hed = 1/tau_e/(1j*w + 1/tau_e*He) 
    Hed = local_alpha/tau_e/(1j*w + local_alpha/tau_e*He)
    #Hid = 1/tau_i/(1j*w + 1/tau_i*Hi)
    Hid = local_alpha/tau_i/(1j*w + local_alpha/tau_i*Hi)

    Heid = gei*He*Hi/(1+gei*He*Hi)
    Htotal = a*Hed + (1-a)/2*Hid + (1-a)/2*Heid

    q1 = 1/local_alpha*tauC*(1j*w + local_alpha/tauC*He*ev)
    #q1 = tauC*(1j*w + 1/tauC*He*ev)
    qthr = zero_thr*np.abs(q1[:]).max()
    magq1 = np.maximum(np.abs(q1), qthr)
    angq1 = np.angle(q1)
    q1 = np.multiply(magq1, np.exp(1j*angq1))
    freqresp = np.divide(Htotal, q1)

    freqresp_out = 0
    for k in range(1, K):
        freqresp_out += freqresp[k] * Vv[:, k]

    FCmodel = np.matmul(np.matmul(Vv[:, 1:K],
                        np.diag(freqresp[1:K]**2)),
                        np.transpose(Vv[:, 1:K]))

    den = np.sqrt(np.abs(freqresp_out))
    FCmodel = np.matmul(np.matmul(np.diag(1/den), FCmodel), np.diag(1/den))
    return freqresp, ev, Vv, freqresp_out, FCmodel

# def paramlist_todict(paramlist):


def network_transfer_local_alpha(brain, parameters, w):
    """Network Transfer Function for spectral graph model.

    Args:
        brain (Brain): specific brain to calculate NTF
        parameters (dict): parameters for ntf. We shall keep this separate from Brain
        for now, as we want to change and update according to fitting.
        frequency (float): frequency at which to calculate NTF

    Returns:
        freqresp (numpy asarray):
        ev (numpy asarray): Eigen values
        Vv (numpy asarray): Eigen vectors
        freqresp_out (numpy asarray):  Each region's frequency response for
        the given frequency (w)
        FCmodel (numpy asarray): Functional connectivity - still in the works

    """
    C = brain.reducedConnectome
    D = brain.distance_matrix

    tau_e = parameters['tau_e']
    tau_i = parameters['tau_i']
    speed = parameters['speed']
    gei   = parameters['gei'  ]
    gii   = parameters['gii'  ]
    tauC  = parameters['tauC' ]
    alpha = parameters['alpha']
    #local_alpha  = 1

    # Not being used: Pin = 1 and tau_syn = 0.002
    # Defining some other parameters used:
    zero_thr = 0.05
    use_smalleigs = True  # otherwise uses full eig()
    numsmalleigs = np.round(2/3*C.shape[0])  # 2/3
    a = 0.5  # fraction of signal at a node that is recurrent excitatory
    #  gei = 4 # excitatory-inhibitory synaptic conductance as ratio of E-E syn
    #  gii = 1 # inhibitory-inhibitory synaptic conductance as ratio of E-E syn
    #  tauC = 0.5*tau_e

    rowdegree = np.transpose(np.sum(C, axis=1))
    coldegree = np.sum(C, axis=0)

    qind = rowdegree + coldegree < 0.2*np.mean(rowdegree + coldegree)
    rowdegree[qind] = np.inf
    coldegree[qind] = np.inf

    nroi = C.shape[0]
    if use_smalleigs is True:
        K = numsmalleigs
        K = K.astype(int)
    else:
        K = nroi

    Tau = 0.001*D/speed
    Cc = C*np.exp(-1j*Tau*w)

    # Eigen Decomposition of Complex LaPlacian Here
    L1 = 0.8*np.identity(nroi) #0.8I in matlab
    L2 = np.divide(1, np.sqrt(np.multiply(rowdegree,coldegree))+np.spacing(1))
    L = L1 - alpha*np.matmul(np.diag(L2), Cc)

    d, v = np.linalg.eig(L) # decomposition with scipy.linalg.eig
    eig_ind = np.argsort(np.abs(d)) # sorting in ascending order and absolute value
    eig_vec = v[:, eig_ind] # re-indexing eigen vectors according to sorted index
    eig_val = d[eig_ind] # re-indexing eigen values with same sorted index

    ev = np.transpose(eig_val)
    Vv = eig_vec[:, 0:K]

    # Cortical model
    He = np.divide(1/tau_e**2, (1j*w+1/tau_e)**2)
    Hi = np.divide(gii*1/tau_i**2, (1j*w+1/tau_i)**2)

    #Hed = 1/tau_e/(1j*w + 1/tau_e*He) 
    Hed = alpha/tau_e/(1j*w + alpha/tau_e*He)
    #Hid = 1/tau_i/(1j*w + 1/tau_i*Hi)
    Hid = alpha/tau_i/(1j*w + alpha/tau_i*Hi)

    Heid = gei*He*Hi/(1+gei*He*Hi)
    Htotal = a*Hed + (1-a)/2*Hid + (1-a)/2*Heid

    q1 = 1/alpha*tauC*(1j*w + alpha/tauC*He*ev)
    #q1 = tauC*(1j*w + 1/tauC*He*ev)
    qthr = zero_thr*np.abs(q1[:]).max()
    magq1 = np.maximum(np.abs(q1), qthr)
    angq1 = np.angle(q1)
    q1 = np.multiply(magq1, np.exp(1j*angq1))
    freqresp = np.divide(Htotal, q1)

    freqresp_out = 0
    for k in range(1, K):
        freqresp_out += freqresp[k] * Vv[:, k]

    FCmodel = np.matmul(np.matmul(Vv[:, 1:K],
                        np.diag(freqresp[1:K]**2)),
                        np.transpose(Vv[:, 1:K]))

    den = np.sqrt(np.abs(freqresp_out))
    FCmodel = np.matmul(np.matmul(np.diag(1/den), FCmodel), np.diag(1/den))
    return freqresp, ev, Vv, freqresp_out, FCmodel