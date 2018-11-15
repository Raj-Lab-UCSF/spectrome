"""Module for computing basic quantities from a spectral graph model."""

import numpy as np
from scipy.io import loadmat
import os
from scipy.stats import pearsonr

def network_transfer_function(brain=brain, frequency=w):
    """Network Transfer Function for spectral graph model.

    Args:
        brain (Brain): specific brain to calculate NTF
        frequency (float): frequency at which to calculate NTF

    Returns:
        freqresp (numpy asarray):
        ev (numpy asarray): Eigen values
        Vv (numpy asarray): Eigen vectors
        freqresp_out (numpy asarray):  Each region's frequency response for
        the given frequency (w)
        FCmodel (numpy asarray): Functional connectivity - still in the works

    """
    ntf_params = brain.ntf_parameters()
    tau_e = ntf_params['tau_e']
    tau_i = ntf_params['tau_i']
    alpha = ntf_params['alpha']
    speed = ntf_params['speed']
    gei   = ntf_params['gei'  ]
    gii   = ntf_params['gii'  ]
    tauC  = ntf_params['tauC' ]

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

    # Cc = np.real(C*np.exp(-1j*Tau*w)).astype(float)
    Cc = C*np.exp(-1j*Tau*w)

    L1 = 0.8*np.identity(nroi)
    L2 = np.divide(1, np.sqrt(rowdegree*coldegree)+np.spacing(1))
    # diag(1./(sqrt(rowdegree.*coldegree)+eps));
    L = L1 - np.matmul(np.diag(L2), Cc)
    # L = np.array(L,dtype=np.float64)

    # try scipy.sparse.linalg.eigs next
    if use_smalleigs is True:
        d, v = np.linalg.eig(L)
        eig_ind = np.argsort(np.real(d))
        eig_vec = v[:, eig_ind]
        eig_val = d[eig_ind]
    else:
        d, v = np.linalg.eig(L)
        eig_ind = np.argsort(np.abs(d))
        eig_vec = v[:, eig_ind]
        eig_val = d[eig_ind]

    ev = np.transpose(eig_val[0:K])
    Vv = eig_vec[:, 0:K]  # why is eigv 1 all the same numbers?

    # Cortical model
    He = np.divide(1/tau_e**2, (1j*w+1/tau_e)**2)
    Hi = np.divide(gii*1/tau_i**2, (1j*w+1/tau_i)**2)

    Hed = alpha/tau_e/(1j*w + alpha/tau_e*He)
    Hid = alpha/tau_i/(1j*w + alpha/tau_i*Hi)

    Heid = gei*He*Hi/(1+gei*He*Hi)
    Htotal = a*Hed + (1-a)/2*Hid + (1-a)/2*Heid

    q1 = 1/alpha*tauC*(1j*w + alpha/tauC*He*ev)
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


def network_transfer_cost(params, C, D, lpf, FMEGdata, frange,
                          rois_with_MEG=np.arange(0, 68)):
    """Cost function for optimization of the model.

    Currently using negative of Pearson's correlation as cost metric.

    Args:
        params (numpy arr): 1x7 array of parameters: tau_e, tau_i, alpha,
        speed, gei, gii, tauC.
        C (numpy arr): Connectivity matrix.
        D (numpy arr): Distance matrix (introduce delay).
        lpf (numpy arr): low pass filter, designed before computing PSD.
        FMEGdata (numpy arr): Empirical data.
        frange (numpy arr): Vector of frequency bins for which the model
        compute a frequency response.
        rois_with_MEG (numpy arr): Regions with MEG. Usually excludes
        subcortical regions (default: np.arange(0,68)).

    Returns:
        err_out (float): Objective function evaluation result, negative of
        Pearson's correlation between empirica MEG and model result.

    """
    # network_transfer_function current inputs
    # (C, D, w, tau_e = 0.012, tau_i = 0.003, alpha = 1, speed = 5, gei = 4,\
    # gii = 1, tauC = 0.006)
    # defining parameters for the optimizer
    tau_e = params[0]
    tau_i = params[1]
    alpha = params[2]
    speed = params[3]
    gei = params[4]
    gii = params[5]
    tauC = params[6]

    # Computing model's frequency profiles
    freq_model = []
    err_min = np.zeros(rois_with_MEG.shape)
    for i in frange:
        w = 2*np.pi*i
        _, _, _, freqresp_out, _ = network_transfer_function(
                                                             C,
                                                             D,
                                                             w,
                                                             tau_e=tau_e,
                                                             tau_i=tau_i,
                                                             alpha=alpha,
                                                             speed=speed,
                                                             gei=gei,
                                                             gii=gii,
                                                             tauC=tauC
                                                             )
        freq_model.append(freqresp_out)

    freq_model = np.asarray(freq_model)
    freq_model = freq_model[:, rois_with_MEG].transpose()
    for n in rois_with_MEG:
        qdata = FMEGdata[n, :]
        if np.sum(qdata[:]) != 0:
            qdata = mag2db(qdata)
            qdata = qdata - np.mean(qdata)

        qmodel = np.abs(np.squeeze(freq_model[n, :]))
        qmodel = mag2db(np.convolve(qmodel, lpf, mode='same'))
        qmodel = qmodel - np.mean(qmodel)
        if np.sum(qmodel) == 0 or np.sum(qdata) == 0:
            err_min[n] = 0
        else:
            err_min[n] = pearsonr(qdata, qmodel)[0]

    err_out = -np.mean(err_min)
    return err_out
