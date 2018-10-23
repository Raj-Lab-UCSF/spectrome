import numpy as np
from scipy.io import loadmat
import os
from scipy.stats import pearsonr


"""
Functions used in the spectral graph model
"""

def getHCPconn(hcp_dir, conmat_in = 'mean80_fibercount.csv', dmat_in = 'mean80_fiberlength.csv'):
    '''[summary]
    
    Arguments:
        hcp_dir {str} -- [directory to HCP connectome. Example: '/home/user/data/HCPconnectomefolder']
    
    Keyword Arguments:
        conmat_in {str} -- [name of connectivity csv file] (default: {'mean80_fibercount.csv'})
        dmat_in {str} -- [name of fiber distance csv file] (default: {'mean80_fiberlength.csv'})
    
    Returns:
        [Cdk_conn] -- [Connectivity matrix arranged in appropriate order according to permHCP]
        [Ddk_conn] -- [Distance matrix arranged in permHCP order]
        [permHCP] -- [ordering of brain regions in DK86 atlas]
    '''

    cdk_hcp = np.genfromtxt(os.path.join(hcp_dir, conmat_in), delimiter=',', skip_header = 1)
    ddk_hcp = np.genfromtxt(os.path.join(hcp_dir, dmat_in), delimiter = ',', skip_header = 0)
    permHCP = np.concatenate([np.arange(18,52), np.arange(52,86), np.arange(0,9), np.arange(9,18)])
    Cdk_conn = cdk_hcp[permHCP,][:,permHCP]
    Ddk_conn = ddk_hcp[permHCP,][:,permHCP]
    return Cdk_conn, Ddk_conn, permHCP

def Julia_order():
    '''[Get brain region order for Julia Owen's data set, specific for DK86 atlas]
    
    Returns:
        [permJulia, emptyJulia, cortJulia] -- [Brain region orders for all regions, 
        regions with no MEG, and cortical regions respectively]
    '''

    cortJulia_lh = np.array([0, 1, 2, 3, 4,
                   6, 7, 8,
                   10, 11, 12, 13, 14, 15,
                   17, 16,
                   18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
                   5, 32, 33, 9])
    qsubcort_lh = np.array([0, 40, 36, 39, 38, 37, 35, 34, 0])
    qsubcort_rh = qsubcort_lh + 34 + 1
    cortJulia = np.concatenate([cortJulia_lh, 34 + cortJulia_lh])
    cortJulia_rh = cortJulia_lh + 34 + 7
    permJulia = np.concatenate([cortJulia_lh, cortJulia_rh, qsubcort_lh, qsubcort_rh])
    emptyJulia = np.array([68, 77, 76, 85])
    return permJulia, emptyJulia, cortJulia

def getMEGdata(sub_name, cortJulia, MEGfolder):
    '''[Get source localized MEG data and arrange data in correct ordering]
    
    Arguments:
        sub_name {str} -- [name of subject]
        cortJulia {numpy array} -- [Cortical region ordering]
        MEGfolder {str} -- [directory for MEG data]
    
    Returns:
        [MEGdata, coords] -- [MEG data as numpy array, coords = ?]
    '''

    S = loadmat(os.path.join(MEGfolder,sub_name,'DK_timecourse_20.mat'))
    MEGdata = S['DK_timecourse']
    MEGdata = MEGdata[cortJulia,]
    C = loadmat(os.path.join(MEGfolder,sub_name,'DK_coords_meg.mat'))
    coords = C['DK_coords_meg']
    coords = coords[cortJulia,]
    del S, C
    return MEGdata, coords

def mag2db(y):
    '''[convert magnitude response to decibels]
    
    Arguments:
        y {numpy array} -- [Power spectrum, raw magnitude response]
    
    Returns:
        [numpy array] -- [Power spectrum in dB]
    '''

    dby = 20*np.log10(y)
    return dby

def bi_symmetric_c(Cdk_conn, linds, rinds):
    q = np.maximum(Cdk_conn[linds,:][:,linds], Cdk_conn[rinds,:][:,rinds])
    q1 = np.maximum(Cdk_conn[linds,:][:,rinds], Cdk_conn[rinds,:][:,linds])
    Cdk_conn[np.ix_(linds,linds)] = q
    Cdk_conn[np.ix_(rinds,rinds)] = q
    Cdk_conn[np.ix_(linds,rinds)] = q1
    Cdk_conn[np.ix_(rinds,linds)] = q1
    return Cdk_conn

def reduce_extreme_dir(Cdk_conn, max_dir = 0.95, f=7):
    thr = f*np.mean(Cdk_conn[Cdk_conn > 0])
    C = np.minimum(Cdk_conn, thr)
    C = max_dir * C + (1-max_dir) * C
    return C

# def getMEGspectra(use_allJuliaMEG = False):
#    if use_allJuliaMEG:
#        S = loadmat(os.path.join(MEGfolder,'allFMEGdata.mat'))
#        FMEGdata = np.mean(S['allFMEGdata'], axis = 2)
#    else:
#        # Calculate spectra from source localized MEG
#        for row in MEGdata:

def NetworkTransferFunction(C, D, w, tau_e = 0.012, tau_i = 0.003, alpha = 1, speed = 5, gei = 4, gii = 1, tauC = 0.006):
    '''Network Transfer Function for spectral graph model.

    Arguments:
        C {numpy array} -- Connectivity matrix
        D {numpy array} -- Distance matrix
        w {float} -- frequency input. 

    Keyword Arguments:
        tau_e {float} -- Excitatory time constant parameter (default: {0.012})
        tau_i {float} -- Inhibitory time cosntant paramter (default: {0.003})
        alpha {float} -- [description] (default: {1})
        speed {float} -- Transmission velocity (default: {5})
        gei {float} -- Gain parameter (default: {4})
        gii {float} -- Gain parameter (default: {1})
        tauC {float} -- [description] (default: {0.006})

        Returns: freqresp, ev, Vv, freqresp_out, FCmodel
        freqresp {numpy array} -- 
        ev {numpy array} -- Eigen values
        Vv {numpy array} -- Eigen vectors
        freqresp_out {numpy array} -- Each region's frequency response for the given frequency (w)
        FCmodel {numpy array} -- Functional connectivity - still in the works
    '''

    # Not being used: Pin = 1 and tau_syn = 0.002
    
    # Defining some other parameters used: 
    zero_thr = 0.05
    use_smalleigs = True #otherwise uses full eig()
    numsmalleigs = np.round(2/3*C.shape[0]) #2/3
    a = 0.5 # fraction of signal at a node that is recurrent excitatory
    #gei = 4 # excitatory-inhibitory synaptic conductance as a ratio of E-E synapse
    #gii = 1 # inhibitory-inhibitory synaptic conductance as a ratio of E-E synapse
    #tauC = 0.5*tau_e
            
    rowdegree = np.transpose(np.sum(C,axis=1))
    coldegree = np.sum(C,axis=0)
    
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
    
    #Cc = np.real(C*np.exp(-1j*Tau*w)).astype(float)
    Cc = C*np.exp(-1j*Tau*w)
    
    L1 = 0.8*np.identity(nroi)
    L2 = np.divide(1,np.sqrt(rowdegree*coldegree)+np.spacing(1)) #diag(1./(sqrt(rowdegree.*coldegree)+eps));
    L = L1 - np.matmul(np.diag(L2),Cc)
    #L = np.array(L,dtype=np.float64)

    # try scipy.sparse.linalg.eigs next
    if use_smalleigs is True:
        d, v = np.linalg.eig(L)
        eig_ind = np.argsort(np.real(d))
        eig_vec = v[:,eig_ind]
        eig_val = d[eig_ind]
    else:
        d, v = np.linalg.eig(L)
        eig_ind = np.argsort(np.abs(d))
        eig_vec = v[:,eig_ind]
        eig_val = d[eig_ind]
        
    ev = np.transpose(eig_val[0:K])
    Vv = eig_vec[:,0:K] #why is eigv 1 all the same numbers?

    # Cortical model
    He = np.divide(1/tau_e**2,(1j*w+1/tau_e)**2)
    Hi = np.divide(gii*1/tau_i**2,(1j*w+1/tau_i)**2)

    Hed = alpha/tau_e/(1j*w + alpha/tau_e*He)
    Hid = alpha/tau_i/(1j*w + alpha/tau_i*Hi)

    Heid = gei*He*Hi/(1+gei*He*Hi)
    Htotal = a*Hed + (1-a)/2*Hid + (1-a)/2*Heid

    q1 = 1/alpha*tauC*(1j*w + alpha/tauC*He*ev)
    qthr = zero_thr*np.abs(q1[:]).max()
    magq1 = np.maximum(np.abs(q1),qthr)
    angq1 = np.angle(q1)
    q1 = np.multiply(magq1,np.exp(1j*angq1))
    freqresp = np.divide(Htotal,q1)

    freqresp_out = 0
    for k in range(1,K):
        freqresp_out += freqresp[k] * Vv[:,k]
    
    FCmodel = np.matmul(np.matmul(Vv[:,1:K], np.diag(freqresp[1:K]**2)), np.transpose(Vv[:,1:K]))

    den = np.sqrt(np.abs(freqresp_out))
    FCmodel = np.matmul(np.matmul(np.diag(1/den), FCmodel), np.diag(1/den))
    return freqresp, ev, Vv, freqresp_out, FCmodel

def networktransfer_costfun(params, C, D, lpf, FMEGdata, frange, rois_with_MEG = np.arange(0,68)):
    '''[Costfunction for optimization of the model, currently using pearon's correlation as metric]
    
    Arguments:
        params {numpy array} -- Currently using 7 parameters (1x7 array): tau_e, tau_i, alpha, speed, gei, gii, tauC
        C {numpy array} -- Connectivity matrix
        D {numpy array} -- Distance matrix to introduce delay
        lpf {numpy array} -- low pass filter, designed before computing PSD
        FMEGdata {numpy array} -- Empirical data 
        frange {numpy array} -- Vector of frequency bins for which the model computes a frequency response.
    
    Keyword Arguments:
        rois_with_MEG {numpy array} -- which regions has MEG, only optimize these regions. Usually exclude 
        subcortical regions (default: {np.arange(0,68)})
    
    Returns:
        float -- Should be a float? Objective function evaluation result, negative of Pearson's correlation between empirical
        MEG and model result because we want to maximize. 
    '''

    # NetworkTransferFunction current inputs
    # (C, D, w, tau_e = 0.012, tau_i = 0.003, alpha = 1, speed = 5, gei = 4, gii = 1, tauC = 0.006)
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
        _, _, _, freqresp_out, _ = NetworkTransferFunction(C, D, w, tau_e = tau_e, tau_i = tau_i, 
                                                          alpha = alpha, speed = speed, gei = gei, 
                                                          gii = gii, tauC = tauC)
        freq_model.append(freqresp_out)
        
    freq_model = np.asarray(freq_model)
    freq_model = freq_model[:,rois_with_MEG].transpose()
    for n in rois_with_MEG:
        qdata = FMEGdata[n,:]
        if np.sum(qdata[:]) != 0:
            qdata = mag2db(qdata)
            qdata = qdata - np.mean(qdata)
            
        qmodel = np.abs(np.squeeze(freq_model[n,:]))
        qmodel = mag2db(np.convolve(qmodel, lpf, mode = 'same'))
        qmodel = qmodel - np.mean(qmodel)
        if np.sum(qmodel) == 0 or np.sum(qdata) == 0:
            err_min[n] = 0
        else:
            err_min[n] = pearsonr(qdata,qmodel)[0]
            
    err_out = -np.mean(err_min)
    return err_out