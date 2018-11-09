import os
import matplotlib.pyplot as mpl
import numpy as np
import nitime.algorithms as tsa
#from scipy.io import loadmat, savemat
from scipy.signal import lfilter, firls, decimate
from scipy.stats import pearsonr
from scipy.optimize import basinhopping
from spectralgraph import *
"""
Script for running spectral graph model on MEG data from Ashish's initial code transfer.
"""
############################ Set Up ###############################
root_wd = '/home/axiezai/lab/spectral/'
#cwd = os.getcwd()
hcp_dir = os.path.join(root_wd,'data','HCPconnectomes')
# Subject name, data folder directories named after each subject
# MEG folder directory, and grab Julia's cortical brain region orders
sub_name = '8002.101'
MEGfolder = '/home/axiezai/lab/spectral/data/'

# Get structural connectivity matrix + distance between regions
Cdk_conn, Ddk_conn, permHCP = getHCPconn(hcp_dir)
permJulia, emptyJulia, cortJulia = Julia_order()

# Some other ordering that was in the original code:
linds = np.concatenate([np.arange(0,34), np.arange(68,77)]) #matlab: [0:34, 68:77]
rinds = np.concatenate([np.arange(34,68), np.arange(77,86)]) #matlab: [35:68 78:86]
#permHCP in matlab: [19:52, 53:86, 1:9, 10:18]  # permutes 86 x 86 HCP conn matrix so that subcorts are at the end, as in previous work

q = np.maximum(Cdk_conn[linds,:][:,linds], Cdk_conn[rinds,:][:,rinds])
q1 = np.maximum(Cdk_conn[linds,:][:,rinds], Cdk_conn[rinds,:][:,linds])

Cdk_conn = bi_symmetric_c(Cdk_conn, linds, rinds)
C = reduce_extreme_dir(Cdk_conn)

######################### Setting up MEG data stuff ################################
fs = 600 #sampling frequency
fmin = 2 # 2Hz - 45Hz signal range, filter for this with hbp
fmax = 45
fvec = np.linspace(fmin,fmax,40)
hbp = firls(101, np.array([0, 0.2*fmin, 0.9*fmin, fmax-2, fmax+5, 100])*2/fs,
           desired = np.array([0, 0, 1, 1, 0, 0])) #for detrending, a bandpass
lpf = np.array([1, 2, 5, 2, 1])
lpf = lpf/np.sum(lpf)
ind_del = hbp.size #number of coefficients in hbp. Delete that number in beginning of signal due to filtering

#Load MEG data
MEGdata, coords = getMEGdata(sub_name, cortJulia, MEGfolder)

# Get multi-taper spectra
print('Computing multi-taper power spectrum for all MEG data...')
FMEGdata = []
for row in MEGdata:
    q = lfilter(hbp, 1, row)
    q = q[ind_del:-1]
    ds_q = decimate(q, 4, axis = 0)
    f, psd, nu = tsa.multi_taper_psd(ds_q, Fs = fs/4, NW = 3, BW = 1, adaptive = False, jackknife = False)
    Fdata = np.convolve(psd, lpf, mode = 'same')
    FMEGdata.append(Fdata)
    
FMEGdata=np.asarray(FMEGdata)
assert FMEGdata.shape[0] == 68 #make sure we have 68 regions spectra

# plot between 2Hz and 45Hz, which was previously defined by fmin/fmax. Used for filtering
ind_fmin = np.abs(f-fmin).argmin()
ind_fmax = np.abs(f-fmax).argmin()
frange = f[ind_fmin:ind_fmax]
FMEGrange = FMEGdata[:,ind_fmin:ind_fmax]
#fig1, ax1 = mpl.subplots(1, 1)
# Plotting source localized MEG data
mpl.figure(num = 1, figsize=[6.5,3.8], dpi = 300)
mpl.xlabel('Frequency (Hz)')
mpl.ylabel('Magnitude (dB)')
for g in range(len(FMEGdata)):
    mpl.plot(frange,mag2db(FMEGrange[g,:]))
#mpl.show()

# Compute for all frequencies in fvec - this is the debug part in Ashish's SCFC_onJuliaMEG3.m
#freqresp = []
evec = []
Vvec =[]
fqall = []
freq_model = []

print('Computing model power spectrum for all given frequency (w)...')
for i in fvec:
    w = 2*np.pi*i
    fq, ev, Vv, freqresp_out, _ = NetworkTransferFunction(C, Ddk_conn, w)
    fqall.append(fq)
    evec.append(ev)
    Vvec.append(Vv)
    freq_model.append(freqresp_out)
    
fqall=np.asarray(fqall)
evec = np.asarray(evec)
Vvec = np.asarray(Vvec)
freq_model = np.asarray(freq_model)
np.moveaxis(Vvec,1,0).shape #86x40x57 just like matlab

# Plotting eigen vector's frequency response
ev_freqresp = np.abs(np.transpose(fqall))
fig2 = mpl.figure(num = 2, dpi = 300)
#fig_ev, ax_ev = mpl.subplot(1,1,1)
ax_ev = mpl.gca()
im1 = mpl.imshow(mag2db(ev_freqresp), extent = [0, 40, 0, 57], aspect='auto')
fig2.colorbar(im1)
ax_ev.set_title('Frequency Response of Each Eigen Vector')
#mpl.show()


############################# OPTIMIZATION TO MEG DATA ################################
## define upper and lower bounds for parameters:
# tau_e; tau_i; alpha; speed; gei; gii; tauC
# from scipy.optimize import basinhopping
VLB = np.array([0.005, 0.005, 0.1, 5.0, 0.5, 0.5,0.005])
VUB = np.array([0.020, 0.020, 1.0, 20, 5.0, 5.0, 0.020])
x0 = np.array([0.012, 0.003, 1, 5, 4, 1, 0.006])
minimizer_kwargs = {"method":"L-BFGS-B","args":(C, Ddk_conn, lpf, FMEGrange, frange)}
print('Starting basin hopping optimization...')

opt_output = basinhopping(networktransfer_costfun, x0, niter = 2, minimizer_kwargs = minimizer_kwargs, disp = True)
#save output from basinhopping
np.save('optm_params.npy',opt_output)

## Visualize optimized model spectra with MEG source localized spectra