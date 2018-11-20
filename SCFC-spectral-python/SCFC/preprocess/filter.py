'''Functions for filtering, eg. removing low frequency component drift in data'''
import sys, os
sys.path.append("..")
from scipy.signal import lfilter, firls, decimate
import nitime.algorithms as tsa
import numpy as np
import matplotlib.pyplot as plt
from utils.functions import mag2db


def filter_MEG(MEGdata, fsampling=600, fmin=2, fmax=45, regions =68, plot = True):
    fs = fsampling
    fvec = np.linspace(fmin,fmax,40)
    hbp = firls(101, np.array([0, 0.2*fmin, 0.9*fmin, fmax-2, fmax+5, 100])*2/fs,
           desired = np.array([0, 0, 1, 1, 0, 0])) #for detrending, a bandpass
    lpf = np.array([1, 2, 5, 2, 1])
    lpf = lpf/np.sum(lpf)
    ind_del = hbp.size #number of coefficients in hbp. Delete that number in beginning of signal due to filtering

    FMEGdata = []
    for row in MEGdata:
        q = lfilter(hbp, 1, row)
        q = q[ind_del:-1]
        ds_q = decimate(q, 4, axis = 0)
        f, psd, nu = tsa.multi_taper_psd(ds_q, Fs = fs/4, NW = 3, BW = 1, adaptive = False, jackknife = False)
        Fdata = np.convolve(psd, lpf, mode = 'same')
        FMEGdata.append(Fdata)

    FMEGdata=np.asarray(FMEGdata)
    assert FMEGdata.shape[0] == regions #make sure we have 68 regions spectra

    if plot == True:
        ind_fmin = np.abs(f-fmin).argmin()
        ind_fmax = np.abs(f-fmax).argmin()
        frange = f[ind_fmin:ind_fmax]
        FMEGrange = FMEGdata[:,ind_fmin:ind_fmax]
        #fig1, ax1 = mpl.subplots(1, 1)
        # Plotting source localized MEG data
        plt.figure(num = 1, figsize=[6.5,3.8], dpi = 300)
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Magnitude (dB)')
        for g in range(len(FMEGdata)):
            plt.plot(frange,mag2db(FMEGrange[g,:]))
    else:
        break

    return FMEGdata, fvec
