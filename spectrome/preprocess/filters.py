'''Functions for filtering, eg. removing low frequency component drift in data'''
################## MAY DEPRECATE THIS FILE IN THE FUTURE, get_freq_spectrum is moved
################## to preprocess.py, saving nodict version just in case we may need it.
from scipy.signal import lfilter, firls, decimate
import nitime.algorithms as tsa
import numpy as np
import matplotlib.pyplot as plt
from ..utils.functions import mag2db

def get_freq_spectrum(timeseries_data, fsampling, fmin= 2, fmax = 45, N = 68, downsample_factor = 4):
    """[Filters timeseries_data with a band pass filter designed with cutoff frequencies [fmin, fmax],
        the filtered time series will be downsampled by a factor of downsample_factor with a low-pass
        filter. The downsampled time series is then transformed into the frequency domain using the
        multi-taper method]

    Args:
        timeseries_data ([type]): [N x t time series data, with N brain regions for source localized data or
                                   N channels for sensor data. Duration = t time points, or t/fs seconds]
        fsampling ([type]): [sampling frequency of timeseries_data, no default given because this will vary]
        downsample_factor ([type]): Defaults to 4. [The ratio to downsample data, for example, 600 Hz MEG data
                                    will be downsampled to 150 Hz. This will be used in the decimate function,
                                    which has a low-pass filter built in to eliminate harmonic contamination]
        fmin (int, optional): Defaults to 2. [low cutoff frequency]
        fmax (int, optional): Defaults to 45. [high cutoff frequency]
        N (int, optional): Defaults to 68. [number of regions or number of channels]

    Returns:
        Freq_data[type]: [power spectrum for all input regions/channels]
        f : frequency vector/bins for power spectrum
    """

    fs = fsampling
    fvec = np.linspace(fmin,fmax,40)
    hbp = firls(101, np.array([0, 0.2*fmin, 0.9*fmin, fmax-2, fmax+5, 100])*2/fs,
           desired = np.array([0, 0, 1, 1, 0, 0])) #for detrending, a bandpass
    lpf = np.array([1, 2, 5, 2, 1])
    lpf = lpf/np.sum(lpf)
    ind_del = hbp.size #number of coefficients in hbp. Delete that number in beginning of signal due to filtering

    Freq_data = {}

    for key in list(timeseries_data.keys()):
        data = timeseries_data[key]
        data = data.astype(float)
        row = np.asarray(data)
        q = lfilter(hbp, 1, row)
        q = q[ind_del:-1] # delete transient portions
        ds_q = decimate(q, downsample_factor, axis = 0)
        f, psd, nu = tsa.multi_taper_psd(ds_q, Fs = fs/downsample_factor, NW = 3, BW = 1, adaptive = False, jackknife = False)
        Fdata = np.convolve(psd, lpf, mode = 'same')
        Freq_data[key] = Fdata

    assert len(Freq_data) == N #make sure we have correct number of regions/channels in the spectra

    # ind_fmin = np.abs(f-fmin).argmin()
    # ind_fmax = np.abs(f-fmax).argmin()
    # frange = f[ind_fmin:ind_fmax]
    # FMEGrange = Freq_data[:,ind_fmin:ind_fmax]
    return Freq_data, f

def get_freq_spectrum_notdict(timeseries_data, fsampling, fmin = 2, fmax = 45, N = 68, downsample_factor = 4, plot = True):
    """[When input data is not in a dict format. Filters timeseries_data with a band pass filter designed
        with cutoff frequencies [fmin, fmax], the filtered time series will be downsampled by a factor of
        downsample_factor with a low-pass filter. The downsampled time series is then transformed into the
        frequency domain using the multi-taper method]

    Args:
        timeseries_data ([type]): [N x t time series data, with N brain regions for source localized data or
                                   N channels for sensor data. Duration = t time points, or t/fs seconds]
        fsampling ([type]): [sampling frequency of timeseries_data, no default given because this will vary]
        downsample_factor ([type]): Defaults to 4. [The ratio to downsample data, for example, 600 Hz MEG data
                                    will be downsampled to 150 Hz. This will be used in the decimate function,
                                    which has a low-pass filter built in to eliminate harmonic contamination]
        fmin (int, optional): Defaults to 2. [low cutoff frequency]
        fmax (int, optional): Defaults to 45. [high cutoff frequency]
        N (int, optional): Defaults to 68. [number of regions or number of channels]
        plot (boolean, optional): Defaults to True. [Plot the spectra?]

    Returns:
        Freq_data[type]: [power spectrum for all input regions/channels]
        f : frequency vector/bins for power spectrum
    """
    fs = fsampling
    fvec = np.linspace(fmin,fmax,40)
    hbp = firls(101, np.array([0, 0.2*fmin, 0.9*fmin, fmax-2, fmax+5, 100])*2/fs,
           desired = np.array([0, 0, 1, 1, 0, 0])) #for detrending, a bandpass
    lpf = np.array([1, 2, 5, 2, 1])
    lpf = lpf/np.sum(lpf)
    ind_del = hbp.size #number of coefficients in hbp. Delete that number in beginning of signal due to filtering

    Freq_data = []
    for row in timeseries_data:
        q = lfilter(hbp, 1, row)
        q = q[ind_del:-1]
        ds_q = decimate(q, downsample_factor, axis = 0)
        f, psd, nu = tsa.multi_taper_psd(ds_q, Fs = fs/downsample_factor, NW = 3, BW = 1, adaptive = False, jackknife = False)
        Fdata = np.convolve(psd, lpf, mode = 'same')
        Freq_data.append(Fdata)

    Freq_data=np.asarray(Freq_data)
    assert Freq_data.shape[0] == N #make sure we have N regions/channels spectra

    ind_fmin = np.abs(f-fmin).argmin()
    ind_fmax = np.abs(f-fmax).argmin()
    frange = f[ind_fmin:ind_fmax]
    Freq_range = Freq_data[:,ind_fmin:ind_fmax]

    if plot == True:
        #fig1, ax1 = mpl.subplots(1, 1)
        # Plotting source localized MEG data
        plt.figure(num = 1, figsize=[3,1.5], dpi = 300)
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Magnitude (dB)')
        for g in range(len(Freq_data)):
            plt.plot(frange,mag2db(Freq_range[g,:]))
    else:
        print('No plot')

    output = (Freq_data, fvec, Freq_range, frange)

    return output
