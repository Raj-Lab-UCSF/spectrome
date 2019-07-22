"""functions to denoise, downsample, re-order, label, and transform the input functional
time series data. These are applied prior to any further processing."""

import csv
import numpy as np
import nitime.algorithms as tsa
import matplotlib
matplotlib.use('PS')
import matplotlib.pyplot as plt

from ..utils import path as pth
from ..utils.functions import mag2db
from scipy.io import loadmat
from scipy.signal import lfilter, firls, decimate


def add_key_to_matdata(label_filepath, data):
    """Add dictionary keys (brain regions) to MEG data from raw MAT files
    (note that this currently works for one particularly style of brain region input,
    and will be generalised in future). The ordering used here is the desikan ordering.

    Args:
        label_filepath (type): full path and filename of the file containing the
        correct order for the data being passed.
        data (array): numpy array of the data in non-dictionary form.

    Returns:
        type: a dictionary of the data, keyed by brain region.

    """
    label_file = open(label_filepath, "r")
    lines = label_file.readlines()
    label_file.close()

    # hack for Chang's data -- cleaning up ROIs list format -- can change for other versions
    # of data
    i = 0
    for line in lines:
        index_stop = line.find(".")
        ind_newline = line.find("\n")
        lines[i] = line[0:2].upper() + line[index_stop + 1 : ind_newline].lower()
        i += 1

    # import the data and apply the list members as keys, resave data in better format

    data_dict = {}
    for keys, values in zip(lines, data[:, :]):
        data_dict[keys] = values

    return data_dict


def add_key_data(label_filepath, data_filepath):
    """Add dictionary keys (brain regions) to MEG data from raw MAT files
    (note that this currently works for one particularly style of brain region input,
    and will be generalised in future). The ordering used here is the desikan ordering.

    Args:
        label_filepath (type): full path and filename of the file containing the
        correct order for the data being passed.
        data_filepath (type): full path and filename of data MAT file.

    Returns:
        type: a dictionary of the data, keyed by brain region.

    """

    label_file = open(label_filepath, "r")
    lines = label_file.readlines()
    label_file.close()

    # hack for Chang's data -- cleaning up ROIs list format -- can change for other versions
    # of data
    i = 0
    for line in lines:
        index_stop = line.find(".")
        ind_newline = line.find("\n")
        lines[i] = line[0:2].upper() + line[index_stop + 1 : ind_newline].lower()
        i += 1

    # import the data and apply the list members as keys, resave data in better format
    data = loadmat(data_filepath)
    data = data["DK_timecourse"]

    data_dict = {}
    for keys, values in zip(lines, data[:, :]):
        data_dict[keys] = values

    return data_dict


def add_key_coords(label_filepath, coordfile):
    """Add dictionary keys (brain regions) to brain region coordinates from raw MAT files
    (note that this currently works for one particularly style of brain region input,
    and will be generalised in future).

    Args:
        label_filepath (type): full path and filename of the file containing the
        correct order for the data being passed.
        coordfile (type): full path and filename of coord MAT file.

    Returns:
        type: dictionary of coordinates, keyed by brain region.

    """
    label_file = open(label_filepath, "r")
    lines = label_file.readlines()
    label_file.close()

    i = 0
    for line in lines:
        index_stop = line.find(".")
        ind_newline = line.find("\n")
        lines[i] = line[0:2].upper() + line[index_stop + 1 : ind_newline].lower()
        i += 1

    coords = loadmat(coordfile)
    coords = coords["DK_coords_meg"]

    coord_dict = {}
    for keys, values in zip(lines, coords):
        coord_dict[keys] = values
    return coord_dict


def denoise_timeseries(timeseries_data, fsampling, fmin=2, fmax=45):
    """Filters timeseries_data with a band pass filter designed with cutoff frequencies [fmin, fmax]

    Args:
        timeseries_data (dict): [N x t time series data, with N brain regions for source localized data or
                                   N channels for sensor data. Duration = t time points, or t/fs seconds]
        fsampling ([type]): [sampling frequency of timeseries_data, no default given because this will vary]
        fmin (int, optional): Defaults to 2. [low cutoff frequency]
        fmax (int, optional): Defaults to 45. [high cutoff frequency]

    Returns:
        clean_data (dict): [denoised time series data]
    """

    # fvec = np.linspace(fmin,fmax,Nbins)
    hbp = firls(
        101,
        np.array([0, 0.2 * fmin, 0.9 * fmin, fmax - 2, fmax + 5, 100]) * 2 / fsampling,
        desired=np.array([0, 0, 1, 1, 0, 0]),
    )  # for detrending, a bandpass
    lpf = np.array([1, 2, 5, 2, 1])
    lpf = lpf / np.sum(lpf)
    ind_del = (
        hbp.size
    )  # number of coefficients in hbp. Delete that number in beginning of signal due to filtering
    clean_data = {}
    for key in list(timeseries_data.keys()):
        data = timeseries_data[key]
        data = data.astype(float)
        row = np.asarray(data)
        q = lfilter(hbp, 1, row)
        q = q[ind_del:-1]  # delete transient portions from filter
        clean_data[key] = q

    return clean_data


def downsample_data(data, fsampling, downsample_factor=4):
    """Using the decimate function to low-pass filter and downsample

    Args:
        data (dict): [Data to be downsampled, doesn't necessarily have to be time course]
        fsampling (int): [sampling frequency of data]
        downsample_factor (int, optional): Defaults to 4. [The ratio to downsample data, for example, 600 Hz MEG data
                                    will be downsampled to 150 Hz. This will be used in the decimate function,
                                    which has a low-pass filter built in to eliminate harmonic contamination]

    Returns:
        DS_data (dict): [downsampled data]
        DS_fs (int): [sampling frequency for downsampled data]
    """

    DS_data = {}
    DS_fs = fsampling / downsample_factor

    for key in list(data.keys()):
        heavy_data = np.asarray(data[key].astype(float))
        DS_data[key] = decimate(heavy_data, downsample_factor, axis=0)

    return DS_data, DS_fs


def get_freq_spectrum(timeseries_data, fsampling, fmin, fmax, plot=True):
    """Transform a clean, downsampled time series data into frequency domain
    using the multi-taper method

    Args:
        timeseries_data ([type]): [N x t time series data, with N brain regions for source localized data or
                                   N channels for sensor data. Duration = t time points, or t/fs seconds]
        fsampling ([type]): [sampling frequency of timeseries_data, no default given because this will vary]
        fmin (int, optional): Defaults to 2. [low cutoff frequency]
        fmax (int, optional): Defaults to 45. [high cutoff frequency]
        plot (boolean, optional): Defaults to True. [Plot the spectra?]

    Returns:
        freq_data (dict): [power spectrum for all input regions/channels]
        frange: frequency vector ranging from fmin to fmax
        Freq_range: frequency magnitudes within frange
    """

    freq_data = {}
    lpf = np.array([1, 2, 5, 2, 1])
    lpf = lpf / np.sum(lpf)
    freq_data = {}

    for key in list(timeseries_data.keys()):
        timedata = np.asarray(timeseries_data[key].astype(float))
        f, psd, nu = tsa.multi_taper_psd(
            timedata, Fs=fsampling, NW=3, BW=1, adaptive=False, jackknife=False
        )
        Fdata = np.convolve(psd, lpf, mode="same")
        freq_data[key] = Fdata

    ind_fmin = np.abs(f - fmin).argmin()
    ind_fmax = np.abs(f - fmax).argmin()
    frange = f[ind_fmin:ind_fmax]
    Freq_range = Freq_data[:, ind_fmin:ind_fmax]
    if plot == True:
        # fig1, ax1 = mpl.subplots(1, 1)
        # Plotting source localized MEG data
        plt.figure(num=1, figsize=[3, 1.5], dpi=300)
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Magnitude (dB)")
        for g in range(len(Freq_data)):
            plt.plot(frange, mag2db(Freq_range[g, :]))
    else:
        print("No plot")

    return freq_data, frange, Freq_range


## May get moved to another folder in near future :)
def get_desikan(filepath):
    """Parse the desikan atlas ordering for the 68 cortical regions from the csv file,
    or indeed any other csv dictionary.

    Args:
        filepath (type): full file path.

    Returns:
        regions: list of region names in order of the desikan atlas.
        coords: the coordinates of the desikan atlas, in the same order.

    """
    with open(filepath) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=",")
        line_count = 0
        regions = []
        coords = []
        for row in csv_reader:
            if line_count == 0:
                print(f'Column names are {", ".join(row)}')
                line_count += 1
            else:
                regions.append(row[1])
                x = float(row[2])
                y = float(row[3])
                z = float(row[4])
                coords.append([x, y, z])
    return regions, coords
