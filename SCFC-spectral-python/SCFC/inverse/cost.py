from scipy.stats import pearsonr
import numpy as np
import math

def pearson(data, freq_model):
    """Pearson. Calculate Pearson r correlation for each member of data dict against freq model.

    Args:
        data (dict): Data dictionary.
        freq_model (array): frequency model data in corresponding order

    Returns:
        float: Pearsonr correlation mean, and list of values over regions.

    """
    i = 0
    pearsonr_list = []
    for key in data.keys():
        modregion = freq_model[i,:]
        region = data[key]
        region = [float(x) for x in region]#we need to do this to make sure things are in the correct form for scipy.stats
        modregion = [float(x) for x in modregion]
        pearson = pearsonr(region,modregion)[0]
        pearsonr_list.append(pearson)
        i += 1
    pearson = sum(pearsonr_list)
    return pearson, pearsonr_list


def pearson_cost_oneminus(data, freq_model):
    """pearson_cost_oneminus. A cost function based on the pearson r correlation. Designed to
    be minimised in a fitting process. Uses '1-r' as a measure of data-model error.

    Args:
        data (dict): Data dictionary.
        freq_model (array): frequency model data in corresponding order

    Returns:
        floats: mean error based on pearson r, and a list of values over regions.

    """
    i = 0
    err_list = []
    for key in data.keys():
        modregion = freq_model[i,:]
        region = data[key]
        region = [float(x) for x in region]#we need to do this to make sure things are in the correct form for scipy.stats
        modregion = [float(x) for x in modregion]
        err = pearsonr(region,modregion)[0]
        err = 1 - err #to convert to an error to minimise.
        err_list.append(err)
        i += 1
    err = sum(err_list)
    return err, err_list

def pearson_cost_exp(data, freq_model):
    """pearson_cost. A cost function based on the pearson r correlation. Designed to
    be minimised in a fitting process, using e^(-r) as a measure of data-model distance.

    Args:
        data (dict): Data dictionary.
        freq_model (array): frequency model data in corresponding order

    Returns:
        floats: mean error based on pearson r, and a list of values over regions.

    """
    i = 0
    err_list = []
    for key in data.keys():
        modregion = freq_model[i,:]
        region = data[key]
        region = [float(x) for x in region]#we need to do this to make sure things are in the correct form for scipy.stats
        modregion = [float(x) for x in modregion]
        err = pearsonr(region,modregion)[0]
        err = math.exp(-1*err) #to convert to an error to minimise.
        err_list.append(err)
        i += 1
    err = sum(err_list)
    return err, err_list

def pearson_cost_dB(data, model):
    i = 0
    err_list = []
    for key in data.keys():
        modregion = functions.mag2db(freq_model[i,:])
        region = functions.mag2db(data[key])
        region = [float(x) for x in region]#we need to do this to make sure things are in the correct form for scipy.stats
        modregion = [float(x) for x in modregion]
        err = pearsonr(region,modregion)[0]
        err = 1 - err #to convert to an error to minimise.
        err_list.append(err)
        i += 1
    err = sum(err_list)
    return err, err_list

def network_transfer_cost(params, C, D, lpf, FMEGdata, frange,
                          rois_with_MEG=np.arange(0, 68)):
    """Cost function for optimization of the model (old form).

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
        Pearson's correlation between empirical MEG and model result.

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

    #demean data
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
