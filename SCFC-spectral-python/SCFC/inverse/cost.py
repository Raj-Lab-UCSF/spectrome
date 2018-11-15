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
