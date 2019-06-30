""" running the ntf over a range of frequencies."""
from ..forward import network_transfer as nt
import time as time
import numpy as np


def run_forward(brain, params, freqs):
    """run_forward. Function for running the forward model over the passed in range of frequencies,
    for the handed set of parameters (which must be passed in as a dictionary)

    Args:
        brain (Brain): An instance of the Brain class.
        params (dict): Dictionary of a setting of parameters for the NTF model.
        freqs (array): Array of freqencies for which the model is to be calculated.

    Returns:
        array: Model values for each frequency, for each region of the brain, ordered as according to HCP
        (as in Brain class ordering).

    """

    evec = []
    Vvec = []
    fqall = []
    freq_model = []

    # start = time.time()
    for freq in freqs:
        w = 2 * np.pi * freq
        frequency_response, ev, Vv, model_out, _ = nt.network_transfer_function(brain, params, w)
        fqall.append(frequency_response)
        evec.append(ev)
        Vvec.append(Vv)
        freq_model.append(model_out)

    frequency_response = np.asarray(fqall)
    evec = np.asarray(evec)
    Vvec = np.asarray(Vvec)
    freq_model = np.asarray(freq_model)
    freq_model = np.transpose(freq_model)
    # stop = time.time()
    # duration = stop - start

    # print('Computation time = ', duration)

    return freq_model, frequency_response, evec, Vvec


def run_local_coupling_forward(brain, params, freqs):
    """run_forward. Function for running the forward model over the passed in range of frequencies,
    for the handed set of parameters (which must be passed in as a dictionary)

    Args:
        brain (Brain): An instance of the Brain class.
        params (dict): Dictionary of a setting of parameters for the NTF model.
        freqs (array): Array of freqencies for which the model is to be calculated.

    Returns:
        array: Model values for each frequency, for each region of the brain, ordered as according to HCP
        (as in Brain class ordering).

    """

    evec = []
    Vvec = []
    fqall = []
    freq_model = []

    # start = time.time()
    for freq in freqs:
        w = 2 * np.pi * freq
        frequency_response, ev, Vv, model_out, _ = nt.network_transfer_local_alpha(brain, params, w)
        fqall.append(frequency_response)
        evec.append(ev)
        Vvec.append(Vv)
        freq_model.append(model_out)

    frequency_response = np.asarray(fqall)
    evec = np.asarray(evec)
    Vvec = np.asarray(Vvec)
    freq_model = np.asarray(freq_model)
    freq_model = np.transpose(freq_model)
    # stop = time.time()
    # duration = stop - start

    # print('Computation time = ', duration)

    return freq_model, frequency_response, evec, Vvec

def run_HM_forward(brain, params, freqs):
    """Run the forward calculations with HM's local oscillators
    
    Args:
        same as above!
    """
    evec = []
    Vvec = []
    fqall = []
    freq_model = []
    Htotal_allfreq = []

    # start = time.time()
    for freq in freqs:
        w = 2 * np.pi * freq
        frequency_response, ev, Vv, model_out, Htotal = nt.network_transfer_HM(brain, params, w)
        fqall.append(frequency_response)
        evec.append(ev)
        Vvec.append(Vv)
        freq_model.append(model_out)
        Htotal_allfreq.append(Htotal)

    frequency_response = np.asarray(fqall)
    evec = np.asarray(evec)
    Vvec = np.asarray(Vvec)
    freq_model = np.asarray(freq_model)
    freq_model = np.transpose(freq_model)
    Htotal_out = np.asarray(Htotal_allfreq)

    return freq_model, frequency_response, evec, Vvec, Htotal_out