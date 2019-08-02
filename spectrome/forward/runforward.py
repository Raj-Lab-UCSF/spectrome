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

    eigenvalues = []
    eigenvectors = []
    frequency_response = []
    model_out = []

    for freq in freqs:
        w = 2 * np.pi * freq
        freq_resp, eig_val, eig_vec, freq_model, _ = nt.network_transfer_function(brain, params, w)
        frequency_response.append(freq_resp)
        eigenvalues.append(eig_val)
        eigenvectors.append(eig_vec)
        model_out.append(freq_model)

    frequency_response = np.asarray(frequency_response)
    eigenvalues = np.asarray(eigenvalues)
    eigenvectors = np.asarray(eigenvectors)
    model_out = np.transpose(np.asarray(model_out))

    return model_out, frequency_response, eigenvalues, eigenvectors


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

    eigenvalues = []
    eigenvectors = []
    frequency_response = []
    model_out = []

    for freq in freqs:
        w = 2 * np.pi * freq
        freq_resp, eig_val, eig_vec, freq_model, _ = nt.network_transfer_local_alpha(brain, params, w)
        frequency_response.append(freq_resp)
        eigenvalues.append(eig_val)
        eigenvectors.append(eig_vec)
        model_out.append(freq_model)

    frequency_response = np.asarray(frequency_response)
    eigenvalues = np.asarray(eigenvalues)
    eigenvectors = np.asarray(eigenvectors)
    model_out = np.transpose(np.asarray(model_out))

    return model_out, frequency_response, eigenvalues, eigenvectors

def run_HM_forward(brain, params, freqs):
    """Run the forward calculations with HM's local oscillators
    
    Args:
        same as above!
    """
    eigenvalues = []
    eigenvectors = []
    frequency_response = []
    model_out = []
    Htotal_allfreq = []

    for freq in freqs:
        w = 2 * np.pi * freq
        freq_resp, eig_val, eig_vec, freq_model, Htotal = nt.network_transfer_HM(brain, params, w)
        frequency_response.append(freq_resp)
        eigenvalues.append(eig_val)
        eigenvectors.append(eig_vec)
        model_out.append(freq_model)
        Htotal_allfreq.append(Htotal)

    frequency_response = np.asarray(frequency_response)
    eigenvalues = np.asarray(eigenvalues)
    eigenvectors = np.asarray(eigenvectors)
    model_out = np.transpose(np.asarray(model_out))
    Htotal_frequency_response = np.asarray(Htotal_allfreq)

    return model_out, frequency_response, eigenvalues, eigenvectors, Htotal_frequency_response