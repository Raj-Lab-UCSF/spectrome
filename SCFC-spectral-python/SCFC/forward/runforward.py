''' running the ntf over a range of frequencies.'''
import sys, os
sys.path.append("..")
from forward import network_transfer as nt
import time as time
import numpy as np

def run_forward(brain, params, freqs):

    evec = []
    Vvec =[]
    fqall = []
    freq_model = []

    #start = time.time()
    for freq in freqs:
        w = 2*np.pi*freq
        fq, ev, Vv, freqresp_out, _ = nt.network_transfer_function(brain, params, w)
        fqall.append(fq)
        evec.append(ev)
        Vvec.append(Vv)
        freq_model.append(freqresp_out)

    fqall=np.asarray(fqall)
    evec = np.asarray(evec)
    Vvec = np.asarray(Vvec)
    freq_model = np.asarray(freq_model)
    freq_model = np.transpose(freq_model)
    #stop = time.time()
    #duration = stop - start

    #print('Computation time = ', duration)

    return freq_model
