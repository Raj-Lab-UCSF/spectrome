import os, sys, argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as mpl
import nitime.algorithms as tsa
from spectralgraph import getHCPconn, Julia_order, mag2db, bi_symmetric_c, reduce_extreme_dir, NetworkTransferFunction
from scipy.signal import lfilter, firls, decimate
from scipy.io import savemat
'''
    Script that parses input parameters and inputs them into 
    Network transfer function to compute Spectral graph model
    outputs

'''
# creating an input argument parser object and parsing arguments
parser = argparse.ArgumentParser(description='Process spectral graph inputs')
#parser.add_argument("sub_name", type = str, help = "Anonymous subject name")
parser.add_argument("--tau_e", type = float, default = 0.012, help = "Excitatory time constant")
parser.add_argument("--tau_i", type = float, default = 0.003, help = "Inhibitory time constant")
parser.add_argument("--alpha", type = float, default = 1.0, help = "Alpha parameter")
#parser.add_argument("speed", type = float, default = 5.0, help = "Transmission velocity")
#parser.add_argument("gei", type = float, default = 4.0, help = "Parameter gain E/I")
#parser.add_argument("gii", type = float, default = 1.0, help = "Parameter gain I/I")
#parser.add_argument("tauC", type = float, default = 0.006, help = "Parameter Tau C")
args = parser.parse_args()
#sub_name = args.sub_name
tau_e = args.tau_e
tau_i = args.tau_i
alpha = args.alpha
#speed = args.speed
#gei = args.gei
#gii = args.gii
#tauC = args.tauC

## SET UP ## CHANGE THIS DIRECTORY FOR DIFFERENT PATHS ##
root_wd = '/home/axiezai/lab/spectral/'
hcp_dir = os.path.join(root_wd,'data','HCPconnectomes')
#sub_name = '8002.101'
#MEGfolder = '/home/axiezai/lab/spectral/data/'


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
fvec = np.linspace(fmin,fmax,200)
hbp = firls(101, np.array([0, 0.2*fmin, 0.9*fmin, fmax-2, fmax+5, 100])*2/fs,
           desired = np.array([0, 0, 1, 1, 0, 0])) #for detrending, a bandpass
lpf = np.array([1, 2, 5, 2, 1])
lpf = lpf/np.sum(lpf)
ind_del = hbp.size #number of coefficients in hbp. Delete that number in beginning of signal due to filtering

#Load MEG data
#MEGdata, coords = getMEGdata(sub_name, cortJulia, MEGfolder)

# DEFINE PARAMETER RANGE FOR SPEED, GEI, GII, AND TAUC
num_steps = 10
speed = np.linspace(5, 20, num_steps)
gei = np.linspace(0.5, 5, num_steps)
gii = np.linspace(0.5, 5.0, num_steps)
tauC = np.linspace(0.005, 0.020, num_steps)

#output = pd.DataFrame([], columns = ['Tau_e','Tau_i','Alpha','Speed','Gei','Gii','TauC','FreqModel'])
output = pd.DataFrame()
df_ind = 0
for speediter in speed:
    for gei_iter in gei:
        for gii_iter in gii:
            for tauc_iter in tauC:
                freq_model = []
                for i in fvec:
                    w = 2*np.pi*i
                    _, _, _, freqresp_out, _ = NetworkTransferFunction(C, Ddk_conn, w)
                    freq_model.append(freqresp_out)
                
                freq_model = np.asarray(freq_model)
                df = pd.DataFrame({"Tau_e":tau_e, 'Tau_i':tau_i, 'Alpha':alpha, 'Speed':speediter, 'Gei':gei_iter, 'Gii':gii_iter, 'TauC':tauc_iter, 'FreqModel':[freq_model]}, index = [df_ind])
                output = output.append(df, ignore_index = True)
                df_ind += 1

print(output.shape)