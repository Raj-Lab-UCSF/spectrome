#generic modules
import numpy as np
import time
import argparse

#SCFC modules
from ..forward import network_transfer as nt
from ..utils import functions
from ..brain import Brain
from ..read import data_reader as dr
from ..preprocess import permute as perm
from ..utils import path as pth
from ..preprocess import filters
from scipy.signal import lfilter, firls, decimate

'''
    Script that parses input parameters and inputs them into
    Network transfer function to compute Spectral graph model
    outputs

'''



# define a brain from Brain() class and get connectome ordered as HCP
mybrain = Brain.Brain()
directory = pth.get_sibling_path('data')
mybrain.add_connectome(directory)
mybrain.reorder_connectome(mybrain.connectome, mybrain.distance_matrix)

# reordering to HCP
label_path = pth.get_sibling_path('atlases')
HCP_order = os.path.join(label_path, 'HCP_list.h5')
mybrain.add_ordering(HCP_order)

# reducing extreme values
mybrain.bi_symmetric_c()
mybrain.reduce_extreme_dir()


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

# DEFINE PARAMETER RANGE FOR SPEED, GEI, GII, AND TAUC
num_steps = 10
speed = np.linspace(5, 20, num_steps)
gei = np.linspace(0.5, 5, num_steps)
gii = np.linspace(0.5, 5.0, num_steps)
tauC = np.linspace(0.005, 0.020, num_steps)

#output = pd.DataFrame([], columns = ['Tau_e','Tau_i','Alpha','Speed','Gei','Gii','TauC','FreqModel'])

def make_filename(tau_e, tau_i, alpha, speediter, gei_iter, gii_iter, tauc_iter):
    filename  = "tauE_" +   '{0:.3f}'.format(tau_e, 3)
    filename += "_tauI_" +  '{0:.3f}'.format(tau_i, 3)
    filename += "_alpha_" + '{0:.3f}'.format(alpha, 3)
    filename += "_speed_" + '{0:.3f}'.format(speediter, 3)
    filename += "_gei_" +   '{0:.3f}'.format(gei_iter, 3)
    filename += "_gii_" +   '{0:.3f}'.format(gii_iter, 3)
    filename += "_tauC_" +  '{0:.3f}'.format(tauc_iter, 3) + '.h5'
    return filename

df_ind = 0
for speediter in speed:
    for gei_iter in gei:
        for gii_iter in gii:
            for tauc_iter in tauC:
                t0 = time.time()
                file_name = make_filename(tau_e,
                                          tau_i,
                                          alpha,
                                          speediter,
                                          gei_iter,
                                          gii_iter,
                                          tauc_iter)

                freq_model = []
                for i in fvec:
                    w = 2*np.pi*i
                    fq, ev, Vv, freqresp_out, _ = nt.network_transfer_function(mybrain,
                                                           parameters = mybrain.ntf_params,
                                                           w=2*np.pi)
                    freq_model.append(freqresp_out)

                freq_model = np.asarray(freq_model)

                # export
                data_dict = {"tau_e":tau_e, 'tau_i':tau_i, 'alpha':alpha, 'speed':speediter, 'gei':gei_iter, 'gii':gii_iter, 'tau_c':tauc_iter, 'freq':[freq_model]}
                file_path = os.path.join(directory, 'SCFC-data', file_name)
                pth.save_hdf5(file_path, data_dict)
                t1 = time.time()
                print('saved new file in:', str(t1-t0), 'seconds')
