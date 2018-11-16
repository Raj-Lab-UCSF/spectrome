import sys
sys.path.append("..")

import numpy as np
from forward import network_transfer as nt
from utils import path
from brain import brain

new_brain = brain.Brain()
hcp_dir = path.get_data_path()
new_brain.set_hcp_connectome(hcp_dir)

# Get structural connectivity matrix + distance between regions
Cdk_conn = new_brain.Cdk_conn
Ddk_conn = new_brain.Ddk_conn
perm_HCP = new_brain.permHCP

# Some other ordering that was in the original code:
linds = np.concatenate([np.arange(0,34), np.arange(68,77)])
rinds = np.concatenate([np.arange(34,68), np.arange(77,86)])

q = np.maximum(Cdk_conn[linds,:][:,linds], Cdk_conn[rinds,:][:,rinds])
q1 = np.maximum(Cdk_conn[linds,:][:,rinds], Cdk_conn[rinds,:][:,linds])

#not sure what this does.
new_brain.bi_symmetric_c(linds, rinds)
new_brain.reduce_extreme_dir(Cdk_conn)

C = new_brain.reducedC

default_ntf_params = {'tau_e':0.012, 'tau_i':0.003, 'alpha':1.0,
                       'speed':5.0, 'gei':4.0, 'gii':1.0, 'tauC':0.006}


fq, ev, Vv, freqresp_out, _ = nt.network_transfer_function(C=C, D=Ddk_conn, parameters = default_ntf_params, w=2*np.pi)


freqresp_out[0]





sg.network_transfer_function(C, Ddk_conn, w)
