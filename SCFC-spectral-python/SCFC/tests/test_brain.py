import sys
sys.path.append("..")

from brain import Brain
from utils import path
import numpy as np

def test_ntf_parameters():
    """Check whether a Brain class is being initialized correctly."""
    new_brain = Brain.Brain()
    default_ntf_params = {'tau_e':0.012, 'tau_i':0.003, 'alpha':1.0,
                           'speed':5.0, 'gei':4.0, 'gii':1.0, 'tauC':0.006}
    assert new_brain.ntf_params == default_ntf_params

def test_set_hcp_connectome():
    new_brain = Brain.Brain()
    hcp_dir = path.get_data_path()
    new_brain.add_connectome(hcp_dir)
    new_brain.reorder_connectome(new_brain.connectome, new_brain.distance_matrix)
    assert np.array_equal(new_brain.connectome[0][:3],  np.array([0.   , 0.125, 3.15 ]))
