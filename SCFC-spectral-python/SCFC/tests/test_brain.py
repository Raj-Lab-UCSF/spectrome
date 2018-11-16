import sys
sys.path.append("..")

from brain import brain
from utils import path

def test_ntf_parameters():
    """Check whether a Brain class is being initialized correctly."""
    new_brain = brain.Brain()
    default_ntf_params = {'tau_e':0.012, 'tau_i':0.003, 'alpha':1.0,
                           'speed':5.0, 'gei':4.0, 'gii':1.0, 'tauC':0.006}
    assert new_brain.ntf_parameters() == default_ntf_params

def test_set_hcp_connectome():
    new_brain = brain.Brain()
    hcp_dir = path.get_data_path()
    new_brain.set_hcp_connectome(hcp_dir)
    assert new_brain.Cdk_conn[0][:3] ==  np.array([0.   , 0.125, 3.15 ])
