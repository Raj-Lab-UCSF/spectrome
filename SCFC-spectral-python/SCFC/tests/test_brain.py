import sys
sys.path.append("..")

from brain import brain
from utils import path

def test_ntf_parameters():
    """Check whether a Brain class is being initialized correctly."""
    new_brain = Brain()
    default_ntf_params = {'tau_e':0.012, 'tau_i':0.003, 'alpha':1.0,
                           'speed':5.0, 'gei':4.0, 'gii':1.0, 'tauC':0.006}
    assert new_brain.ntf_parameters() == default_ntf_params

# def test_get_hcp_connectome():
#     """Check whether get_hcp_connectome returns expected values"""
#     new_brain = Brain()
#
#     hcp_dir = path.get_data_path()
#     conn = brain.Connectome()
#     get_hcp_connectome(hcp_dir)
#
#     CONN_ARR = np.array([18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
#     33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
#     50, 51,52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66,
#     67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83,
#     84, 85, 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14,
#     15, 16, 17])
#
#     assert conn[-1].all() == CONN_ARR.all()
