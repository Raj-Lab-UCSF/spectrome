# SCFC modules
from ..forward import network_transfer as nt
from ..utils import functions
from ..brain import Brain
from ..read import data_reader as dr
from ..preprocess import permute as perm
from ..utils import path as path
from scipy.signal import lfilter, firls, decimate

# generic modules
import numpy as np
import time


my_brain = Brain.Brain()

data_dir = path.get_sibling_path("data")

my_brain.add_connectome(data_dir)

my_brain.reorder_connectome(my_brain.connectome, my_brain.distance_matrix)


def test_distance_matrix():
    assert my_brain.distance_matrix.shape == (86, 86)


my_brain.bi_symmetric_c()
my_brain.reduce_extreme_dir()


def test_weird_functions():
    assert my_brain.reducedConnectome[0][1] == 0.125


def test_ntf():
    fq, ev, Vv, f, _ = nt.network_transfer_function(
        my_brain, parameters=my_brain.ntf_params, w=2 * np.pi
    )


def test_run_forward():
    parameters = my_brain.ntf_params
    fvec = np.array([40.0])
    freq_model, frequency_response, evec, Vvec = rf.run_forward(
        my_brain, parameters, fvec
    )
    assert freq_model.shape == (86, 1)
    assert frequency_response.shape == (86, 1)
    assert evec.shape == (1, 86)
    assert Vvec.shape == (1, 86, 57)
