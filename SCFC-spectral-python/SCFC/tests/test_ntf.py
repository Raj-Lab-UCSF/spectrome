import sys
sys.path.append("..")

#SCFC modules
from forward import network_transfer as nt
from utils import functions
from brain import Brain
from read import data_reader as dr
from preprocess import permute as perm
from utils import path as path
# from preprocess import filters

#generic modules
import numpy as np
from scipy.signal import lfilter, firls, decimate
import time


my_brain = Brain.Brain()

data_dir = path.get_sibling_path('data')

my_brain.add_connectome(data_dir)

my_brain.reorder_connectome(my_brain.connectome, my_brain.distance_matrix)

def test_distance_matrix():
    assert my_brain.distance_matrix.shape == (86, 86)
