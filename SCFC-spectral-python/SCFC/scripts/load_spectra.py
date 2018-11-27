import sys,os
sys.path.append("..")

# #SCFC modules
# from forward import network_transfer as nt
# from utils import functions
# from brain import Brain
# from read import data_reader as dr
# from preprocess import permute as perm
from utils import path as pth
# from preprocess import filters
#
# #generic modules
# import numpy as np
# from scipy.signal import lfilter, firls, decimate
# import time
# import argparse


directory = pth.get_sibling_path('data')
file_path = os.path.join(directory, 'SCFC-data')
file_name = os.path.join(file_path, 'tauE_0.012_tauI_0.003_alpha_1.000_speed_5.000_gei_0.500_gii_0.500_tauC_0.005.h5')
