import sys
sys.path.append("..")

import numpy as np
import os

from brain import Brain
from utils import path as pth
from preprocess import preprocess as pp
from preprocess import permute as pm

## This bit processes the data to produce dictionary format versions.
datapath = pth.get_sibling_path('data')
MEG_path = os.path.join(datapath, 'MEG-chang')

directories = pth.walk_tree(MEG_path)
coord_filename = 'DK_coords_meg.mat'
data_filename = 'DK_timecourse_20.mat'
out_coords = 'DK_coords_meg.h5'
out_data = 'DK_timecourse_20.h5'

labelfile = 'OrderingAlphabetical_68ROIs.txt'
label_path = pth.get_sibling_path('dictionaries')
label_filename = os.path.join(label_path, labelfile)

for dir in directories:
    abspath = os.path.join(MEG_path,dir)
    coord_path = os.path.join(abspath, coord_filename)
    data_path = os.path.join(abspath, data_filename)

    if os.path.isfile(data_path):
        data_dict = pp.add_key_data(label_filename, data_path)
        pth.save_hdf5(os.path.join(abspath, out_data), data_dict)

        coord_dict = pp.add_key_coords(label_filename, coord_path)
        pth.save_hdf5(os.path.join(abspath, out_coords), coord_dict)

    else:
        print('folder ', abspath, ' is empty.')
