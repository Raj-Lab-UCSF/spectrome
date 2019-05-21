""" Script to apply the filtering step to the data and save out to dictoinary files --
do prior to further analysis as this is a slow step """

import numpy as np
import os
from ..preprocess import filters as ft
from ..utils import path as pth

#datapath = pth.get_sibling_path('data')
#MEG_path = os.path.join(datapath, 'MEG-chang')
MEG_path = '/Users/Megan/RajLab/MEG-chang'
directories = pth.walk_tree(MEG_path)

data_filename = 'DK_timecourse_20.h5'
out_data = 'DK_timecourse_20_filtered.h5'
out_freq = 'DK_timecourse_20_filterfreqs.h5'

for dir in directories:
    abspath = os.path.join(MEG_path,dir)
    data_path = os.path.join(abspath, data_filename)

    if os.path.isfile(data_path):
        #load data
        MEGdata = pth.read_hdf5(data_path)
        FMEGdata, f = ft.filter_MEG(MEGdata)

        pth.save_hdf5(os.path.join(abspath, out_data), FMEGdata)
        pth.save_hdf5(os.path.join(abspath, out_freq), f)

    else:
        print('folder ', abspath, ' is empty.')
