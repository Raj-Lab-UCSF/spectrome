import sys
sys.path.append("..")

import numpy as np
import os

from brain import Brain
from utils import path as pth
from preprocess import preprocess as pp
from preprocess import permute as pm

def test_get_desikan():
    labelfile = 'desikan_atlas_68.csv'
    label_path = pth.get_sibling_path('dictionaries')
    label_filename = os.path.join(label_path, labelfile)
    regions, coords = pp.get_desikan(label_filename)
    assert regions[-1] == 'fusiform'

def test_get_HCP_order():
    labelfile = 'mean80_fibercount.csv'
    outlist = 'HCP_list.h5'
    label_path = pth.get_sibling_path('data')
    label_filename = os.path.join(label_path, labelfile)
    regions = pm.get_HCP_order(label_filename, save = False)
    assert regions[0] == 'LHbankssts'

def test_reorder_connectome():
    outlist = 'HCP_list.h5'
    path = pth.get_sibling_path('dictionaries')
    datapath = pth.get_sibling_path('data')
    confile = os.path.join(datapath, 'mean80_fibercount.csv')
    distfile = os.path.join(datapath, 'mean80_fiberlength.csv')
    filename = os.path.join(path, outlist)

    conn = np.genfromtxt(confile, delimiter=',', skip_header=1)
    dist = np.genfromtxt(distfile, delimiter=',', skip_header=0)

    _, _, perm = pm.reorder_connectome(conn, dist)
    assert np.array_equal(perm[:4], np.array([18, 19,20, 21]))
