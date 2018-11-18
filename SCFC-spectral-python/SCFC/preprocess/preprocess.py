'''functions to make sure that data is in the right order, or that data is given dictionary
keys labelling the data according to the brain regions'''

import sys, os
sys.path.append("..")
from utils import path
from scipy.io import loadmat
import re
import h5py
import deepdish as dd

def key_data(label_filepath, data_filepath):
    ''' Transfer unlabelled matrix data in to a dictionary format, using
    labels taken from the label file (ideally this should already be the case,
    this is a desperate hack.)
    '''

    label_file = open(label_filepath, "r")
    lines = label_file.readlines()
    label_file.close()

    #hack for Chang's data -- cleaning up ROIs list format -- can change for other versions
    #of data
    i = 0
    for line in lines:
        index_stop = line.find('.')
        ind_newline = line.find('\n')
        lines[i] = line[index_stop+1:ind_newline].lower()
        i += 1

    #import the data and apply the list members as keys, resave data in better format
    data = loadmat(data_filepath)
    data = data['DK_timecourse']

    data_dict ={}
    for keys, values in zip(lines, data[:,:]):
        data_dict[keys] = values

    return data_dict

def save_hdf5(path, dict):
    '''save out in to a dictionary format, using HDF5 format and deepdish package'''
    dd.io.save(path, dict)

def read_hdf5(path):
    '''read in a hdf5 file'''
    dict = dd.io.load(path)
    return dict


if __name__ == "__main__":

    data_path = path.get_sibling_path('data')
    for (data_path, dirs, files) in os.walk(data_path):
        print(files)

    # labelfile = 'OrderingAlphabetical_68ROIs.txt'
    # datafile = '8002.101/DK_timecourse_20.mat'
    # outfile = '8002.101/DK_timecourse_20_test.h5'
    #
    # data_path = path.get_sibling_path('data')
    # data_filename = os.path.join(data_path, datafile)
    #
    # label_path = path.get_sibling_path('dictionaries')
    # label_filename = os.path.join(label_path, labelfile)
    #
    # outdict = key_data(label_filename, data_filename)
    #
    # # print(outdict.keys())
    #
    # outpath = os.path.join(data_path, outfile)
    # # save_dict(outpath, outdict)
    #
    # print(read_hdf5(outpath)['bankssts'][0])
