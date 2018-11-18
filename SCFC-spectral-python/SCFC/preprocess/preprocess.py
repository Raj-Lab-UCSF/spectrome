'''functions to make sure that data is in the right order, or that data is given dictionary
keys labelling the data according to the brain regions'''

import sys, os
sys.path.append("..")
from utils import path
from scipy.io import loadmat
import re
import h5py
import deepdish as dd

def add_key_data(label_filepath, data_filepath):
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

def add_key_coords(label_filepath, coordfile):
    label_file = open(label_filepath, "r")
    lines = label_file.readlines()
    label_file.close()

    i = 0
    for line in lines:
        index_stop = line.find('.')
        ind_newline = line.find('\n')
        lines[i] = line[index_stop+1:ind_newline].lower()
        i += 1

    coords = loadmat(coordfile)
    coords = coords['DK_coords_meg']

    coord_dict ={}
    for keys, values in zip(lines, coords):
        coord_dict[keys] = values
    return coord_dict




def save_hdf5(path, dict):
    '''save out in to a dictionary format, using HDF5 format and deepdish package'''
    dd.io.save(path, dict)

def read_hdf5(path):
    '''read in a hdf5 file'''
    dict = dd.io.load(path)
    return dict

def walk_tree(datapath):
    directories = []
    for (path, dirs, files) in os.walk(datapath):
        directories.append(dirs)
    return directories[0]



if __name__ == "__main__":

    datapath = '/Users/Megan/RajLab/MEG-chang'
    directories = walk_tree(datapath)
    coord_filename = 'DK_coords_meg.mat'
    data_filename = 'DK_timecourse_20.mat'
    out_coords = 'DK_coords_meg.h5'
    out_data = 'DK_timecourse_20.h5'

    labelfile = 'OrderingAlphabetical_68ROIs.txt'
    label_path = path.get_sibling_path('dictionaries')
    label_filename = os.path.join(label_path, labelfile)

    for dir in directories:
        abspath = os.path.join(datapath,dir)
        coord_path = os.path.join(abspath, coord_filename)
        data_path = os.path.join(abspath, data_filename)

        data_dict = add_key_data(label_filename, data_path)
        save_hdf5(os.path.join(abspath, out_data), data_dict)

        coord_dict = add_key_coords(label_filename, coord_path)
        save_hdf5(os.path.join(abspath, out_coords), coord_dict)



    
