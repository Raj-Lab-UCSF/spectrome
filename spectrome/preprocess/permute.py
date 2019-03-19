'''Functions for rearranging matrices/data files'''

import sys, os
sys.path.append("..")
from utils import path as pth
import csv
import numpy as np

def order_dict(data, orderfile):
    """order_dict. Reorders a dictionary according to the order of keys passed in the
    list in orderfile. Returns a reordered dictionary.

    Args:
        data (type): Input data dictionary.
        orderfile (type): Name of file containing list of keys in the desired standard order.

    Returns:
        type: Dictionary of data, reordered to match the input standard of the orderfile.

    """
    dataorder = data.keys()
    standardlist = pth.read_hdf5(orderfile)
    newdata = {}
    loc_in_standard = []

    for key in standardlist:
        if key in dataorder:
            newdata[key] = data[key]
            loc_in_standard.append(standardlist.index(key))
        else:
            break
            # print('Skipping region of brain -- not in data')

    return newdata

def get_HCP_order(filepath, save=False, fileout = None, cortexstart = 18):
    """Import the HCP connectome, and create a list with the same order so
    that input data can be rearranged to compare. The dictionary keys are standardised to single
    words, lower case only. The HCP is also rearranged so that the cortex comes first,
    non-cortex parts of the brain are placed after.

    Args:
        filepath (string): Path to HCP connectome file.
        save (Boolean): Save output list to file?
        fileout (string): Location of output list.
        cortexstart (int): Index of start of cortex in original HCP ordering.

    Returns:
        List: List of brain regions; HCP rearranged so that the cortex comes first,
        non-cortex parts of the brain are placed after.


    """

    with open(filepath) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        region_order = []
        for row in csv_reader:
            if line_count == 0:
                for cell in row:
                    index = cell.find('ctx')
                    if index != -1:
                        cell = cell[4:6].upper() + cell[index+7:].lower()
                        region_order.append(cell)
                    else:
                        cell = cell
                        region_order.append(cell)
                line_count += 1
            else:
                break

    #put the non-cortex to the end of the order list
    order = region_order[cortexstart:]
    for item in region_order[0:cortexstart]:
        order.append(item)

    if save:
        pth.save_hdf5(fileout, order)

    return order

def reorder_connectome(conmat, distmat, save = False, cortexstart=18):
    """A function to rearrange matrices by a cyclical permutation (no rearranging of order).
    This is the equivalent of perm_HCP in the first code version: 
    np.concatenate([np.arange(18, 52),
                              np.arange(52, 86),
                              np.arange(0, 9),
                              np.arange(9, 18)])

    Args:
        conmat (numpy array): Direct input connectome
        distmat (numpy array): Direct input distance matrix
        save (bool): Whether to save out to files.
        cortexstart (int): Index of the first point in the cortex, eg. LHbankssts.

    Returns:
        numpy arrays: Connectome, distance matrix, and the permutation used on them.

    """
    Connectome = conmat
    Distance_matrix = distmat
    permutation = np.concatenate([np.arange(cortexstart, 86), np.arange(0, cortexstart)])
    Connectome = Connectome[permutation,][:, permutation]
    Distance_matrix = Distance_matrix[permutation, ][:, permutation]

    return Connectome, Distance_matrix, permutation
