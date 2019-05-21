'''Collection of functions to manipulate atlases and ROIs'''

import numpy as np
import os

from ..utils.path import get_root_path

def get_DK_names(cortical_only=False):
    """Return the DK regions in alphabetical ordering.

    Args:
        cortical_only (Bool): If true, returns 68 cortical regions instead of 86.

    Returns:
        array: names of DK regions.

    """
    if cortical_only:
        names_file = 'OrderingAlphabetical_68ROIs.txt'
    else:
        names_file = 'OrderingAlphabetical_86ROIs.txt'

    root_path = get_root_path()
    file_path = os.path.join(root_path,'atlases/')
    file_path = os.path.join(file_path, names_file)
    names_array = np.genfromtxt(file_path,dtype='str')
    return(names_array)

def get_HCP_names(cortical_only=False):
    """Return the DK regions following yet another naming convention.

    Args:
        cortical_only (Bool): If true, returns 68 cortical regions instead of 86.

    Returns:
        array: names of DK regions following convention 2.

    """
    if cortical_only:
        names_file = 'OrderingAlphabetical_68ROIs.txt'
    else:
        names_file = 'OrderingAlphabetical_86ROIs.txt'

    root_path = get_root_path()
    file_path = os.path.join(root_path,'atlases/')
    file_path = os.path.join(file_path, names_file)
    names_array = np.genfromtxt(file_path,dtype='str')
    return(names_array)
