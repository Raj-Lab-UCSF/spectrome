'''Functions to read in the data and return required inputs to the other analysis functions,
including reordering of the data and matrices according to standard dictionary orders.'''
import sys, os
sys.path.append("..")
from utils import path as pth

def read_dict(filename):
    """read_MEG_dict. Reads a HDF5 file containing a dictionary of MEG data, keyed
    by brain region. Convention for the cortex is 'R/LH[lowercasebrainregion]' where
    spaces are not permitted.

    Args:
        filename (type): filename.

    Returns:
        type: Returns dictionary, keyed by brain region.

    """
    data = pth.read_hdf5(filename)
    return data


def get_MEG_data_frommat(sub_name, ordering, MEGfolder):
    """Get source localized MEG data and arrange it following ordering method
    --using the matrix file (older function, deprecated).

    Args:
        sub_name (str): Name of subject.
        ordering (arr): Cortical region ordering (e.g. `cortJulia`).
        MEGfolder (str): Directory for MEG data.

    Returns:
        MEGdata (arr): MEG data.
        coords (type):

    """
    S = loadmat(os.path.join(MEGfolder, sub_name, 'DK_timecourse_20.mat'))
    MEGdata = S['DK_timecourse']
    MEGdata = MEGdata[ordering, ]
    C = loadmat(os.path.join(MEGfolder, sub_name, 'DK_coords_meg.mat'))
    coords = C['DK_coords_meg']
    coords = coords[ordering, ]
    del S, C

    return MEGdata, coords
