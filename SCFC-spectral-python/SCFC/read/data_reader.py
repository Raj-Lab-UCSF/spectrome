'''Functions to read in the data and return required inputs to the other analysis functions'''
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
    print(standardlist)

    for key in standardlist:
        if key in dataorder:
            newdata[key] = data[key]
            loc_in_standard.append(standardlist.index(key))
        else:
            print('Skipping region of brain -- not in data')

    return newdata


def get_MEG_data_frommat(sub_name, ordering, MEGfolder):
    """Get source localized MEG data and arrange it following ordering method
    --using the matrix file (older function).

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

if __name__ == "__main__":
    data = read_MEG_dict('/Users/Megan/RajLab/MEG-chang/8002.101/DK_timecourse_20.h5')
    # print(data['LHbankssts'].shape)
    HCP_path = pth.get_sibling_path('dictionaries')
    HCPfilename = os.path.join(HCP_path, 'HCP_order.h5')
    order_MEG(data, HCPfilename)
