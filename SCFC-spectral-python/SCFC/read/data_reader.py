'''Functions to read in the data and return required inputs to the other analysis functions'''

def sorting(data, standard):
    '''If your data is in some stupid order, this should reorder it according to the
    standard you are using -- assuming you've passed in data as a dicitonary'''


def get_MEG_data(sub_name, ordering, MEGfolder):
    """Get source localized MEG data and arrange it following ordering method.

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
