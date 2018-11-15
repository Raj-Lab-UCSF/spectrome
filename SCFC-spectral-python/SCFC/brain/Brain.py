import numpy as np

class Brain:
    """Short summary.

    Attributes:
        atlas (type): Description of parameter `atlas`.
        connectome (type): Description of parameter `connectome`.
        ntf_parameters (type): Description of parameter `ntf_parameters`.

    """
    def __init__(self):
        #Body variables
        self.atlas = None
        self.connectome = None
        self.ordering = None
        self.ntf_parameters = {'tau_e':0.012,
                               'tau_i':0.003,
                               'alpha':1.0,
                               'speed':5.0,
                               'gei':4.0,
                               'gii':1.0,
                               'tauC':0.006
                               }

def get_Julia_order():
    """Get Julia Owen's brain region ordering (specific for DK86 atlas).

    Args:

    Returns:
        permJulia (type): Brain region orders for all regions
        emptyJulia (type): Brain regions with no MEG
        cortJulia (type): Brain cortical regions.

    """
    cortJulia_lh = np.array([0, 1, 2, 3, 4, 6, 7, 8, 10, 11, 12, 13, 14,
                             15, 17, 16, 18, 19, 20, 21, 22, 23, 24, 25,
                             26, 27, 28, 29, 30, 31, 5, 32, 33, 9])
    qsubcort_lh = np.array([0, 40, 36, 39, 38, 37, 35, 34, 0])
    qsubcort_rh = qsubcort_lh + 34 + 1
    cortJulia = np.concatenate([cortJulia_lh, 34 + cortJulia_lh])
    cortJulia_rh = cortJulia_lh + 34 + 7
    permJulia = np.concatenate([cortJulia_lh, cortJulia_rh,
                                qsubcort_lh, qsubcort_rh])
    emptyJulia = np.array([68, 77, 76, 85])
    return permJulia, emptyJulia, cortJulia


def get_HCP_connectome(hcp_dir,
                       conmat_in='mean80_fibercount.csv',
                       dmat_in='mean80_fiberlength.csv'):
    """Short summary.

    Args:
        hcp_dir (str): directory to HCP connectome.
        conmat_in (type): name of connectivity csv file.
        dmat_in (type): name of fiber distance csv file.

    Returns:
        Cdk_conn(arr): Connectivity matrix oredered according to permHCP
        Ddk_conn(arr): Distance matrix ordered by permHCP
        permHCP(arr): Ordering of brain regions in DK86 atlas

    """
    cdk_hcp = np.genfromtxt(os.path.join(hcp_dir, conmat_in),
                            delimiter=',', skip_header=1)

    ddk_hcp = np.genfromtxt(os.path.join(hcp_dir, dmat_in),
                            delimiter=',', skip_header=0)

    permHCP = np.concatenate([np.arange(18, 52),
                              np.arange(52, 86),
                              np.arange(0, 9),
                              np.arange(9, 18)])

    Cdk_conn = cdk_hcp[permHCP, ][:, permHCP]
    Ddk_conn = ddk_hcp[permHCP, ][:, permHCP]

    return Cdk_conn, Ddk_conn, permHCP

def bi_symmetric_c(Cdk_conn, linds, rinds):
    """Short summary.

    Args:
        Cdk_conn (type): Description of parameter `Cdk_conn`.
        linds (type): Description of parameter `linds`.
        rinds (type): Description of parameter `rinds`.

    Returns:
        type: Description of returned object.

    """
    q = np.maximum(Cdk_conn[linds, :][:, linds], Cdk_conn[rinds, :][:, rinds])
    q1 = np.maximum(Cdk_conn[linds, :][:, rinds], Cdk_conn[rinds, :][:, linds])
    Cdk_conn[np.ix_(linds, linds)] = q
    Cdk_conn[np.ix_(rinds, rinds)] = q
    Cdk_conn[np.ix_(linds, rinds)] = q1
    Cdk_conn[np.ix_(rinds, linds)] = q1
    return Cdk_conn


def reduce_extreme_dir(Cdk_conn, max_dir=0.95, f=7):
    """Short summary.

    Args:
        Cdk_conn (type): Description of parameter `Cdk_conn`.
        max_dir (type): Description of parameter `max_dir`.
        f (type): Description of parameter `f`.

    Returns:
        type: Description of returned object.

    """
    thr = f*np.mean(Cdk_conn[Cdk_conn > 0])
    C = np.minimum(Cdk_conn, thr)
    C = max_dir * C + (1-max_dir) * C
    return C
