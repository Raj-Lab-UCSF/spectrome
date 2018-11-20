import os, sys
sys.path.append("..")
import numpy as np
from utils import path as pth
import read.data_reader as dr
import preprocess.permute as perm

class Brain:
    """A class containing data that represents a single brain.

    Attributes:
        connectome (array): Array of the connectome.
        reducedConnectome (array): Connectome with extreme components culled (?? Add doc please).
        distance_matrix (array): Matrix of distances between brain regions.
        permutation (array): The permutation applied to the data file connectome.
        ordering (type): Description of parameter `ordering`.
        ntf_params (dict): Parameters for the network transfer model.

    """

    def __init__(self):
        #Body variables
        self.connectome = None
        self.reducedConnectome = None
        self.distance_matrix = None
        self.permutation  = None
        self.ordering = {'permJulia':None,
                         'emptyJulia':None,
                         'cortJulia':None}

        self.ntf_params = {'tau_e':0.012,
                           'tau_i':0.003,
                           'alpha':1.0,
                           'speed':5.0,
                           'gei':4.0,
                           'gii':1.0,
                           'tauC':0.006
                           }

    def add_MEG(self, filename):
        '''Importing MEG data for this brain'''
        self.MEGdata = dr.read_dict(filename)

    def order_MEG(self, orderfile):
        '''Reordering the MEG data dictionary to match the standard given by the
        list in orderfile, for instance the HCP_list.h5 file in 'dictionaries'.'''
        self.MEGdata = perm.order_dict(self.MEGdata, orderfile)


    def add_connectome(self, hcp_dir,
                           conmat_in='mean80_fibercount.csv',
                           dmat_in='mean80_fiberlength.csv'):

        self.connectome = np.genfromtxt(os.path.join(hcp_dir, conmat_in),
                                delimiter=',', skip_header=1)

        self.distance_matrix = np.genfromtxt(os.path.join(hcp_dir, dmat_in),
                                delimiter=',', skip_header=0)

    def add_ordered_connectome(self, confile, distfile):
        '''add a connectome and distance matrix using ordering directly'''
        con, dist, permutation = perm.reorder_connectome(conmatfile = confile, distmatfile = distfile)
        self.connectome = con
        self.distance_matrix = dist
        self.permutation = permutation

    def reorder_connectome(self, connectome, distancematrix):
        '''re-order the present connectome and distance matrix -- note that this is
        a first iteration and some work needs to be done to make it flexible with regards
        the specific ordering.'''
        con, dist, permutation = perm.reorder_connectome(conmat = connectome, distmat = distancematrix)
        self.connectome = con
        self.distance_matrix = dist
        self.permutation = permutation

#########THe functions below this all become a bit mysterious! More explanation please!


    def set_julia_order(self):
        """Set Julia Owen's brain region ordering (specific for DK86 atlas).

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
        cortJulia_rh = cortJulia_lh + 34 + 7


        self.ordering['emptyJulia'] = np.array([68, 77, 76, 85])
        self.ordering['cortJulia'] = np.concatenate([cortJulia_lh, 34 + cortJulia_lh])
        self.ordering['permJulia'] = np.concatenate([cortJulia_lh, cortJulia_rh,
                                    qsubcort_lh, qsubcort_rh])


    def bi_symmetric_c(self):
        """Short summary.

        Args:
            Cdk_conn (type): Description of parameter `Cdk_conn`.
            linds (type): Description of parameter `linds`.
            rinds (type): Description of parameter `rinds`.

        Returns:
            type: Description of returned object.

        """
        # Some other ordering that was in the original code:
        linds = np.concatenate([np.arange(0,34), np.arange(68,77)])
        rinds = np.concatenate([np.arange(34,68), np.arange(77,86)])

        q = np.maximum(self.connectome[linds, :][:, linds], self.connectome[rinds, :][:, rinds])
        q1 = np.maximum(self.connectome[linds, :][:, rinds], self.connectome[rinds, :][:, linds])
        self.connectome[np.ix_(linds, linds)] = q
        self.connectome[np.ix_(rinds, rinds)] = q
        self.connectome[np.ix_(linds, rinds)] = q1
        self.connectome[np.ix_(rinds, linds)] = q1


    def reduce_extreme_dir(self, max_dir=0.95, f=7):
        """Short summary.

        Args:
            Cdk_conn (type): Description of parameter `Cdk_conn`.
            max_dir (type): Description of parameter `max_dir`.
            f (type): Description of parameter `f`.

        Returns:
            type: Description of returned object.

        """
        thr = f*np.mean(self.connectome[self.connectome > 0])
        C = np.minimum(self.connectome, thr)
        C = max_dir * C + (1-max_dir) * C
        self.reducedConnectome = C
