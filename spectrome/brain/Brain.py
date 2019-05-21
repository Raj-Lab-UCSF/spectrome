import os, sys
sys.path.append("..")
import numpy as np
from utils import path as pth
from forward import get_complex_laplacian as fwd
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
        self.ordering = None
        self.laplacian = None
        self.eigenvalues = None
        self.norm_eigenmodes = None

        self.ntf_params = {'tau_e':0.012,
                           'tau_i':0.003,
                           'alpha':1.0,
                           'speed':5.0,
                           'gei':4.0,
                           'gii':1.0,
                           'tauC':0.006
                           }


    def add_ordering(self, filename):
        standard_list = pth.read_hdf5(filename)
        self.ordering = standard_list


    def add_functional_data(self, filename):
        '''Importing functional data for this brain'''
        self.timeseries_data = dr.read_dict(filename)


    def order_functional_data(self, orderfile):
        '''Reordering the functional data dictionary to match the standard given by the
        list in orderfile, for instance the HCP_list.h5 file in 'dictionaries'.'''
        self.timeseries_data = perm.order_dict(self.timeseries_data, orderfile)


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

    def add_laplacian_eigenmodes(self, connectome, distancematrix, w, speed = 10, num_ev = 10):
        "add complex Laplacian `L` and selected eigen modes and eigen values"
        L, selected_Evec, sorted_Eval = fwd.get_complex_laplacian(C = connectome, D = distancematrix, w = w, speed = speed, num_ev = num_ev)
        # Normalize eigen vectors for better visualization
        norm_eigs = np.zeros(selected_Evec.shape)
        for i in np.arange(0,num_ev):
            vdata = np.maximum(selected_Evec[:,i], 
                   np.mean(selected_Evec[:,i])-np.std(selected_Evec[:,i]))
            vdata = vdata - np.amin(vdata)
        
            vdata = np.minimum(vdata, np.mean(vdata)+np.std(vdata))
            vdata = vdata/np.amax(vdata)
            norm_eigs[:,i] = vdata

        self.laplacian = L
        self.raw_eigenvectors = selected_Evec
        self.norm_eigenmodes = norm_eigs
        self.eigenvalues = sorted_Eval

    def bi_symmetric_c(self):
        """Short summary.

        Args:
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
            max_dir (type): Description of parameter `max_dir`.
            f (type): Description of parameter `f`.

        Returns:
            type: Description of returned object.

        """
        thr = f*np.mean(self.connectome[self.connectome > 0])
        C = np.minimum(self.connectome, thr)
        C = max_dir * C + (1-max_dir) * C
        self.reducedConnectome = C
