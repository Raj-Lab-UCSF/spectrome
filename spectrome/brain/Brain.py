import numpy as np
import os
from sklearn.preprocessing import minmax_scale

from ..read import data_reader as dr
from ..preprocess import permute as perm
from ..utils import path as pth
from ..forward import laplacian as fwd


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
        # Body variables
        self.connectome = None
        self.reducedConnectome = None
        self.distance_matrix = None
        self.permutation = None
        self.ordering = None
        self.laplacian = None
        self.eigenvalues = None
        self.norm_eigenmodes = None
        self.regular_eigenvalues = None
        self.regular_laplacian = None
        self.norm_regular_eigenmodes = None
        self.raw_regular_eigenvectors = None

        self.ntf_params = {
            "tau_e": 0.012,
            "tau_i": 0.003,
            "alpha": 1.0,
            "speed": 5.0,
            "gei": 4.0,
            "gii": 1.0,
            "tauC": 0.006,
        }

    def add_ordering(self, filename):
        standard_list = pth.read_hdf5(filename)
        self.ordering = standard_list

    def add_functional_data(self, filename):
        """Importing functional data for this brain"""
        self.timeseries_data = dr.read_dict(filename)

    def order_functional_data(self, orderfile):
        """Reordering the functional data dictionary to match the standard given by the
        list in orderfile, for instance the HCP_list.h5 file in 'atlases'."""
        self.timeseries_data = perm.order_dict(self.timeseries_data, orderfile)

    def add_connectome(
        self,
        hcp_dir,
        conmat_in="mean80_fibercount.csv",
        dmat_in="mean80_fiberlength.csv",
    ):

        self.connectome = np.genfromtxt(
            os.path.join(hcp_dir, conmat_in), delimiter=",", skip_header=1
        )

        self.distance_matrix = np.genfromtxt(
            os.path.join(hcp_dir, dmat_in), delimiter=",", skip_header=0
        )

    def add_ordered_connectome(self, confile, distfile):
        """add a connectome and distance matrix using ordering directly"""
        con, dist, permutation = perm.reorder_connectome(
            conmatfile=confile, distmatfile=distfile
        )
        self.connectome = con
        self.distance_matrix = dist
        self.permutation = permutation

    def reorder_connectome(self, connectome, distancematrix):
        """re-order the present connectome and distance matrix -- note that this is
        a first iteration and some work needs to be done to make it flexible with regards
        the specific ordering."""
        con, dist, permutation = perm.reorder_connectome(
            conmat=connectome, distmat=distancematrix
        )
        self.connectome = con
        self.distance_matrix = dist
        self.permutation = permutation

    def decompose_complex_laplacian(self, alpha, k, f = None, speed = None, num_ev=86, vis = False, take_abs = True):
        "Add complex laplacian `L` and selected eigenmodes and eigen values based on 2 parameters alpha and k"
        L, selected_Evec, sorted_Eval = fwd.decompose_complex_laplacian(
            C=self.reducedConnectome,
            D=self.distance_matrix,
            alpha=alpha,
            k=k,
            f = f,
            speed = speed,
            num_ev=num_ev,
            take_abs=take_abs
        )

        # Normalize eigen vectors for better visualization
        if vis == True:
            norm_eigs = np.zeros(selected_Evec.shape)
            for i in np.arange(0, num_ev):
                vdata = np.maximum(
                    selected_Evec[:, i],
                    np.mean(selected_Evec[:, i]) - np.std(selected_Evec[:, i]),
                )
                vdata = vdata - np.amin(vdata)

                vdata = np.minimum(vdata, np.mean(vdata) + np.std(vdata))
                vdata = vdata / np.amax(vdata)
                norm_eigs[:, i] = vdata
        else:
            if take_abs is True:
                norm_eigs = minmax_scale(selected_Evec, axis=1)
            elif take_abs is False:
                norm_eigs = selected_Evec

        self.complex_laplacian = L
        self.raw_eigenvectors = selected_Evec
        self.norm_eigenmodes = norm_eigs
        self.eigenvalues = sorted_Eval

    def decompose_regular_laplacian(self, alpha=1, num_ev=86, vis=True):
        L, selected_Evec, sorted_Eval = fwd.decompose_regular_laplacian(
            C=self.reducedConnectome, alpha=alpha, num_ev=num_ev
        )
        # Normalize eigen vectors
        if vis == True:
            norm_eigs = np.zeros(selected_Evec.shape)
            for i in np.arange(0, num_ev):
                vdata = np.maximum(
                    selected_Evec[:, i],
                    np.mean(selected_Evec[:, i]) - np.std(selected_Evec[:, i]),
                )
                vdata = vdata - np.amin(vdata)

                vdata = np.minimum(vdata, np.mean(vdata) + np.std(vdata))
                vdata = vdata / np.amax(vdata)
                norm_eigs[:, i] = vdata
        else:
            norm_eigs = minmax_scale(selected_Evec, axis=1)
        self.regular_laplacian = L
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
        linds = np.concatenate([np.arange(0, 34), np.arange(68, 77)])
        rinds = np.concatenate([np.arange(34, 68), np.arange(77, 86)])

        q = np.maximum(
            self.connectome[linds, :][:, linds], self.connectome[rinds, :][:, rinds]
        )
        q1 = np.maximum(
            self.connectome[linds, :][:, rinds], self.connectome[rinds, :][:, linds]
        )
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
        thr = f * np.mean(self.connectome[self.connectome > 0])
        C = np.minimum(self.connectome, thr)
        C = max_dir * C + (1 - max_dir) * C
        self.reducedConnectome = C
