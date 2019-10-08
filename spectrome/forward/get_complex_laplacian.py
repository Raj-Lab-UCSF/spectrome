"""Extract Laplacian from structural connectome based on 
transmission speed and oscillatory frequency"""

import numpy as np


def get_complex_laplacian(C, D, w, alpha=1, speed=10, num_ev=86):
    """ Extract complex laplacian based on frequency "omega", returns number of eigen
        vectors. This function sorts eigen vectors by ascending order, meaning the sorting
        begins with the eigen vectors associated with smallest absolute value of eigen values. 
        
        Args:
        - C (array): connectome
        - D (array): Distance matrix for C
        - omega (float): frequency  
        - alpha (float): coupling/diffusivity parameter
        - speed (int): default 10 m/s, transmission velocity
        - num_ev (int): number of eigen vectors you want as output
        Output:
        - L: complex laplacian
        - selected_Evec: eigen vectors
        - sorted_Evals: eigen values
    """
    nroi = C.shape[0]  # Get number of ROIs in connectome

    # Get row degree and col degree, make extremely big connections inf
    rowdegree = np.transpose(np.sum(C, axis=1))  # row degrees
    coldegree = np.sum(C, axis=0)  # column degrees
    qind = rowdegree + coldegree < 0.2 * np.mean(rowdegree + coldegree)
    rowdegree[qind] = np.inf
    coldegree[qind] = np.inf

    Tau = 0.001 * D / speed  # delay as a function of distance and transmission speed
    Cc = C * np.exp(-1j * Tau * w)  # Complex Laplacian

    # Compute Laplacian and decompose into eigen vectors
    # L1 = 0.8*np.identity(nroi)
    L1 = np.identity(nroi)
    L2 = np.divide(
        1, np.sqrt(np.multiply(rowdegree, coldegree)) + np.spacing(1)
    )  # diag(1./(sqrt(rowdegree.*coldegree)+eps));
    L = L1 - alpha * np.matmul(
        np.diag(L2), Cc
    )  # Final Laplacian = I - alpha*1/diag(sqrt(deg)

    # decomposition with linalg.eig
    K = nroi
    d, v = np.linalg.eig(L)
    # sorting in ascending & absolute value
    eig_ind = np.argsort(np.abs(d))
    eig_vec = v[:, eig_ind]
    # abseiv = np.abs(eig_vec[:,0])
    sorted_Evals = d[eig_ind]

    # defining some intermediate variables
    ev = np.transpose(sorted_Evals[0:K])
    Vv = eig_vec[:, 0:K]
    Vvec = np.asarray(Vv)

    # Select eigen modes based on num_ev as output...
    selected_Evec = []
    for k in np.arange(0, num_ev):
        abs_Vvec = np.abs(Vvec[:, k])
        selected_Evec.append(abs_Vvec)  # this appends as rows

    selected_Evec = np.transpose(
        np.asarray(selected_Evec)
    )  # transpose back to column vectors
    return L, selected_Evec, sorted_Evals
