import numpy as np


def random_Cij_und(V, E):
    """Generate undirected binary random network:
    
    Args:
        V (int): number of vertices
        E (int): number of edges
    
    Returns:
        Cij [array]: undirected, binary random network with 0 diagonal
    """
    indmat = np.triu(np.logical_not(np.eye(V)))
    i = np.asarray(
        indmat.ravel().nonzero()
    )  # linear indices of upper triangle elements
    RP = np.random.permutation(i.shape[1])
    random_indices = i[:, RP]  # vertices randomly permuted

    # select E edges out of random_indices:
    Cij = np.zeros([V, V])
    V_indices = np.squeeze(np.asarray(np.unravel_index(random_indices, indmat.shape)))
    Cij[V_indices[0, 0:E], V_indices[1, 0:E]] = 1  # every selected Edge gets 1
    Cij = Cij + np.transpose(Cij)
    return Cij


def random_Cij_dir(V, E):
    """Generate directed binary random network:
    
    Args:
        V (int): number of vertices
        E (int): number of edges
    
    Returns:
        Cij [array]: directed, binary random network with 0 diagonal, unsymmetric
    """
    indmat = np.logical_not(np.eye(V))
    i = np.asarray(
        indmat.ravel().nonzero()
    )  # linear indices of upper triangle elements
    RP = np.random.permutation(i.shape[1])
    random_indices = i[:, RP]  # vertices randomly permuted

    # select E edges out of random_indices:
    Cij = np.zeros([V, V])
    V_indices = np.squeeze(np.asarray(np.unravel_index(random_indices, indmat.shape)))
    Cij[V_indices[0, 0:E], V_indices[1, 0:E]] = 1  # every selected Edge gets 1
    return Cij


def add_weights(Cij, u, s):
    """Add weights to a binary network based on mean and standard deviation
    
    Args:
        Cij (array): input matrix
        u ([type]): edge mean
        s ([type]): edge variance
    
    Returns:
       wCij (array): Cij with weights added to non-zero vertices
    """
    # find vertices in upper triangle:
    inds = np.squeeze(
        np.asarray(np.where(np.triu(Cij)))
    )  # find upper triangle vertices
    # create new weighted Cij
    wCij = np.zeros([86, 86])
    # assign value to inds based on mean u and variance s
    for i in np.arange(0, inds.shape[1]):
        wCij[inds[0, i], inds[1, i]] = u + np.random.random(1) * s

    wCij = wCij + np.transpose(wCij)
    return wCij


def exp_neg_dist_Cij(distance_matrix, sparsity = 0.8):
    """Generate network weighted by exponent of negative distance given sparsity:
    
    Args:
        distance_matrix (array): input distance matrix
    
    Returns:
        negdist [array]: symmetric distance weighted network with 0 diagonal
    """
    negdist = np.exp(-distance_matrix)
    negdist = negdist - np.diag(np.diag(negdist))

    V = len(negdist)  # number of vertices
    #compute current sparsity:
    current_sparsity = np.count_nonzero(negdist == 0)/(len(negdist)**2)

    # Remove lowest distances to achieve sparsity:
    while current_sparsity < sparsity:
        triu_negdist = np.triu(negdist)
        triu_nonzeros = np.asarray(np.triu(negdist).ravel().nonzero())
        # replace lowest distances with 0's
        min_loc = np.where(triu_negdist == np.amin(triu_negdist[np.unravel_index(triu_nonzeros, negdist.shape)]))
        triu_negdist[min_loc] = 0
        #update sparsity:
        negdist = triu_negdist + np.transpose(triu_negdist)
        current_sparsity = np.count_nonzero(negdist == 0)/(len(negdist)**2)

    return negdist


def distance_matrix(C, distances):
    """Generate distance matrix for a given matrix, the values of distances in output is sampled from a given distance matrix
    
    Args:
        C ([array]): [connectome you want to pair a random distance matrix for]
        distances ([array]): [description]
    
    Returns:
        [array]: [distance matrix]
    """
    V = len(C)  # number of vertices
    Drand = np.zeros([V, V])

    # find index of upper triangle where edges have non-zero value:
    triu_inds = np.triu(C)
    lin_i = np.asarray(triu_inds.nonzero())

    # create distance distribution to sample from
    dist_ind = np.asarray(np.nonzero(np.triu(distances)))
    dist_dist = distances[dist_ind[0, :], dist_ind[1, :]]

    # assign randomly sampled values to lin_i
    for i in np.arange(0, lin_i.shape[1]):
        Drand[lin_i[0, i], lin_i[1, i]] = np.random.choice(dist_dist, replace=False)

    Drand = Drand + np.transpose(Drand)
    return Drand
