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


def exp_neg_dist_Cij(distance_matrix):
    """Generate network weighted by exponent of negative distance:
    
    Args:
        distance_matrix (array): input distance matrix
    
    Returns:
        dist_wCij [array]: symmetric distance weighted network with 0 diagonal
    """
    negdist = np.exp(-distance_matrix)
    negdist = negdist - np.diag(np.diag(negdist))

    V = len(negdist)  # number of vertices
    E = np.floor(np.count_nonzero(negdist) / 2).astype(int)  # number of edges

    linear_indices = np.squeeze(
        np.asarray(np.where((negdist < 1) & (negdist > 0.0001)))
    )  # get rid of extremes
    dist_distribution = negdist[linear_indices[0, :], linear_indices[1, :]]

    # Sample from dist_distribution into a matrix:
    triu_indices = np.triu(np.logical_not(np.eye(V)))
    i = np.asarray(
        triu_indices.ravel().nonzero()
    )  # linear indices of upper triangle elements

    RP = np.random.permutation(i.shape[1])  # randomly perm the indices
    random_indices = i[:, RP]  # vertices randomly permuted
    random_indices = np.squeeze(
        np.asarray(np.unravel_index(random_indices, negdist.shape))
    )

    dist_wCij = np.zeros([V, V])  # V-vertices matrix
    for n in np.arange(0, E):
        dist_wCij[random_indices[0, n], random_indices[1, n]] = np.random.choice(
            dist_distribution, replace=False
        )

    dist_wCij = dist_wCij + np.transpose(dist_wCij)
    return dist_wCij


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
