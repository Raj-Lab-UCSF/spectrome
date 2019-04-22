''' functions to sort and compute stats on eigen modes'''
import pandas as pd
import numpy as np

def eig_fc_get_standardz(x,y, binary_thresh = 0, binarize = False):
    """Permutes both canonical networks and input eigenmodes 1000 times
    and calculates the standardized overlap score between each eigen mode
    and each canonical network, standardize by mean and standard deviation

    
    Args:
        x (array): [eigenmode vector]
        y (array): [canonical functional network vector]
        binary_thresh (int, optional): Defaults to 0. threshold for binarizing regional values
        binarize (bool, optional): Defaults to False. Do you want to binarize or not
    
    Returns:
        zxy [array]: standardized overlapscores
        sxy [array]: non-standardized overlap scores
    """

    np.random.seed(24)
    if binarize is True:
        eig_thresh_ind = x > binary_thresh
        fc_thresh_ind = y > binary_thresh
        x[eig_thresh_ind] = 1
        y[fc_thresh_ind] = 1
        
    mean_std_perm = np.zeros(1000)
    #std_perm = np.zeros(1000)
    for i in np.arange(1000):
        xperm = np.random.permutation(x)
        yperm = np.random.permutation(y)
        mean_std_perm[i] = get_sxy(xperm,yperm)
        
    sxy = get_sxy(x,y)
    # compute mean and standard deviation:
    zxy = (sxy - np.mean(mean_std_perm))/np.std(mean_std_perm)
    return zxy, sxy

def get_overlap_score_dfs(all_eig, all_fc_networks, binarize = False):
    """[summary]
    
    Args:
        all_eig (array): eigen modes that you want to compare
        all_fc_networks (array): all canonical networks being compared to
    
    Returns:
        df_overlap_score [Pandas DataFrame]: DataFrame with all overlap scores
        df_sxy [Pandas DataFrame]: DataFrame with 
    """

    df_cols = all_fc_networks.index
    df_ind = ['Eig #%d' % x for x in np.arange(all_eig.shape[0])+1]
    df_overlap_score = pd.DataFrame([], index = df_ind, columns = df_cols)
    df_sxy = pd.DataFrame([], index = df_ind, columns = df_cols)
    eigcounter = 0
    for eignum in df_sxy.index:
        for name in all_fc_networks.index:
            fc_vector = all_fc_networks.loc[name].values
            eig_vector = all_eig[eigcounter,:]
            df_overlap_score.at[eignum, name], df_sxy.at[eignum, name] = eig_fc_get_standardz(eig_vector, fc_vector, binarize = binarize)
        eigcounter += 1
    return df_overlap_score, df_sxy

def get_sxy(x,y):
    """Get s(x,y) where
    Inputs:
        x - module/eigen mode
        y - canoical network 
    """
    ind_eig = np.where(x > 0)
    ind_fc = np.where(y > 0)
    # find intersection and union
    fc_eig_intersect = np.intersect1d(ind_eig, ind_fc)
    fc_eig_union = np.union1d(ind_eig, ind_fc)
    #s(x,y):
    overlap_sxy = np.sum(x[fc_eig_intersect])/np.sum(x[fc_eig_union])
    return overlap_sxy

def get_purity_score(df_sxy):
    """[summary]
    
    Args:
        df_sxy (Pandas DataFrame): DataFrame of either overlap scores or standardized overlap scores
    
    
    Returns:
        purity_score [array]: [purity score for each eigenmode]
    """

    #df_sxy = df_sxy.div(df_sxy.sum(axis=1), axis = 0)
    overlap_scores = df_sxy.values
    #ind_zeros = overlap_scores == 0
    #overlap_scores[ind_zeros] = np.spacing(1)
    
    purity_vec = []
    for row in overlap_scores:
        #purity_vec.append(np.multiply(row, np.log2(row.astype(np.float64))))
        purity_vec.append(np.multiply(row, np.log2(np.abs(row.astype(np.float64))+1)))
    
    purity_score = -np.sum(purity_vec, axis = 1)
    return np.asarray(purity_score), overlap_scores