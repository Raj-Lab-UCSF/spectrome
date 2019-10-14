""" functions to sort and compute stats on eigen modes"""
import pandas as pd
import numpy as np

from scipy.spatial import distance
from scipy.stats import entropy, spearmanr, pearsonr

def eig_fc_get_standardz(x, y, binary_thresh=0.1, nreps=1000):
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
    ub, lb = 1, 0  # 1 or 0s after thresholding
    x = np.where(x > binary_thresh, ub, lb)
    y = np.where(y > binary_thresh, ub, lb)

    np.random.seed(24)

    mean_std_perm = np.zeros(nreps)
    # std_perm = np.zeros(1000)
    for i in np.arange(nreps):
        xperm = np.random.permutation(x)
        yperm = np.random.permutation(y)
        mean_std_perm[i] = get_sxy(xperm, yperm)

    sxy = get_sxy(x, y)
    # compute mean and standard deviation:
    zxy = (sxy - np.mean(mean_std_perm)) / np.std(mean_std_perm)
    return zxy, sxy


def get_overlap_score_dfs(all_eig, all_fc_networks, threshold=0.1):
    """[summary]
    
    Args:
        all_eig (array): eigen modes that you want to compare
        all_fc_networks (array): all canonical networks being compared to
    
    Returns:
        df_overlap_score [Pandas DataFrame]: DataFrame with all overlap scores
        df_sxy [Pandas DataFrame]: DataFrame with 
    """

    df_cols = all_fc_networks.index
    df_ind = ["Eig #%d" % x for x in np.arange(all_eig.shape[1]) + 1]
    df_overlap_score = pd.DataFrame([], index=df_ind, columns=df_cols)
    df_sxy = pd.DataFrame([], index=df_ind, columns=df_cols)
    eigcounter = 0
    for eignum in df_sxy.index:
        for name in all_fc_networks.index:
            fc_vector = all_fc_networks.loc[name].values
            eig_vector = all_eig[:, eigcounter]
            df_overlap_score.at[eignum, name], df_sxy.at[
                eignum, name
            ] = eig_fc_get_standardz(eig_vector, fc_vector, binary_thresh=threshold)
        eigcounter += 1
    return df_overlap_score, df_sxy


def get_sxy(x, y):
    """Get s(x,y) where
    Inputs:
        x - module/eigen mode
        y - canoical network 
    """
    ind_eig = np.where(x > 0.1)
    ind_fc = np.where(y > 0.1)
    # find intersection and union
    fc_eig_intersect = np.intersect1d(ind_eig, ind_fc)
    fc_eig_union = np.union1d(ind_eig, ind_fc)
    # s(x,y):
    overlap_sxy = np.sum(x[fc_eig_intersect]) / np.sum(x[fc_eig_union])
    # normalize by number of regions thats in the canonical network:

    return overlap_sxy


def get_dice_df(x, y):
    """[Dice similarity score betewen two boolean arrays, then make it into a Pandas Dataframe
    translate into dataframes after we compute all the arrays?]
    
    Args:
        x ([array]): [eigen mode array with boolean array values]
        y ([Dataframe]): [canonical network data frame with boolean array vaues]

    Returns:
        dice_df [array]: dice scores between each eigen mode and each canonical network    
    """
    df_cols = y.index
    df_ind = ["Eig #%d" % x for x in np.arange(x.shape[1]) + 1]
    df_dice = pd.DataFrame([], index=df_ind, columns=df_cols)
    eigcounter = 0
    for eignum in df_dice.index:
        for name in y.index:
            em = x[:, eigcounter]  # eigen mode values
            fc = y.loc[name].values
            # df_dice.at[eignum, name] = (distance.dice(em,fc)*np.count_nonzero(fc))
            df_dice.at[eignum, name] = distance.dice(
                em, fc
            )  # without multiplying by support
        eigcounter += 1
    return df_dice

def get_corr_df(x,y, method = 'spearman'):
    """Pearson's or Spearman's correlation between arrays
    """
    df_cols = y.index
    df_ind = ["Eig #%d" % x for x in np.arange(x.shape[1]) + 1]
    df_corr = pd.DataFrame([], index=df_ind, columns=df_cols)
    eigcounter = 0
    for eignum in df_corr.index:
        for name in y.index:
            em = x[:,eigcounter] # eigenmode
            fc = np.nan_to_num(y.loc[name].values)
            if method == 'pearson':
                df_corr.at[eignum, name] = pearsonr(em, fc)[0]
            elif method == 'spearman':
                df_corr.at[eignum, name] = spearmanr(em,fc)[0]
        eigcounter += 1
    return df_corr

def get_coupled_corr()

"""
def get_jaccard_df(x, y):
    #[Jaccard similarity score betewen two boolean arrays, then make it into a Pandas Dataframe
    #translate into dataframes after we compute all the arrays?]
    
    #Args:
    #    x ([array]): [eigen mode array with boolean array values]
    #    y ([Dataframe]): [canonical network data frame with boolean array vaues]

    #Returns:
    #    jaccard_df [array]: dice scores between each eigen mode and each canonical network    

    df_cols = y.index
    df_ind = ["Eig #%d" % x for x in np.arange(x.shape[1]) + 1]
    jaccard_df = pd.DataFrame([], index=df_ind, columns=df_cols)
    eigcounter = 0
    for eignum in jaccard_df.index:
        for name in y.index:
            em = x[:, eigcounter]  # eigen mode values
            fc = y.loc[name].values
            jaccard_df.at[eignum, name] = jaccard_score(
                em, fc
            )  # without multiplying by support
        eigcounter += 1
    return jaccard_df
"""

def get_entropy_score(df_sxy):
    """[Compute purity score as a measure of entropy]
    
    Args:
        df_sxy (Pandas DataFrame): DataFrame of either overlap scores or standardized overlap scores
    
    
    Returns:
        purity_score [array]: [purity score for each eigenmode]
    """

    # df_sxy = df_sxy.div(df_sxy.sum(axis=1), axis = 0)
    overlap_scores = df_sxy.values
    # ind_zeros = overlap_scores == 0
    # overlap_scores[ind_zeros] = np.spacing(1)

    entropy_vec = []
    for row in overlap_scores:
        # purity_vec.append(np.multiply(row, np.log2(row.astype(np.float64))))
        entropy_vec.append(entropy(row.astype("float")))
        # purity_vec.append(np.multiply(row, np.log2(np.abs(row.astype(np.float64))+1)))

    # purity_score = -np.sum(purity_vec, axis = 1)
    return np.asarray(entropy_vec)
