"""
Shared utility functions for the model cluster correction notebooks
"""
import numpy as np
import pandas as pd
import h5py
from scipy.ndimage import label

def load_ieeg_data(filepath):
    """
    Load iEEG data from MATLAB v7.3 .mat file (HDF5 format).
    
    Parameters
    ----------
    filepath : str or Path
        Path to .mat file
        
    Returns
    -------
    dict
        Dictionary containing:
        - sid: subject IDs (array)
        - ch: channel labels (array)
        - hit_miss: hit (1) or miss (0) labels (array)
        - region: brain region labels (array) - only if present in file
        - data: timeseries data (channels × timepoints)
    """
    with h5py.File(str(filepath), 'r') as f:
        # load data arrays
        data = np.array(f['data']).T
        hit_miss = np.array(f['hit_miss']).flatten()
        
        # handle string arrays (stored as references in HDF5)
        sid = []
        for i in range(f['sid'].shape[1]):
            ref = f['sid'][0, i]
            sid.append(''.join(chr(int(c)) for c in f[ref][:].flatten()))
        sid = np.array(sid)
        
        ch = []
        for i in range(f['ch'].shape[1]):
            ref = f['ch'][0, i]
            ch.append(''.join(chr(int(c)) for c in f[ref][:].flatten()))
        ch = np.array(ch)
        
        # put in dictionary
        data_dict = {
            'sid': sid,
            'ch': ch,
            'hit_miss': hit_miss,
            'data': data
        }
        
        # include region only if it's actual data (not internal refs)
        try:
            if 'region' in f.keys() and hasattr(f['region'], 'shape'):
                region = []
                for i in range(f['region'].shape[1]):
                    ref = f['region'][0, i]
                    region.append(''.join(chr(int(c)) for c in f[ref][:].flatten()))
                region = np.array(region)
                data_dict['region'] = region
        except:
            pass
    
    return data_dict

def create_df(data_dict: dict) -> pd.DataFrame:
    """
    Create a structured DataFrame from the loaded data.
    
    Parameters
    ----------
    data_dict : dict
      Dictionary from load_ieeg_data()
    
    Returns
    -------
    pd.DataFrame
      Wide-format DataFrame with:
      - Rows: channels (one per channel per condition)
      - Columns: sid, channel, hit_miss, region (if present), + timepoint columns (t0, t1, ..., t138)
    """    
    n_ch, n_time = data_dict['data'].shape
    
    # create base df with metadata
    df = pd.DataFrame({
        'sid': data_dict['sid'],
        'channel': data_dict['ch'],
        'hit_miss': data_dict['hit_miss'].astype(int),
    })

    # add region if present
    if 'region' in data_dict:
        df['region'] = data_dict['region']
    
    # add timepoint columns
    time_cols = [f't{i}' for i in range(n_time)]
    time_data = pd.DataFrame(
        data_dict['data'],
        columns=time_cols
    )
    
    # concatenate metadata and timeseries data
    df = pd.concat([df, time_data], axis = 1)
    
    # convert categorical variables
    df['sid'] = df['sid'].astype('category')
    df['hit_miss'] = df['hit_miss'].astype('category')
    if 'region' in df.columns:
        df['region'] = df['region'].astype('category')
    
    return df

def cluster_test(datobs, datrnd, tail = 0, alpha = 0.05, clusteralpha = 0.05, clusterstat = 'sum'):
    """
    Perform cluster-based permutation test.
    
    Parameters
    ----------
    datobs : array
        Observed data (MxNx...xZ)
    datrnd : array
        Null distribution (MxNx...xZxPerm)
    tail : int, optional
        Test direction: -1 (datobs < null), 1 (datobs > null), 0 (two-tailed). Default is 0
    alpha : float, optional
        Critical level. Default is 0.05
    clusteralpha : float, optional
        Nonparametric threshold for cluster candidates. Default is 0.05
    clusterstat : str, optional
        How to combine statistics in clusters ('sum' or 'size'). Default is 'sum'
    
    Returns
    -------
    h : array
        Logical matrix indicating where significant clusters were found
    p : array
        Matrix of p-values associated with clusters
    clusterinfo : dict
        Extra cluster info
    """    
    # which dimension contains the randomizations
    rndsiz = datrnd.shape
    rnddim = len(rndsiz) - 1
    numrnd = rndsiz[rnddim]
    
    # handle cluster stat option
    cluster_stat_sum = (clusterstat == 'sum')
    cluster_stat_size = (clusterstat == 'size')
    
    # determine thresholds for cluster candidates
    if tail == 0:
        clusteralpha = clusteralpha / 2
    
    cluster_threshold_neg = np.quantile(datrnd, clusteralpha, axis = rnddim)
    cluster_threshold_pos = np.quantile(datrnd, 1 - clusteralpha, axis = rnddim)
    
    def compute_cluster_stats(dat, clus_cand):
        """
        Compute cluster statistics for binary mask.
        """
        labeled_array, num_clusters = label(clus_cand)
        clus_stats = []
        inds = []
        
        for l in range(1, num_clusters + 1):
            idx = np.where(labeled_array == l)
            
            if cluster_stat_sum:
                # ignore single-pixel clusters
                if len(idx[0]) > 1:
                    clus_stats.append(np.sum(dat[idx]))
                    inds.append(idx)
                elif cluster_stat_size:
                    clus_stats.append(len(idx[0]))
                    inds.append(idx)
        
        return np.array(clus_stats), inds
    
    def find_and_characterize_clusters(dat):
        """
        Find and characterize clusters in data.
        """
        clus_stats_pos, pos_inds = [], []
        clus_stats_neg, neg_inds = [], []
        
        if tail >= 0:
            clus_stats_pos, pos_inds = compute_cluster_stats(dat, dat >= cluster_threshold_pos)
        if tail <= 0:
            clus_stats_neg, neg_inds = compute_cluster_stats(dat, dat <= cluster_threshold_neg)
        
        if tail == -1:
            clus_stats_pos, pos_inds = np.array([]), []
        elif tail == 1:
            clus_stats_neg, neg_inds = np.array([]), []
        
        return clus_stats_pos, clus_stats_neg, pos_inds, neg_inds
    
    # cluster candidates for observed data
    clus_observed_pos, clus_observed_neg, pos_inds, neg_inds = find_and_characterize_clusters(datobs)
    
    # maximum and minimum cluster statistics for random data
    null_pos = np.full(numrnd, np.nan)
    null_neg = np.full(numrnd, np.nan)
    
    for k in range(numrnd):
        if k % round(numrnd / 10) == 0:
            print(f'Processing permutation {k+1} of {numrnd}...')
    
        clus_rnd_pos, clus_rnd_neg, _, _ = find_and_characterize_clusters(datrnd[..., k])
    
        if len(clus_rnd_pos) > 0:
            null_pos[k] = np.max(clus_rnd_pos)
        if len(clus_rnd_neg) > 0:
            null_neg[k] = np.min(clus_rnd_neg)
    
    null_pos = null_pos[~np.isnan(null_pos)]
    null_neg = null_neg[~np.isnan(null_neg)]

    # sort null distributions
    null_pos = np.sort(null_pos)[::-1]  # descending
    null_neg = np.sort(null_neg)  # ascending
    
    # compare observed clusters to max/min random clusters to obtain p-value
    clus_p_pos = np.ones(len(clus_observed_pos))
    for k in range(len(clus_observed_pos)):
        clus_p_pos[k] = (np.sum(null_pos > clus_observed_pos[k]) + 1) / (numrnd + 1)
    
    clus_p_neg = np.ones(len(clus_observed_neg))
    for k in range(len(clus_observed_neg)):
        clus_p_neg[k] = (np.sum(null_neg < clus_observed_neg[k]) + 1) / (numrnd + 1)
    
    # post-processing of output
    clusterinfo = {}
    
    # matrix of p-values
    p = np.ones(datobs.shape)
    
    if tail >= 0:
        for k in range(len(clus_p_pos)):
            p[pos_inds[k]] = clus_p_pos[k]
    
    if tail <= 0:
        for k in range(len(clus_p_neg)):
            if clus_p_neg[k] < p[neg_inds[k]][0]:
                p[neg_inds[k]] = clus_p_neg[k]
    
    if tail == 0:
        # for two-tailed, multiply p-values by 2
        p = np.minimum(1, p * 2)
    
    # result of the hypothesis test
    h = p < alpha
    
    return h, p, clusterinfo
    