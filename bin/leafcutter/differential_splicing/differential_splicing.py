import pandas as pd
import numpy as np
import torch
from collections import defaultdict, OrderedDict
from leafcutter.differential_splicing.dm_glm import dirichlet_multinomial_anova, SimpleGuide, CleverGuide
from leafcutter.differential_splicing import bayes_glm
from timeit import default_timer as timer
from tqdm import tqdm
from multiprocessing import Pool # 33s to fit. 
from functools import partial
import pyro.distributions as dist

#import concurrent.futures as cnc # not compatible with pyro
#from pathos.multiprocessing import ProcessingPool as Pool # 24s
#from joblib import Parallel, delayed # 60s

try: 
    from scipy.stats import false_discovery_control
except: 
    print("Warning: your scipy version is older than 1.11.0 so scipy.stats.false_discovery_control isn't available.")
    def false_discovery_control(ps, method = "bh"):
        order = np.argsort(ps) # put p values in ascending order
        ps = ps[order]
        m = len(ps)
        i = np.arange(1, m+1)
        ps *= m / i # adjust
        if method == 'by': ps *= np.sum(1. / i)
        ps = ps[np.argsort(order)] # put back in original order
        return np.clip(ps, 0, 1)
    
def robust_fdr(ps, method = "bh"): 
    mask = ~np.isnan(ps)
    p_adjust = ps.copy()
    p_adjust[mask] = false_discovery_control(ps[mask], method = method)
    return p_adjust
    
def task(inp, x, torch_types, kwargs, confounders = None, max_cluster_size=10, min_samples_per_intron=5, min_samples_per_group=4, min_coverage=0, min_unique_vals = 10): 

    normalize = lambda g: g/g.sum()

    clu, cluster_counts, idx = inp
    
    cluster_start_time = timer()

    cluster_size = cluster_counts.shape[1]

    if cluster_size>max_cluster_size: 
        return(["Too many introns in cluster"])
    if cluster_size <= 1: 
        return(["<=1 junction in cluster"])

    sample_totals=cluster_counts.sum(1)
    samples_to_use=sample_totals>0
    if samples_to_use.sum()<=1: 
        return(["<=1 sample with coverage>0"])
    sample_totals=sample_totals[samples_to_use]
    if (sample_totals>=min_coverage).sum()<=1:
        return(["<=1 sample with coverage>min_coverage"])
    # this might be cleaner using anndata
    x_subset=x[samples_to_use] # assumes one covariate? ah no, covariates handled later 
    cluster_counts=cluster_counts[samples_to_use,]
    introns_to_use=(cluster_counts>0).sum(0)>=min_samples_per_intron # only look at introns used by at least 5 samples
    if introns_to_use.sum()<2:
        return(["<2 introns used in >=min_samples_per_intron samples"])
    cluster_counts=cluster_counts[:,introns_to_use]
    
    # this is the only part that depends on x
    unique_vals, ta = np.unique(x_subset[sample_totals>=min_coverage], return_counts = True)
    if x_subset.dtype.kind in 'OUS': # categorical x
        if (ta >= min_samples_per_group).sum()<2: # at least two groups have this
            return(["Not enough valid samples"])
        x_subset = pd.get_dummies(x_subset, drop_first = True)
        to_drop = (x_subset == 0).all()
        x_subset = x_subset.loc[:, ~to_drop ] # remove empty groups
    else: # continuous x
        if len(unique_vals) < min_unique_vals:
            return(["Not enough valid samples"])
        x_subset = pd.DataFrame( {"x" : x_subset} )
    
    x_only = torch.tensor(x_subset.to_numpy(), **torch_types)
    intercept = torch.ones(x_only.shape[0], 1, **torch_types) 
    if confounders is None:
        x_full = torch.cat((intercept, x_only), axis = 1)
        x_null = intercept
    else:
        these_confounders = torch.tensor(confounders.iloc[samples_to_use].to_numpy(), **torch_types)
        #filter out confounders with no standard deviation
        these_confounders = these_confounders[:,torch.std(these_confounders, dim = 0) != 0.]
        x_full = torch.cat((intercept, these_confounders, x_only), axis = 1)
        x_null = torch.cat((intercept, these_confounders), axis = 1)
    
    y = torch.tensor(cluster_counts, **torch_types)

    #fitting_start_time = timer()
    #cluster_time += (fitting_start_time - cluster_start_time)
    loglr, df, lrtp, null_fit, full_fit, refit_null_flag = dirichlet_multinomial_anova(
        x_full, 
        x_null, 
        y, 
        **kwargs)
    #fitting_end_time = timer()
    #fitting_time += (fitting_end_time - fitting_start_time)

    #results_start_time = timer()


    x_dim = x_subset.shape[1]
    P_null = null_fit.beta.shape[0]
    
    # extract effect sizes
    logef = pd.DataFrame( full_fit.beta[-x_dim:,:].cpu().numpy().T, columns = x_subset.columns ).add_prefix('logef_')
    
    # calculate model-based PSI for each category. For continuous x this will correspond to being +1s.d. from the mean. 
    perturbed = torch.stack( [ normalize( (full_fit.beta[0,:] + full_fit.beta[P_null+i,:]).softmax(0) * full_fit.conc) for i in range(x_dim) ]).T.cpu().numpy()
    perturbed = pd.DataFrame( perturbed, columns = x_subset.columns ).add_prefix('psi_')
    
    junc_results = pd.concat(
        [pd.DataFrame({
            'cluster' : [clu] * y.shape[1], 
            'intron' : idx[idx][introns_to_use].index, 
            'psi_0' : normalize(full_fit.beta[0,:].softmax(0) * full_fit.conc).cpu().numpy()}), 
        logef, 
        perturbed], axis = 1)

    return [
        "Success", 
        {
            'loglr': loglr,
            'null_ll': -null_fit.loss,
            'null_exit_status': null_fit.exit_status,
            'full_ll': -full_fit.loss,
            'full_exit_status': full_fit.exit_status,
            'df': df,
            'p': lrtp
        }, 
        junc_results]


def differential_splicing(counts, x, confounders = None, max_cluster_size=10, min_samples_per_intron=5, min_samples_per_group=4, min_coverage=0, min_unique_vals = 10, device = "cpu", timeit = False, num_cores = 1,  **kwargs):
    '''Perform pairwise differential splicing analysis.

    counts: An [introns] x [samples] dataframe of counts. The rownames must be of the form chr:start:end:cluid. If the counts file comes from the leafcutter clustering code this should be the case already.
    x: A [samples] numeric pandas Series, should typically be 0s and 1s, although in principle scaling shouldn't matter.
    confounders: A [samples] x [confounders] pandas dataframe to be controlled for in the GLM. Factors should already have been converted to a 1-of-(K-1) encoding
    ####, e.g. using model.matrix (see scripts/leafcutter_ds.R for how to do this). Can be None, implying no covariates are controlled for.
    max_cluster_size: Don't test clusters with more introns than this (Default = 10)
    min_samples_per_intron: Ignore introns used (i.e. at least one supporting read) in fewer than n samples (Default = 5)
    min_samples_per_group: Require this many samples in each group to have at least min_coverage reads (Default = 4)
    min_coverage: Require min_samples_per_group samples in each group to have at least this many reads (Default = 20)
    device: Device for pytorch opperations; can be "cpu" or "gpu" (Default = "cpu")
    timeit: Whether or not to return a dictionary of times for the cluster processing and fitting steps (Default = False)
    kwargs: keyword arguments passed to dirichlet_multinomial_anova
    '''

    junc_meta = counts.index.to_series().str.split(':',expand=True).rename(columns = {0:"chr", 1:"start", 2:"end", 3:"cluster"})
    cluster_ids = junc_meta.cluster.unique()
    
    torch_types = { "device" : device, "dtype" : torch.float } # would we ever want float64? 
    
    cluster_time = 0
    fitting_time = 0
    results_time = 0

    idx = [ clu == junc_meta.cluster for clu in cluster_ids ]
    cluster_counts = [ np.array(counts.loc[ i,: ]).transpose() for i in idx ]
    
    task_task = partial(task, torch_types = torch_types, 
            kwargs = kwargs, 
            x = x, confounders = confounders, max_cluster_size=max_cluster_size, min_samples_per_intron=min_samples_per_intron, min_samples_per_group=min_samples_per_group, min_coverage=min_coverage, min_unique_vals = min_coverage)

    with Pool(processes=num_cores) as pool:
        pool_results = pool.map(task_task, zip(cluster_ids, cluster_counts, idx))

    status_df = pd.DataFrame({ "status" : [ g[0] for g in pool_results ]}) 
    status_df.index = cluster_ids

    results = { k:v[1] for k,v in zip(cluster_ids, pool_results) if v[0] == "Success" }
    cluster_table = pd.DataFrame(results.values()) 
    cluster_table.index = results.keys()
    cluster_table['p.adjust'] = robust_fdr(cluster_table['p'], method = 'bh')

    #add failed clusters
    cluster_table = cluster_table.merge(status_df, how = 'outer', left_index = True, right_index = True)
    
    junc_results = [ v[2] for v in pool_results if v[0] == "Success" ]
    junc_table = pd.concat(junc_results, axis=0) # note this should handle missing categories fine
    
    for group in x.unique():
        # excludes baseline which is 0
        if 'psi_' + group in junc_table.columns:
            junc_table['deltapsi_' + group] = junc_table['psi_' + group] - junc_table['psi_0']
    
    time_dict = dict(zip(['cluster_filtering', 'fitting', 'results_processing'], [cluster_time, fitting_time, results_time]))
    
    if timeit:
        return cluster_table, junc_table, status_df, time_dict
    else:
        return cluster_table, junc_table, status_df

def task_junc(clu, cluster_counts, idx, x, torch_types, kwargs, confounders = None, min_samples_per_intron=5, min_samples_per_group=4, min_coverage=0, min_unique_vals = 10): 

    normalize = lambda g: g/g.sum()

    cluster_start_time = timer()

    cluster_size = cluster_counts.shape[1]

    if cluster_size <= 1: 
        return(["<=1 junction in cluster"])

    sample_totals=cluster_counts.sum(1)
    samples_to_use=sample_totals>0
    if samples_to_use.sum()<=1: 
        return(["<=1 sample with coverage>0"])
    sample_totals=sample_totals[samples_to_use]
    if (sample_totals>=min_coverage).sum()<=1:
        return(["<=1 sample with coverage>min_coverage"])
    # this might be cleaner using anndata
    x_subset=x[samples_to_use] # assumes one covariate? ah no, covariates handled later 
    cluster_counts=cluster_counts[samples_to_use,]
    introns_to_use=(cluster_counts>0).sum(0)>=min_samples_per_intron # only look at introns used by at least 5 samples
    if introns_to_use.sum()<2:
        return(["<2 introns used in >=min_samples_per_intron samples"])
    cluster_counts=cluster_counts[:,introns_to_use]
    
    # this is the only part that depends on x
    unique_vals, ta = np.unique(x_subset[sample_totals>=min_coverage], return_counts = True)
    if x_subset.dtype.kind in 'OUS': # categorical x
        if (ta >= min_samples_per_group).sum()<2: # at least two groups have this
            return(["Not enough valid samples"])
    else: # continuous x
        if len(unique_vals) < min_unique_vals:
            return(["Not enough valid samples"])
    
    return ["Success", introns_to_use]


def differential_splicing_junc(counts, x, confounders = None, min_samples_per_intron=5, min_samples_per_group=4, min_coverage=0, min_unique_vals = 10, device = "cpu", timeit = False, num_cores = 1, learn_conc_prior = False, learn_beta_scale_prior = True,  **kwargs):
    '''Perform pairwise differential splicing analysis.

    counts: An [introns] x [samples] dataframe of counts. The rownames must be of the form chr:start:end:cluid. If the counts file comes from the leafcutter clustering code this should be the case already.
    x: A [samples] numeric pandas Series, should typically be 0s and 1s, although in principle scaling shouldn't matter.
    confounders: A [samples] x [confounders] pandas dataframe to be controlled for in the GLM. Factors should already have been converted to a 1-of-(K-1) encoding
    ####, e.g. using model.matrix (see scripts/leafcutter_ds.R for how to do this). Can be None, implying no covariates are controlled for.
    max_cluster_size: Don't test clusters with more introns than this (Default = 10)
    min_samples_per_intron: Ignore introns used (i.e. at least one supporting read) in fewer than n samples (Default = 5)
    min_samples_per_group: Require this many samples in each group to have at least min_coverage reads (Default = 4)
    min_coverage: Require min_samples_per_group samples in each group to have at least this many reads (Default = 20)
    device: Device for pytorch opperations; can be "cpu" or "gpu" (Default = "cpu")
    timeit: Whether or not to return a dictionary of times for the cluster processing and fitting steps (Default = False)
    kwargs: keyword arguments passed to dirichlet_multinomial_anova
    '''

    junc_meta = counts.index.to_series().str.split(':',expand=True).rename(columns = {0:"chr", 1:"start", 2:"end", 3:"cluster"})
    cluster_ids = junc_meta.cluster.unique()
    normalize = lambda g: g/g.sum()
    
    torch_types = { "device" : device, "dtype" : torch.float } # would we ever want float64? 
    
    cluster_time = 0
    fitting_time = 0
    results_time = 0

    print("Preprocessing data")
    y = []
    n = []
    juncs = []
    for clu in tqdm(cluster_ids):
        idx = clu == junc_meta.cluster
        cluster_counts = np.array(counts.loc[ idx,: ]).transpose()
        res = task_junc(clu, cluster_counts, idx, torch_types = torch_types, kwargs = kwargs, x = x, confounders = confounders, min_samples_per_intron=min_samples_per_intron, min_samples_per_group=min_samples_per_group, min_coverage=min_coverage, min_unique_vals = min_coverage)
        if res[0] == "Success": 
            introns_to_use = res[1]
            y_here = cluster_counts[:,introns_to_use]
            n_here = np.broadcast_to(y_here.sum(1)[:, np.newaxis], y_here.shape)
            n.append(torch.tensor(n_here, **torch_types))
            y.append(torch.tensor(y_here, **torch_types))
            juncs.append(junc_meta.cluster[idx][introns_to_use])

    if len(n)==0: 
        raise ValueError("No testable junctions") 
    y = torch.cat(y, dim = 1)
    n = torch.cat(n, dim = 1)
    juncs = pd.concat(juncs)

    if x.dtype.kind in 'OUS': # categorical x
        x = pd.get_dummies(x, drop_first = True)

    x_only = torch.tensor(x.to_numpy(), **torch_types)
    intercept = torch.ones(x_only.shape[0], 1, **torch_types) 
    if confounders is None:
        x_full = torch.cat((intercept, x_only), axis = 1)
        x_null = intercept
    else:
        these_confounders = torch.tensor(confounders.iloc[samples_to_use].to_numpy(), **torch_types)
        #filter out confounders with no standard deviation
        these_confounders = these_confounders[:,torch.std(these_confounders, dim = 0) != 0.]
        x_full = torch.cat((intercept, these_confounders, x_only), axis = 1)
        x_null = torch.cat((intercept, these_confounders), axis = 1)

    sas = bayes_glm.SpikeAndSlabModel(
        gamma_shape = dist.Gamma(2., 1.) if learn_conc_prior else 2., 
        gamma_rate = dist.Gamma(2., 10.) if learn_conc_prior else 0.2, 
        beta_scale = dist.HalfCauchy(1.) if learn_beta_scale_prior else 2. 
    )
    losses_null, losses_full, losses = sas.fit(x_null, x_full, y, n, alpha = 0., num_particles = 1, **kwargs)
    marg_prob, log_bayes_factor = sas.estimate_marginal_posterior(x_null, x_full, y, n, alpha = 1.)
    
    return losses_null, losses_full, losses, pd.DataFrame({"junc":juncs, "prob_diff":marg_prob, "log_bayes_factor":log_bayes_factor})
