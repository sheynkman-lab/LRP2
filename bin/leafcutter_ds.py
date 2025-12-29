# putting other imports later so help prints quickly
import argparse

# Create a parser
parser = argparse.ArgumentParser(description="LeafCutter differential splicing command line tool. Required inputs:\n <counts_file>: Intron usage counts file. Must be .txt or .txt.gz, output from clustering pipeline.\n <groups_file>: Two+K column file: 1. sample names (must match column names in counts_file), 2. groups. Some samples in counts_file can be missing from this file, in which case they will not be included in the analysis. Additional columns can be used to specify confounders, e.g. batch/sex/age/RIN. Numeric columns will be treated as continuous, so use e.g. batch1, batch2, batch3 rather than 1, 2, 3 if you have a categorical variable.")

# Add command-line arguments
parser.add_argument("counts_file", help="Intron usage counts file. Must be .txt or .txt.gz, output from clustering pipeline.")

parser.add_argument("groups_file", help="Two+K column file: 1. sample names (must match column names in counts_file), 2. groups (currently only two groups, i.e. pairwise, supported. Some samples in counts_file can be missing from this file, in which case they will not be included in the analysis. Additional columns can be used to specify confounders, e.g. batch/sex/age. Numeric columns will be treated as continuous, so use e.g. batch1, batch2, batch3 rather than 1, 2, 3 if you have a categorical variable.")

parser.add_argument("-0", "--baseline_group", default="Control", help="Only for categorical x: the group to which others will be compared [default %(default)s].")
parser.add_argument("-o", "--output_prefix", default="leafcutter_ds", help="The prefix for the output files, <prefix>_cluster_significance.txt (containing test status, log likelihood ratio, degree of freedom, and p-value for each cluster) and <prefix>_effect_sizes.txt (containing the effect sizes and estimated PSI for each intron) [default %(default)s]")
parser.add_argument("-s", "--max_cluster_size", default=float('inf'), type=int, help="Don't test clusters with more introns than this [default %(default)s]")
parser.add_argument("-i", "--min_samples_per_intron", default=5, type=int, help="Ignore introns used (i.e. at least one supporting read) in fewer than n samples [default %(default)s]")
parser.add_argument("-g", "--min_samples_per_group", default=3, type=int, help="Only relevant for categorical x. Require this many samples in each group to have at least min_coverage reads [default %(default)s]")
parser.add_argument("-c", "--min_coverage", default=20, type=int, help="Minimum number of total reads for a cluster to consider it worth testing [default %(default)s]")
parser.add_argument("-u", "--min_unique_vals", default=10, type=int, help="Only relevant for continuous x. Require min_unique_vals unique values after filtering for samples with cluster count > min_coverage [default %(default)s]")
parser.add_argument("-e", "--exon_file", default=None, help="File defining known exons, example in data/gencode19_exons.txt.gz. Columns should be chr, start, end, strand, gene_name. Optional, only just to label the clusters.")
parser.add_argument("--init", default="brr", help="One of One of brr (Bayesian ridge regression), rr (ridge regression), mult (multinomial logistic regression) or `0` (set to 0).")
parser.add_argument("--timeit", default=False, type = bool, help="Whether to print out total time spent at different steps of leafcutter-ds. This is mostly for benchmarking or debugging.")
parser.add_argument("-p", "--num_threads", default=1, type=int, help="Number of threads to use  [default %(default)s]")
  

# Parse the command-line arguments
args = parser.parse_args()
#args = parser.parse_args("-0 male --exon_file gencode.v43.basic.annotation_exons.txt.gz -o python_sex --timeit True Geuvadis_M_vs_F_perind_numers.counts_sample.gz sex_groups.txt".split())

from timeit import default_timer as timer
import_start = timer()
import leafcutter
from leafcutter.differential_splicing.differential_splicing import differential_splicing
import leafcutter.utils

import pandas as pd

import os
import sys
import numpy as np
from sklearn.preprocessing import StandardScaler, OneHotEncoder, scale
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
import_end = timer()

# Access the parsed arguments
print(f"Loading counts from {args.counts_file}")
if not pd.io.common.file_exists(args.counts_file):
    raise FileNotFoundError(f"File {args.counts_file} does not exist")
counts = pd.read_table(args.counts_file, sep = '\s+')

# Loading metadata from groups_file
print(f"Loading metadata from {args.groups_file}")
if not pd.io.common.file_exists(args.groups_file):
    raise FileNotFoundError(f"File {args.groups_file} does not exist")
meta = pd.read_table(args.groups_file, header=None, sep = '\s+')
meta = meta.rename(dict(zip([0, 1], ["sample", "group"])), axis = 1)

# Check if there are more than 2 columns in the metadata DataFrame
confounders = None
if len(meta.columns) > 2:
    # Extract the confounders (columns 3 and onwards)
    confounders = meta.iloc[:, 2:]
    confounders.columns = confounders.columns.astype(str)
    # Initialize a list to store column transformations
    transformations = []
    # Iterate over columns in confounders DataFrame
    for col in confounders.columns:
        if confounders[col].dtype.kind == 'f':
            # Scale continuous confounders
            transformations.append((col, StandardScaler(), [col]))
        else:
            # Treat non-continuous variables as categorical and one-hot encode them
            transformations.append((col, OneHotEncoder(drop='first', sparse_output=False), [col]))
    # Create a column transformer to apply the specified transformations
    column_transformer = ColumnTransformer(transformations)
    # Fit and transform the confounders
    confounders = column_transformer.fit_transform(confounders)
    #remove samples with missing data (not in origial leafcutter as far as I can tell).
    confounders = pd.DataFrame(confounders, index = meta['sample']).dropna()

    all_samples = set(meta.loc[:,'sample'])
    removed_samples = all_samples - set(confounders.index)
    if len(removed_samples) != 0:
        print('Samples removed due to missing values in covariates...')
        print(','.join(list(removed_samples)))

    # Encode the "group" column as numeric (0 and 1)
    meta = meta[meta['sample'].isin(confounders.index)]

    #if permute: numeric_x = np.random.permutation(numeric_x)
    counts = counts[confounders.index]
else:
    # Subset counts to only those in x even without any confounders
    # edge case for if the counts file has more samples than the groups file
    original_count_samples = set(counts.columns)
    if original_count_samples != set(meta['sample']):
        removed_samples = original_count_samples - set(meta['sample'])
        print('Samples removed because they are not in the groups file...')
        print(','.join(list(removed_samples)))
    counts = counts[meta['sample']]
    

scale_factor = 1.

if meta["group"].dtype.kind in 'OUS': # Object Unicode String => Categorical
    # Convert the "group" column in meta to a categorical variable with ordered levels
    unique_values = sorted(set(meta["group"]), key=lambda g: g != args.baseline_group)
    meta["group"] = pd.Categorical(meta["group"], categories=unique_values, ordered=True)

    num_groups_w_enough_samples_1 = (meta["group"].value_counts() >= args.min_samples_per_intron).sum()
    num_groups_w_enough_samples_2 = (meta["group"].value_counts() >= args.min_samples_per_group).sum()
    if (num_groups_w_enough_samples_1 < 2) or (num_groups_w_enough_samples_2 < 2): 
        raise ValueError("There are no groups with enough samples to test. You can reduce min_samples_per_intron (-i option) or min_samples_per_group (-g), but note that we have only carefully checked the calibration of leafcutter p-values down to n=4 samples per group.")

else: # continuous
    scale_factor = np.std(meta["group"], ddof=0)
    meta["group"] = scale(meta["group"])

print("Settings: " + str(args))

print("Running differential splicing analysis...")

setup_end = timer()
res = differential_splicing(counts, meta["group"], confounders = confounders, max_cluster_size = args.max_cluster_size, min_samples_per_intron = args.min_samples_per_intron, min_samples_per_group = args.min_samples_per_group, min_coverage = args.min_coverage, init = args.init, device = "cpu", num_cores = args.num_threads, timeit = args.timeit)

if args.timeit: 
    cluster_table, junc_table, status_df, time_dict = res
else:
    cluster_table, junc_table, status_df = res

print("Fit Time", timer() - setup_end)

wrap_up_start = timer()
cluster_table['cluster'] = cluster_table.index

# get cluster to chromsome mapping example chr2:37649:38658:clu_2_NA
intron_meta = leafcutter.utils.get_intron_meta(counts.index)

# Add gene names to clusters if an exon file is available
exon_file = args.exon_file
if not (exon_file is None):
    print("Loading exons from " + exon_file)
    if os.path.isfile(exon_file):
        exons_table = pd.read_table(exon_file, sep = '\t')
        exons_table['chr'] = leafcutter.utils.add_chr(exons_table.loc[:,'chr'])
        intron_meta['chr'] = leafcutter.utils.add_chr(intron_meta.loc[:,'chr'])
        clu_gene_map = leafcutter.utils.map_clusters_to_genes(intron_meta, exons_table)
        if clu_gene_map is None: 
            cluster_table["genes"] = np.nan
            print("Warning: no clusters matched genes.")
        else: 
            cluster_table = cluster_table.merge(clu_gene_map, on = 'cluster', how = 'left')
    else:
        raise FileNotFoundError(f"File {exon_file} does not exist")
else:
    cluster_table["genes"] = np.nan
    print("No exon_file provided.")

output_prefix = args.output_prefix

clu_to_chrclu = dict(zip(intron_meta["cluster"], intron_meta["chr"] + ':' + intron_meta["cluster"]))
cluster_table['cluster'] = cluster_table['cluster'].map(clu_to_chrclu) 

print('Saving results...')
cluster_table.loc[:,('cluster', 'status', 'loglr', 'df', 'p', 'p.adjust', 'genes')].to_csv(output_prefix + "_cluster_significance.txt", sep = '\t', index = False, na_rep='NA')

#just rearranging junction columns for backwards compatibility
rename_juncs = dict(zip(list(junc_table.columns), list(junc_table.columns)))
rename_juncs['psi_0'] = 'psi_' + args.baseline_group
junc_table.rename(rename_juncs, axis = 1, inplace = True)
logef_cols = list(junc_table.columns[junc_table.columns.str.startswith('logef_')])
psi_cols = list(junc_table.columns[junc_table.columns.str.startswith('psi_')])
deltapsi_cols = list(junc_table.columns[junc_table.columns.str.startswith('deltapsi_')])

junc_table.loc[:,['intron'] + logef_cols + psi_cols + deltapsi_cols].to_csv(output_prefix + "_effect_sizes.txt", sep = '\t', index = False, na_rep='NA')

wrap_up_end = timer()
if args.timeit == True:
    print('Time to import stuff:')
    print(round((import_end - import_start), 3))
    print('Time to setup stuff:')
    print(round((setup_end - import_end), 3))
    print('Time spent filtering clusters:')
    print(round(time_dict['cluster_filtering'], 3))
    print('Time spent fitting model:')
    print(round(time_dict['fitting'], 3))
    print('Time spent summarizing results:')
    print(round(time_dict['results_processing'], 3))
    print('Time spent wrapping up:')
    print(round((wrap_up_end - wrap_up_start), 3))
