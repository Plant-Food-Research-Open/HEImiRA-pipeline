#!/usr/bin/env python3
import pysam
import pandas as pd
import glob
from collections import Counter
from multiprocessing import Pool
import os.path as p

# floats in output tables
float_accuracy = '%.4f'

# Nextflow vars
pool_size = int(!{task.cpus})

meta_csv = '!{metadata}'
bamfiles = glob.glob('./*.bam')
clean_names = int(!{clean_names})
output_prefix = '!{output_prefix}'

print(f"pool size: {pool_size}")
print(f"metadata file: {meta_csv}")
print("bam files:\n", "\n".join(bamfiles))

def lc_suffix(s):
    """
    removes longest common suffix of a list of strings

    returns list with the longest common suffix removed form each string
    """
    min_len = min([len(x) for x in s])
    i = 0
    while i < min_len:
        t = [x[i:] for x in s]
        if len(set(t)) == 1:
            break
        i += 1
    return [x[0:i] for x in s]
    
def su_prefix(s):
    """
    Shortest unique prefix of a list of strings

    returns list of each input string's shortest unique prefix
    """
    slen = len(s)
    print(slen)
    min_len = min([len(x) for x in s])
    i = 1
    while i < min_len:
        t = [x[0:i] for x in s]
        if len(set(t)) == slen:
            break
        i += 1
    return [x[0:i] for x in s]

def get_counts(bf):
    """
    Get counts from bam files
    """
    sam = pysam.AlignmentFile(bf, "rb")
    # Get raw mapped counts for each seq(contig) from index stats
    return pd.Series({
        x.contig : x.mapped for x in sam.get_index_statistics()
    }, name=p.basename(bf)
    )

# Get counts in parallel
with Pool(processes=pool_size) as pool:
    counters = pool.map(get_counts, bamfiles)

# combine the counters from each bam file counts into a dataframe
counts = pd.concat(counters, axis=1, verify_integrity=True)
counts.index.name = 'heimirap_id'
# Clean up column names
if clean_names > 0: # level 1: basic
    counts.columns = lc_suffix(counts.columns)
if clean_names > 1: # level 2: aggressive
    counts.columns = su_prefix(counts.columns)

# Save counts
counts.to_csv(f"{output_prefix}_counts.csv",sep=',', index=True)

meta_df = pd.read_csv(meta_csv)
# we only need the flat part to reference back to the reference table
meta_df = meta_df[['heimirap_id', 'category', 'redundancy', 'LCA', 'seq']].drop_duplicates()
meta_df.set_index('heimirap_id', inplace=True)

# Calculate total and mean for each feature
counts_total = pd.Series(counts.sum(axis=1), name='counts_total')
counts_mean = pd.Series(counts.mean(axis=1), name='counts_mean')

# Normalise to total mapped count
norm = counts.div(counts.sum(axis=0), axis=1).add_prefix('norm_')
norm_total = pd.Series(norm.sum(axis=1), name='norm_total')
norm_mean = pd.Series(norm.mean(axis=1), name='norm_mean')

# get total host-mapped counts
category_counts = pd.concat([counts, meta_df['category']], axis=1, verify_integrity=True)
host_sum = category_counts[category_counts['category'] == 'host'].drop('category', axis=1).sum(axis=0)

# Normalise to total host-mapped count
host_norm = counts.div(host_sum, axis=1).add_prefix('host_norm_')
host_norm_total = pd.Series(host_norm.sum(axis=1), name='host_norm_total')
host_norm_mean = pd.Series(host_norm.mean(axis=1), name='host_norm_mean')

# Concat everything together
data = pd.concat([
    counts_total, counts_mean,
    norm_total, norm_mean,
    host_norm_total, host_norm_mean,
    meta_df,
    ],
    axis=1, verify_integrity=True
)
data.index.name = 'heimirap_id'

# Save summary
data.to_csv(f"{output_prefix}_counts_summary.csv",sep=',', index=True, float_format=float_accuracy)

# Save normalised counts
norm.to_csv(f"{output_prefix}_norm_counts.csv",sep=',', index=True, float_format=float_accuracy)
host_norm.to_csv(f"{output_prefix}_host_norm_counts.csv",sep=',', index=True, float_format=float_accuracy)






