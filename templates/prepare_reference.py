#!/usr/bin/env python3

import os.path
import gzip
import pandas as pd
from collections import Counter
from Bio import SeqIO, bgzf
from Bio.Seq import Seq, back_transcribe
from Bio.SeqRecord import SeqRecord

seq_file  = '!{ref}'
taxa_file = '!{taxa}'
host_organism = '!{host_organism}'
target_organism = '!{target_organism}'
out_fasta = '!{out_fasta}'
out_table = '!{out_table}'

########## HANDLE POSSIBLE GZIPPED FASTA ##########
class GZIPManager(object):
    def __init__(self, file_name):
        self.file_name = file_name
      
    def __enter__(self):
        if self.file_name.endswith('.gz'):
            self.f = gzip.open(self.file_name, 'rt')
        else:
            self.f = open(self.file_name, 'r')
        return self.f
  
    def __exit__(self, *args):
        self.f.close()

#################### FUNCTIONS ####################
def cprefix(x):
    """
    Returns the common prefix of a list of strings (x). In this context, they
    represent the kingdom level and lower phylogeny of the miRNA species.

    If no such common prefix exists, returns the label 'cross-kingdom'.
    """
    return os.path.commonprefix(list(x)) or 'cross-kingdom'

def get_category(org_dict, hostorg, targetorg):
    """
    Returns a category string derived from the presence of the hostorg and
    targetorg keys in the org_dict, a dictionary of organisms as keys.
    """
    cattable = {0: 'env', 1: 'host', 2: 'target', 3: 'ambiguous'}
    catbit = 0
    if hostorg in org_dict: catbit += 1
    if targetorg in org_dict: catbit += 2
    return cattable[catbit]


def get_SeqRecord(row):
    """
    Return  a DNA alphabet Bio.SeqRecord object from a collapsed_df row
    """
    seq = Seq(row.seq)
    return SeqRecord(
        back_transcribe(seq),
        id = row.heimirap_id,
        description = f"duplicates={row.count} category={row.category} LCA={row.LCA}",
        )

def get_organism_code (df, s, verbose=True):
    """Get miRBase code from taxonomy dataframe (df) matching string (s)"""
    res = df.loc[
        ( df.organism.str.match(s, case=False) ) |
        ( df.organism_name.str.contains(s, case=False) )
    ]
    if len(res) == 0:
        print(f"Invalid organism name {s}")
        return None
    if len(res) > 1:
        print(f"Ambiguous organism name '{s}' matches:")
        print(res.to_string(index=False))
        return None

    if verbose:
        print(f"Matched organism name: {s}")
        print(res.to_string(index=False))
    return res.organism.values[0]

##################### TAXONOMY #####################

# Make dataframe from taxonomy data

taxa_df = pd.read_csv(taxa_file, sep="\t")

taxa_df.columns = taxa_df.columns.str.replace('#','')
taxa_df.rename(columns={'name': 'organism_name'}, inplace=True)

# Match organism names with taxonomy data
print("Matching host organism name...")
host_organism = get_organism_code(taxa_df, host_organism)
print("Matching target organism name...")
target_organism = get_organism_code(taxa_df, target_organism)

#################### SEQUENCES ####################
count = 0
sequence_data = []
with GZIPManager(seq_file) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        count += 1

        # Get organism (and name) from id
        organism, name = record.id.split('-', maxsplit=1)
        
        # Get accession, organism_name
        mirid, accession, genus, species, desc_name = record.description.split()
        
        # check id and name from description
        if mirid != record.id :
            raise ValueError(f"ID mismatch: {mirid} <> {record.id}")
        if name != desc_name :
            raise ValueError(f"miRNA Name mismatch: {name} <> {desc_name}")
        
        # Append sequence and metadata
        sequence_data.append(
            {
            'id': record.id,
            'name': name,
            'organism': organism,
            'seq': str(record.seq)
            }
        )

print(f"Read {count} sequence records")

# Make dataframe from seqs
seq_df = pd.DataFrame(sequence_data)

##################### MERGING #####################

# Merge seqs dataframe with taxonomy dataframe
seq_df = seq_df.merge(taxa_df, on='organism', validate='many_to_one')

# Group yb sequence, collapsing redundant seqs but keeping metadata using
# aggregation
collapsed_df = seq_df.groupby('seq', as_index=False).agg({
    'id': set,
    'name': set,
    'organism': Counter,
    'tree': cprefix,
    'NCBI-taxid': set
})

# Add a count of duplicates seqs for each entry
collapsed_df['count'] = seq_df.groupby(seq_df.seq).size().tolist()
# Rename the tree column to common_tree
collapsed_df.rename(inplace=True, columns={
    'tree': 'common_tree',
})

# Add least-common ancestor series from last entry in common_tree
collapsed_df['LCA'] = [ x.strip(';').split(';')[-1] for x in collapsed_df['common_tree'] ]

# Generate a unique sequnce ID for each sequence from dataframe index
collapsed_df['heimirap_id'] = collapsed_df.apply(lambda x: f"mirseq{x.name + 1}", axis=1)

# Add the category
collapsed_df['category'] = [
        get_category(c, host_organism, target_organism) for c in collapsed_df['organism']
        ]

# Create a generator of SeqRecord objects
sequences = ( get_SeqRecord(x) for x in collapsed_df.itertuples(index=True) )

# Base output filename
print(f"Writing output file {out_fasta}")

with bgzf.BgzfWriter(out_fasta, 'wb') as outgz:
    SeqIO.write(sequences=sequences, handle=outgz, format="fasta")

# Output table

# explode collapsed_df into 1 row per id, drop collapsed info except LCA
table = collapsed_df.explode('id').drop(['name', 'organism', 'NCBI-taxid', 'common_tree'], axis=1).merge(seq_df, on=['id', 'seq'], validate='1:1')

# select columns and order for output csv
table = table[[
    'id',
    'category',
    'count',
    'LCA',
    'organism_name',
    'organism',
    'NCBI-taxid',
    'tree',
    'heimirap_id',
    'seq',
    ]].sort_values('id')

# A more explanatory name?
table.rename(columns={'count': 'redundancy'}, inplace=True)

# write csv
table.to_csv(out_table, sep=',', index=False)



