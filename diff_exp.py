import pandas as pd
import csv
import anndata as ad
import scanpy as sc

"""
    Copyright (C) 2023  Ramin Mohammadi

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses>."""

# scanpy gene ranking (DE)
# https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.rank_genes_groups.html

"""
Returns the most differentially expressed genes per meta data grouping of the cells using scanpy's rank_genes_groups() function.
All of the genes are given a score per cell group. Only the genes that are above the score_threshold are kept (default is 0).
Meta data per cell is required for the rank_genes_groups() function to base a varying condition off of among the gene counts
MUST pass log transformed version of the RAW gene count data frame in gene_count_df

input:
    - gene_count_df: log transformed version of the RAW gene count data where cells are rows and genes are columns
    - meta_data_df: data frame containing meta data per cell (cells are rows and columns are different meta data)
    - meta_list: the list of all possible values in the meta data column used
    - meta_col_name: the column name of the meta data being used in the meta data frame
    - score_threshold: genes above this value when ranked will be kept
    - show_top_genes_plot: scanpy plot of the ranked genes per group
    - write_top_genes_to_csv: if True, will write the the genes kept to a .csv file
output:
    - list of genes (gene names) that were kept
"""

def top_diff_exp_genes(gene_count_df: pd.DataFrame, meta_data_df: pd.DataFrame, meta_list: list, 
                        meta_col_name: str, score_threshold=0, show_top_genes_plot=False, write_top_genes_to_csv=False):
    print("Finding most differentially expressed genes...\n")
    # Ann data, .obs should hold the meta data and .var should hold gene names,
    # .X holds log transformed raw gene count data
    adata = ad.AnnData(X=gene_count_df) 
    adata.obs_names = gene_count_df.index.to_list()
    adata.var_names = gene_count_df.columns.to_list()
    adata.obs[meta_col_name] = pd.Categorical(meta_data_df[meta_col_name])
    sc.tl.rank_genes_groups(adata, groupby=meta_col_name) # group by is required argument
        
    if show_top_genes_plot:
        sc.pl.rank_genes_groups(adata)
            
    # genes are ordered/ranked per group of cells by their given score so get top x genes per cell group
    most_diff_exp_genes_set = []
    for cell_group in meta_list:
        dedf = sc.get.rank_genes_groups_df(adata, group=cell_group)
        print("Group: ", cell_group, "\n", dedf)
        for index in dedf.index:
            # only keep genes that were given a score above the threshold to avoid including genes that were not not highly expressed in a group
            if dedf['scores'][index] > score_threshold: 
                most_diff_exp_genes_set.append(dedf['names'][index])

    # Get rid of duplicate genes listed by acquiring the SET of the list of genes
    most_diff_exp_genes_set = list(set(most_diff_exp_genes_set))
 
    if write_top_genes_to_csv:
        with open('most_diff_exp_genes.csv', "w", newline="") as file:
            writer = csv.writer(file)
            writer.writerow(most_diff_exp_genes_set)
            
    # most differentially expressed genes in different groups may overlap (present in more than one group)
    return most_diff_exp_genes_set