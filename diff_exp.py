import pandas as pd
import csv
import anndata as ad
import scanpy as sc

#scanpy gene ranking (DE)
#https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.rank_genes_groups.html

# Returns the x most differentially expressed genes per meta data grouping of the cells using scanpy's rank_genes_groups() function 
# where x=num_top_genes_per_group. Meta data per cell is required for the rank_genes_groups() function to base a 
# varying condition off of among the gene counts

# MUST pass preprocessed version of gene count data frame

def top_diff_exp_genes(preprocessed_count_df: pd.DataFrame, meta_data_df: pd.DataFrame, num_top_genes_per_group: int, meta_list: list, 
                       meta_col_name: str, show_top_genes_plot=False, write_top_genes_to_csv=False):
    print("Finding most differentially expressed genes...\n")
    # Ann data, .obs should hold the meta data and .var should hold gene names,
    # .X holds preprocessed raw gene count data
    ann = ad.AnnData(X=preprocessed_count_df) 
    ann.obs_names = preprocessed_count_df.index.to_list()
    ann.var_names = preprocessed_count_df.columns.to_list()
    ann.obs[meta_col_name] = pd.Categorical(meta_data_df[meta_col_name])
    sc.tl.rank_genes_groups(ann, groupby=meta_col_name) # group by is required argument
    
    if show_top_genes_plot:
        sc.pl.rank_genes_groups(ann)
        
    top_genes_per_cell_type = ann.uns['rank_genes_groups']['names']

    # results are grouped by cell type so get top x genes per cell type and use those genes
    most_diff_exp_genes_set = []
    for cell_type in meta_list:
        most_diff_exp_genes_set.extend(top_genes_per_cell_type[cell_type][:num_top_genes_per_group])
 
    if write_top_genes_to_csv:
        with open('most_diff_exp_genes.csv', "w", newline="") as file:
            writer = csv.writer(file)
            writer.writerow(most_diff_exp_genes_set)
    # most differentially expressed genes in different groups may overlap (present in more than one group)
    # Get rid of duplicate genes listed by returning SET of the list of genes
    return list(set(most_diff_exp_genes_set))