import pandas as pd
import numpy as np
from mca import RunMCA
from hyper import RunCellHGT
from visual import Visualize_Coordinates
from preprocessData import preprocess
from diff_exp import top_diff_exp_genes
from dist import GetDistances
from filter_genes import filter_genes
from visual import Visualize_Distance_GeneExp
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

# Refer to loadAssay.txt to see how to read from a .rds file containting a seurat object in RStudio, then write data to .csv to be read from python. 
# As of 11/15/2023, there are no python packages that are able to read seurat objects (R specific object) from a .rds file

# Step 1 - Store raw gene count data in a pandas data frame.
# raw gene count data
df = pd.read_csv('BaronMatrixFirstHalf.csv')
df.index = df.iloc[:, 0].tolist() # row names
df = df.iloc[:, 1:] # excluding the column of indexes from the .csv that were added after reading from .csv
df = df.astype('float64')

df2 = pd.read_csv('BaronMatrixLastHalf.csv')
df2.index = df2.iloc[:, 0].tolist() # row names
df2 = df2.iloc[:, 1:] # excluding the column of indexes from the .csv that were added after reading from .csv
df2 = df2.astype('float64')

baron = pd.concat([df,df2])
#------------------------

# Optional Step - Restricting genes to protein-coding genes
HgProteinCodingGenes = pd.read_csv('HgProteinCodingGenes.csv')
baron = baron.loc[list(set(HgProteinCodingGenes['x'].tolist()) & set(baron.index))] # inner join of gene lists
#------------------------

# Optional Step - Remove genes that are expressed in less than n cells.
genes = filter_genes(df=baron, n=5) # input dataframe must be genes as rows and cells as columns
baron = baron.loc[genes]
#------------------------

# Optional Step - Store meta data
baron_meta = pd.read_csv('BaronMeta.csv')
baron_meta.index = baron_meta.iloc[:, 0].tolist()
baron_meta = baron_meta.iloc[:, 1:]
#------------------------

# Step 2 (CRUCIAL STEP for determining accuracy of cell type prediction)- Preprocess the raw gene count data through natural 
# log transformation, normalization, and scaling.

# AnnData for top differentially expressed genes reads features/genes in columns and cells in rows. 
# More importantly, the preprocessing function expects cells to be rows and columns to be genes so transpose
baron = baron.transpose()

# There are genes with a count of zero. When preprocessing the raw gene count data, the log transform 
# results in -inf values as the log of 0 is undefined. Add 1 to the raw gene count data to fix. Also leads to
# better results when normalizing then scaling
baron = baron + 1 

# PREPROCESSING THE RAW GENE COUNT DATA. Data frame must have cells as rows and columns as genes
assay_df = preprocess(baron, order=1, log=False)
#------------------------

# Optional step - Keep the most differentially expressed genes
# values in the meta data that the cells will be grouped by
cell_types = ['Acinar cells',
            'Beta cells',      
            'Delta cells',
            'Pancreatic stellate cells',
            'Ductal cells',
            'Alpha cells',
            'Epsilon cells',
            'Gamma (PP) cells',
            'Endothelial cells',
            'Macrophages',
            'Peri-islet Schwann cells',
            'Mast cells',
            'T cells']

# MUST pass log transformed (RAW gene count data + 1) -> +1 was done earlier to avoid taking log of 0 which is undefined
most_diff_exp_genes_set = top_diff_exp_genes(gene_count_df=np.log(baron), meta_data_df=baron_meta, 
                                             meta_list=cell_types, meta_col_name='cell.type', score_threshold=0)
print("Before gene filtering: ", assay_df.shape)
# filter genes
assay_df = assay_df.loc[:, most_diff_exp_genes_set]
print("After gene filtering: ", assay_df.shape)
print(assay_df)
#------------------------

# Step 3 - Multiple Correspondence Analysis (MCA) dimensionality reduction method
# MCA method expects data frame to have cells as rows and genes as columns
# mca_result is an object containing:
    # cellCoordinates: pandas data frame where rows are cells and columns are j (default value for j is 50) different dimensions
    # geneCoordinates: pandas data frame where rows are genes and columns are j (default value for j is 50) different dimensions
    # X: fuzzy-coded indicator matrix
mca_result = RunMCA(assay_df, j=20) # j specifies number of dimensions for cell and gene coordinates (default, j=50)
# NOTE: upon consecutive calls of generating coorddinates, the signs of the coordinate values may differ meaning one time the value will be positive and 
# another time will be negative. BUT, the distance between gene X and cell Y is still always the SAME. Thus, does not affect the next step in cell 
# identification being to find the distance between genes and cells in a j dimensional space.
#------------------------

# Optional Step - Visualize umap of the cell and gene coordinates in 3D space (condenses j dimensions to 3)  
# Visualize_Coordinates(cellCoordinates=mca_result.cellCoordinates, geneCoordinates=mca_result.geneCoordinates, n_neighbors=10, min_dist=0.5)
#------------------------
    
# Step 4: Find the euclidean distance between the coordinates of genes and cells among the j dimensions.
# DT will be a pandas data frame, 2d matrix containing distances between genes and cells (rows are genes, cells are columns)
DT = GetDistances(cellCoordinates=mca_result.cellCoordinates, geneCoordinates=mca_result.geneCoordinates)
# Optionally generate gene coordinates through barycentric relationship between gene and cell coordinates, then find euclidean distance between them
# DT = GetDistances(cellCoordinates=mca_result.cellCoordinates, geneCoordinates=mca_result.geneCoordinates, X=mca_result.X, barycentric=True)
#------------------------

# Step 5: Acquire gene marker dataset (pancreatic cell-type gene signatures)
panglao = pd.read_csv('https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz',sep='\t')
# restricting the analysis to pancreas specific gene signatues
panglao_pancreas = panglao[panglao['organ'] == "Pancreas"]
# restricting to human specific genes
panglao_pancreas = panglao_pancreas[panglao_pancreas['species'].str.contains("Hs") ]
# get gene marker dataset in format where every index is the cell type then the value in every index is a list of the official gene symbols (genes) for that cell type
# panglao_pancreas becomes a pandas Series
panglao_pancreas = panglao_pancreas[['cell type', 'official gene symbol']]
panglao_pancreas = panglao_pancreas.groupby('cell type')['official gene symbol'].apply(list)
print("Gene marker dataset:\n", panglao_pancreas)
# print(panglao_pancreas.index) # .index gives you list of cell type string names
# can get list of genes for a cell type by row index or ['cell type name']: panglao_pancreas[0] or panglao_pancreas[panglao_pancreas.index[0]]    
#------------------------

# Step 6: HYPERGEOMETRIC Testing
HGT = RunCellHGT(DT=DT, gene_sets=panglao_pancreas, n_genes=200, minSize=10, p_adjust=False, log_trans=False) 
HGT.to_csv("test/cellTypePredictions_Baron/test.csv")
#------------------------

# Optional Step: Plot UMAP of the euclidean distance data frame and label each point (cell) with the predicted/actual cell type
# Create dictionary where the keys are the name of the cell types and the value is a color associated with the plotly package
color_scale = {
        'Acinar cells': 'darkmagenta',
        'Beta cells': 'red',      
        'Delta cells': 'green',
        'Pancreatic stellate cells': 'mediumpurple',
        'Ductal cells': 'orange',
        'Alpha cells': '#0000FF',
        'Epsilon cells': 'hotpink',
        'Gamma (PP) cells': 'lightgreen',
        'Endothelial cells': 'pink',
        'Macrophages': 'yellow',
        'Peri-islet Schwann cells': 'khaki',
        'Mast cells': 'darkgoldenrod',
        'T cells': 'firebrick',
        'unassigned': 'grey'
} 
# UMAP labelled with predicted cell types
Visualize_Distance_GeneExp(df=DT.T, cell_types=HGT.to_list(), color_scale=color_scale, umap=True)
# UMAP labelled with actual cell types
actual_cell_types = pd.read_csv("Baron_ActualCellType.csv").iloc[:, 1].to_list() 
Visualize_Distance_GeneExp(df=DT.T, cell_types=actual_cell_types, color_scale=color_scale, umap=True)
#------------------------