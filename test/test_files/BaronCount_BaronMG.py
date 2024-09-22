import pandas as pd
import numpy as np
from Ramin123455.mca import RunMCA
from Ramin123455.hyper import RunCellHGT
from Ramin123455.visual import Visualize_Coordinates
from Ramin123455.preprocessData import preprocess
from Ramin123455.diff_exp import top_diff_exp_genes
from Ramin123455.dist import GetDistances
from Ramin123455.filter_genes import filter_genes
from Ramin123455.visual import Visualize_Distance_GeneExp


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

# Step 1 - Store raw gene count data in a pandas data frame.
# raw gene count data
count = pd.read_csv('C:/RaminMohammadi/TestingCreatedPackage/Datasets_Mostafa/Baron/Baron_pancreas_DF.csv')

count.index = count.iloc[:, 0].tolist() # row names
count = count.iloc[:, 1:] # excluding the column of indexes from the .csv that were added after reading from .csv
count = count.astype('float64')
print(count)
#------------------------


# Optional Step - Remove genes that are expressed in less than n cells.
genes = filter_genes(df=count, n=5) # input dataframe must be genes as rows and cells as columns
count = count.loc[genes]
#------------------------

# Optional Step - Store meta data, ACQUIRE CELL TYPE LABELS
count_meta = pd.read_csv('C:/RaminMohammadi/TestingCreatedPackage/Datasets_Mostafa/Baron/Baron_Labels.csv')
count_meta.index = count_meta.iloc[:, 0].tolist()
count_meta = count_meta.iloc[:, 1:]
print(count_meta)
#------------------------

# Step 2 (CRUCIAL STEP for determining accuracy of cell type prediction)- Preprocess the raw gene count data through natural 
# log transformation, normalization, and scaling.

# AnnData for top differentially expressed genes reads features/genes in columns and cells in rows. 
# More importantly, the preprocessing function expects cells to be rows and columns to be genes so transpose
count = count.transpose()

# There are genes with a count of zero. When preprocessing the raw gene count data, the log transform 
# results in -inf values as the log of 0 is undefined. Add 1 to the raw gene count data to fix. Also leads to
# better results when normalizing then scaling
count = count + 1 

# PREPROCESSING THE RAW GENE COUNT DATA. Data frame must have cells as rows and columns as genes
count = preprocess(count, order=1, log=True)
#------------------------


# Step 3 - Multiple Correspondence Analysis (MCA) dimensionality reduction method
# MCA method expects data frame to have cells as rows and genes as columns
# mca_result is an object containing:
    # cellCoordinates: pandas data frame where rows are cells and columns are j (default value for j is 50) different dimensions
    # geneCoordinates: pandas data frame where rows are genes and columns are j (default value for j is 50) different dimensions
    # X: fuzzy-coded indicator matrix
mca_result = RunMCA(count, j=50) # j specifies number of dimensions for cell and gene coordinates (default, j=50)
# NOTE: upon consecutive calls of generating coorddinates, the signs of the coordinate values may differ meaning one time the value will be positive and 
# another time will be negative. BUT, the distance between gene X and cell Y is still always the SAME. Thus, does not affect the next step in cell 
# identification being to find the distance between genes and cells in a j dimensional space.
#------------------------

# Optional Step - Visualize umap of the cell and gene coordinates in 3D space (condenses j dimensions to 3)  
Visualize_Coordinates(cellCoordinates=mca_result.cellCoordinates, geneCoordinates=mca_result.geneCoordinates, n_neighbors=10, min_dist=0.5)
#------------------------
    
# Step 4: Find the euclidean distance between the coordinates of genes and cells among the j dimensions.
# DT will be a pandas data frame, 2d matrix containing distances between genes and cells (rows are genes, cells are columns)
DT = GetDistances(cellCoordinates=mca_result.cellCoordinates, geneCoordinates=mca_result.geneCoordinates)
# Optionally generate gene coordinates through barycentric relationship between gene and cell coordinates, then find euclidean distance between them
# DT = GetDistances(cellCoordinates=mca_result.cellCoordinates, geneCoordinates=mca_result.geneCoordinates, X=mca_result.X, barycentric=True)
#------------------------

# Step 5: Acquire gene marker / marker gene dataset (pancreatic cell-type gene signatures)
folder_path = 'C:/RaminMohammadi/TestingCreatedPackage/Datasets_Mostafa/Baron/MG_csv'
# IMPORTANT NOTE: the below cell_types list is defined by you and must be in the samee order as the order of the .csv files being read, and names must match labels exactly
MG_dict = {}
cell_types = ['ac_stellate', 'acinar', 'alpha', 'beta', 'delta', 'ductal', 'endothelial', 'epsilon', 'gamma', 'macrophage', 'mast', 'q_stellate', 'schwann', 't_cell']


MG1 = pd.read_csv('C:/RaminMohammadi/TestingCreatedPackage/Datasets_Mostafa/Baron/MG_csv/Baron_Pancreas_GM_3_ac_stellate.csv')
genes = list(MG1['x']) # list of genes
MG_dict.update({'ac_stellate': genes})

MG2 = pd.read_csv('C:/RaminMohammadi/TestingCreatedPackage/Datasets_Mostafa/Baron/MG_csv/Baron_Pancreas_GM_3_acinar.csv')
genes = list(MG2['x']) # list of genes
MG_dict.update({'acinar': genes})

MG2 = pd.read_csv('C:/RaminMohammadi/TestingCreatedPackage/Datasets_Mostafa/Baron/MG_csv/Baron_Pancreas_GM_3_alpha.csv')
genes = list(MG2['x']) # list of genes
MG_dict.update({'alpha': genes})

MG2 = pd.read_csv('C:/RaminMohammadi/TestingCreatedPackage/Datasets_Mostafa/Baron/MG_csv/Baron_Pancreas_GM_3_beta.csv')
genes = list(MG2['x']) # list of genes
MG_dict.update({'beta': genes})

MG2 = pd.read_csv('C:/RaminMohammadi/TestingCreatedPackage/Datasets_Mostafa/Baron/MG_csv/Baron_Pancreas_GM_3_ductal.csv')
genes = list(MG2['x']) # list of genes
MG_dict.update({'ductal': genes})

MG2 = pd.read_csv('C:/RaminMohammadi/TestingCreatedPackage/Datasets_Mostafa/Baron/MG_csv/Baron_Pancreas_GM_3_endothelial.csv')
genes = list(MG2['x']) # list of genes
MG_dict.update({'endothelial': genes})

MG2 = pd.read_csv('C:/RaminMohammadi/TestingCreatedPackage/Datasets_Mostafa/Baron/MG_csv/Baron_Pancreas_GM_3_epsilon.csv')
genes = list(MG2['x']) # list of genes
MG_dict.update({'epsilon': genes})

MG2 = pd.read_csv('C:/RaminMohammadi/TestingCreatedPackage/Datasets_Mostafa/Baron/MG_csv/Baron_Pancreas_GM_3_gamma.csv')
genes = list(MG2['x']) # list of genes
MG_dict.update({'gamma': genes})

MG2 = pd.read_csv('C:/RaminMohammadi/TestingCreatedPackage/Datasets_Mostafa/Baron/MG_csv/Baron_Pancreas_GM_3_macrophage.csv')
genes = list(MG2['x']) # list of genes
MG_dict.update({'macrophage': genes})

MG2 = pd.read_csv('C:/RaminMohammadi/TestingCreatedPackage/Datasets_Mostafa/Baron/MG_csv/Baron_Pancreas_GM_3_mast.csv')
genes = list(MG2['x']) # list of genes
MG_dict.update({'mast': genes})

MG2 = pd.read_csv('C:/RaminMohammadi/TestingCreatedPackage/Datasets_Mostafa/Baron/MG_csv/Baron_Pancreas_GM_3_q_stellate.csv')
genes = list(MG2['x']) # list of genes
MG_dict.update({'q_stellate': genes})

MG2 = pd.read_csv('C:/RaminMohammadi/TestingCreatedPackage/Datasets_Mostafa/Baron/MG_csv/Baron_Pancreas_GM_3_schwann.csv')
genes = list(MG2['x']) # list of genes
MG_dict.update({'schwann': genes})

MG2 = pd.read_csv('C:/RaminMohammadi/TestingCreatedPackage/Datasets_Mostafa/Baron/MG_csv/Baron_Pancreas_GM_3_t_cell.csv')
genes = list(MG2['x']) # list of genes
MG_dict.update({'t_cell': genes})

MG = pd.Series(MG_dict)
MG.index.name = 'cell type'
MG.name = 'official gene symbol'
print(MG)
#------------------------

# Step 6: HYPERGEOMETRIC Testing
HGT = RunCellHGT(DT=DT, gene_sets=MG, n_genes=200, minSize=1, p_adjust=True, log_trans=False) 
HGT.to_csv("cell_type_predictions_BaronCOUNT_BaronMG.csv")
#------------------------

# Optional Step: Plot UMAP of the euclidean distance data frame and label each point (cell) with the predicted/actual cell type
# Create dictionary where the keys are the name of the cell types and the value is a color associated with the plotly package
color_scale = {
        'acinar': 'darkmagenta',
        'ac_stellate': 'red',      
        'alpha': 'green',
        'beta': 'mediumpurple',
        'delta': 'orange',
        'ductal': '#0000FF',
        'endothelial': 'hotpink',
        'epsilon': 'lightgreen',
        'macrophage': 'pink',
        'mast': 'yellow',
        'q_stellate': 'khaki',
        'schwann': 'darkgoldenrod',
        't_cell': 'cyan',
        'gamma': 'black',
        'unassigned': 'grey'
} 
# UMAP labelled with predicted cell types
Visualize_Distance_GeneExp(df=DT.T, cell_types=HGT.to_list(), color_scale=color_scale, umap=True)
# UMAP labelled with actual cell types
actual_cell_types = list(count_meta['Label'])
Visualize_Distance_GeneExp(df=DT.T, cell_types=actual_cell_types, color_scale=color_scale, umap=True)



# UMAP of count
count = pd.read_csv('C:/RaminMohammadi/TestingCreatedPackage/Datasets_Mostafa/Baron/Baron_pancreas_DF.csv')
count.index = count.iloc[:, 0].tolist() # row names
count = count.iloc[:, 1:] # excluding the column of indexes from the .csv that were added after reading from .csv
count = count.astype('float64')

Visualize_Distance_GeneExp(df=count.T, cell_types=HGT.to_list(), color_scale=color_scale, umap=True)
# UMAP labelled with actual cell types
actual_cell_types = list(count_meta['Label'])
Visualize_Distance_GeneExp(df=count.T, cell_types=actual_cell_types, color_scale=color_scale, umap=True)
#------------------------