import pandas as pd
from mca import RunMCA
from hyper import RunCellHGT
from visual import Visualize_Coordinates
from preprocessData import preprocess
from diff_exp import top_diff_exp_genes

"""
Ramin Mohammadi. main.py: File showing how to use the CellId python package to predict cell type 
on raw gene count data using a gene marker database. The raw gene count data contains gene expression
values of gene x for cell y. The gene marker database contains lists of genes for known cell types.
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
    along with this program.  If not, see <https://www.gnu.org/licenses """


# Refer to loadAssay.txt to see how to read from a .rds file containting a seurat object in RStudio, then write data to .csv to be read from python. 
# As of 11/15/2023, there are no python packages that are able to read seurat objects (R specific object) from a .rds file

# Step 1.1 - Store raw gene count data in a pandas data frame.
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

# meta data (optional)
baron_meta = pd.read_csv('BaronMeta.csv')
baron_meta.index = baron_meta.iloc[:, 0].tolist()
baron_meta = baron_meta.iloc[:, 1:]

# Optional Step - Restricting genes to protein-coding genes
HgProteinCodingGenes = pd.read_csv('HgProteinCodingGenes.csv')
baron = baron.loc[list(set(HgProteinCodingGenes['x'].tolist()) & set(baron.index))] # inner join of gene lists

# Step 2 - Preprocess the raw gene count data through normalization, natural log transformation, and scaling.

# AnnData reads features/genes in columns and cells in rows. ALSO, the preprocessing 
# function expects cells to be rows and columns to be genes so transpose
baron = baron.transpose()

# There are genes with a count of zero. When preprocessing the raw gene count data, the log transform 
# results in -inf values as the log of 0 is undefined. Add 1 to the raw gene count data to fix. Also leads to
# better results when normalizing then scaling
baron = baron + 1 

# PREPROCESSING THE RAW GENE COUNT DATA. Data frame must have cells as rows and columns as genes
assay_df = preprocess(baron, scale=False)
assay_df = pd.DataFrame(assay_df, index=baron.index, columns=baron.columns)

# column values in meta data you will group the cells by to filter genes by most differentially expressed
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

# MUST pass preprocessed version of gene count data frame (use preprocess() function beforehand)
most_diff_exp_genes_set = top_diff_exp_genes(preprocessed_count_df=assay_df, meta_data_df=baron_meta, num_top_genes_per_group=1100, 
                                             meta_list=cell_types, meta_col_name='cell.type')

assay_df = assay_df.transpose() # MCA method expects rows to be genes and columns to be cells

filtered_df = assay_df.loc[most_diff_exp_genes_set]
print(filtered_df)
print(filtered_df.shape)
print(assay_df.shape)
#------------


# If want to test with small input, use assay_df_smaller 
# (change the index values as needed). Currently uses first 400 rows and first 200 columns
#assay_df_smaller = assay_df.iloc[:400, :200]

# Perform MCA on input data
# mca_result is an object containing:
    # featuresCoordinates: pandas data frame where rows are genes and columns are nmcs (default value for nmcs is 50) different dimensions
    # cellsCoordinates: pandas data frame where rows are cells and columns are nmcs (default value for nmcs is 50) different dimensions
    # stdev: numpy array containing the singular values created during MCA when performing Singular Value Decomposition
mca_result = RunMCA(assay_df)
# mca_result = RunMCA(filtered_df)

# MCA on smaller input data
# mca_result = RunMCA(assay_df_smaller)

# NOTE: if specifying a value for nmcs in RunMCA(), must provide same value for dims when calling GetDistances()
# Run MCA using only specified features (here is first 60 features) and specifying nmcs
#mca_result = RunMCA(assay_df, features=assay_df.index.tolist()[:60], nmcs=30)


# FUTRUE: THIS CAN BE DONE WITH PARALLEL PROCESSING b/c the plotting of the coordinates takes some time to perform
# and the result is an independent task not need for continuing cell type prediction calculations
#Visualize_Coordinates(mca_result) #-> TAKES LONG TIME
# exit()
# IMPORTANT NOTE: 
    # The orientation of the coordinates plotted is not consistent meaning the absolute values of the coordinates, if ran multiple times, will be the same, but the signs of the coordinate values may differ meaning one time the value will be positive and another time will be negative.
    # BUT, with this in mind, the shape of the scatterplot, no matter the orientaton, is always the same meaning the distance between gene X and cell Y is still always the SAME. Thus, does not affect the next step in cell identification being to find the distance between genes and cells in an nmcs dimensional space.


# Optionally print results
# print("Singular values (stdev): ", mca_result.stdev)
# print("Features coordinates: ", mca_result.featuresCoordinates.shape, mca_result.featuresCoordinates)
# print("Cells coordinates: ", mca_result.cellsCoordinates.shape, mca_result.cellsCoordinates)


# Optionally write mca_result coordinates in form of data frames to CSV files to see entire results
# mca_result.featuresCoordinates.to_csv('Results_csv/featuresCoordinates.csv', index=False)
# mca_result.cellsCoordinates.to_csv('Results_csv/cellsCoordinates.csv', index=False)



# HYPERGEOMETRIC Testing

# get gene marker dataset (pancreatic cell-type gene signatures)
panglao = pd.read_csv('https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz',sep='\t')
# restricting the analysis to pancreas specific gene signatues
panglao_pancreas = panglao[panglao['organ'] == "Pancreas"]
# restricting to human specific genes
panglao_pancreas = panglao_pancreas[panglao_pancreas['species'].str.contains("Hs") ]
# get gene marker dataset in format where every row indexes (row names) are the cell type then every row is a list of the official gene symbols (genes) for that cell type

# panglao_pancreas becomes a pandas Series
panglao_pancreas = panglao_pancreas[['cell type', 'official gene symbol']]
panglao_pancreas = panglao_pancreas.groupby('cell type')['official gene symbol'].apply(list)
print("Gene marker dataset: ")
print(panglao_pancreas.shape)
print(panglao_pancreas)
#print(panglao_pancreas[panglao_pancreas.index[0]])
#print(panglao_pancreas.index) # .index gives you list of cell type string names
# can get list of genes for a cell type by row index or ['cell type name']: panglao_pancreas[0] or panglao_pancreas[panglao_pancreas.index[0]]

HGT = RunCellHGT(mca_result, pathways=panglao_pancreas)

print("Cell type predictions:\n", HGT)
# HGT.to_csv("test/cellTypePredictions_Baron/PartialSVDpredictions.csv")
# HGT.to_csv("test/cellTypePredictions_Baron/Diff_Exp_Py_Predictions.csv")
HGT.to_csv("test/cellTypePredictions_Baron/test.csv")












