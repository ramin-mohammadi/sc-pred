import numpy as np
import pandas as pd
import sys
from scipy.stats import hypergeom, false_discovery_control

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
    along with this program. If not, see <https://www.gnu.org/licenses>."""
    
# Reference: https://github.com/RausellLab/CelliD

"""
Predicts cell type/gene set for each unknown cell using hypergeometric testing. First acquires the gene signature for every cell n 
by keeping the n_genes smallest euclidean distance values per cell in the distance data frame DT. Then finds W_1,W_2, ... W_omega. 
These are the genes in each gene set W_i that are in the genes list P, the genes in the gene expression matrix after any gene filtering 
steps (W_i is a subset of P). The gene sets W_i that have < minSize number of genes are ignored. w is acquired being the intersection of 
the genes between each cell's gene signature and W_i. Then, hypergeometric distribution is performed with a probability mass function 
using variables described on line 120. Probability values (p-values) between each cell and gene set is acquired. Next, optionally adjust 
p-values with Benjamini Hochberg correction on every cell among the gene sets and/or -log base 10 transform p-values (log transform only 
scales the values and has no affect on the prediction results). Finally, cell's gene sets are predictied using rules described on line 155.
input: 
    - DT: pandas data frame, with genes as rows and cells as columns, of the euclidean distances between the coordinates of genes and cells
    - gene_sets: pandas series where each index is labelled with the gene set/cell type name and each value is a list of genes in that gene set/cell type 
    - n_genes: number of closest genes to keep per cell in the euclidean distance data frame, DT, to represent a cell's gene signature
    - minSize: min number of genes that a gene set must have to keep it after taking subset of P
    - p_adjust: if TRUE, apply Benjamini Hochberg correction, adjusting the p-values from hypergeo probability mass function
    - log_trans: if TRUE, apply -log base 10 to p-values after BH correction (not necessary, just upscales small p-values to larger values)
returns a pandas Series:
    - index are the names of the unknown cells and values are the cell type/gene set predictions per cell      
"""

def RunCellHGT(DT: pd.DataFrame, gene_sets: pd.Series, n_genes = 200, minSize = 10, p_adjust = True, log_trans = False):    
    print("\nRanking genes...")
    genes = DT.index
    cells = DT.columns
    
    # Obtain gene signature for each cell by keeping the n_genes closest genes to each cell using the distance matrix DT--------------------------
    # i is a matrix containing the row INDICES of n_genes genes with the smallest distance from each cell (column), sorted least to greatest for each column 
    # -> (GENE RANKINGS FOR EACH CELL by row INDICE)
    i = np.argsort(DT, axis=0)[range(n_genes), :] 
    # gene_signatures contains values 0 and 1 representing if the gene is within the cell's gene signature (top n_genes closest genes to that cell).
    # Each column represents a cell's gene signature with values of 1 to corresponding genes.
    gene_signatures = np.zeros(shape=(len(genes), len(cells))) # create matrix of all zeros
    for column_num in range(i.shape[1]):
        for row_num in range(i.shape[0]):
            gene_signatures[ i[row_num, column_num], column_num ] = 1
    gene_signatures = pd.DataFrame(gene_signatures, index=genes, columns=cells)
    # ------------------------------------------------------------------------
    
    # Gene marker dataset -----------------------------------------------------------------
    # Determine W_1,W_2, ... W_omega. These are the genes in each gene set W_i that are in the genes list P, the genes in the 
    # gene expression matrix after any gene filtering steps (W_i is a subset of P)
    nPathInit = len(gene_sets) # number of cell types/gene sets in gene marker dataset
    cell_type_list = gene_sets.index
    # Row by row (every cell type), filter the genes in the gene marker dataset to be only the genes contained in P
    for i in range(len(gene_sets)): # iterate cell type by cell type
        new_gene_list= []
        for gene in gene_sets.iloc[i]: # iterate gene by gene for ith cell type
            if gene in genes: # check if gene is a gene in the list of unknown cells
                new_gene_list.append(gene)
        gene_sets.iloc[i] = new_gene_list
        
    # GENE MARKER FILTERING STEP: Only keep the cell types/gene sets, where the number of genes 
    # in the cell type W_i after the above step is >= minSize (default is 10)
    # gene_sets now stores data in a dictionary: {"cell type name 1": [list of genes], "cell type name 2": [list of genes], ...}
    dict_gene_sets = {}
    for i in range(len(gene_sets)):
        if len(gene_sets.iloc[i]) >= minSize:
            dict_gene_sets.update({cell_type_list[i] : gene_sets.iloc[i]}) 
    gene_sets = dict_gene_sets
    nPathEnd = len(gene_sets)
    nFiltPath = nPathInit - nPathEnd
    # Stop Hypergeometric testing if there are no cell types in the gene reference dataset that have more than minSize (default val is 10) 
    # genes that are present in the input gene expression data
    if nPathEnd == 0:
        sys.exit(f"All gene sets (possible cell types using the user-chosen gene marker dataset) have less than {minSize} genes in common with the data")
    print(f"{nPathEnd} gene sets kept for hypergeometric testing out of {nPathInit} of the gene sets. {nFiltPath} gene sets were filtered as they had less than {minSize} genes present in the data")
    
    # Place gene sets (W) in columns (cell types as columns and genes in P as rows. Values of 1 represent the genes that are in the cell type/gene set) 
    # gene_sets will contain row INDICES of where the genes in each gene set (W) are located in the distance matrix's (DT) gene list P (since W is a subset of P)
    # W_len contains dictionary with length of gene list for each cell type (practically gene sets W in number format used for the probability mass function)
    W_len = {}         
    genes_list_p = genes.to_list() 
    for key in gene_sets:        
        gene_indice = []
        for gene in gene_sets.get(key):
            gene_indice.append(genes_list_p.index(gene)) # find index of gene in gene list P
        W_len.update({key : len(gene_indice)})
        gene_sets[key] = gene_indice 
    # IMPORTANT: W_binary_matrix contains 0s and 1s representing the genes in each gene set/cell type, W, that are in gene list P.
    # Value is 1 if the gene is present for that cell type. Rows are genes in gene list P, columns are cell types/gene sets
    W_binary_matrix = np.zeros(shape=(len(genes), len(W_len)))
    j = 0
    for key in gene_sets:
        for index in gene_sets.get(key):
            W_binary_matrix[index, j] = 1
        j += 1    
    W_binary_matrix =  pd.DataFrame(W_binary_matrix, index=genes, columns=list(W_len.keys()))
    # ------------------------------------------------------------------------
    
    # gene_signatures contains 0s and 1s representing which of the genes are the top n_genes closest to cell Y (represented by 1)
    # W_binary_matrix contains 0s and 1s representing which of the genes are present (represented by 1) in the cell type (cell types determined from gene marker dataset used)
     
    # Hypergeo ----------------------------------------------------------------
    # Perform probabiliy mass function, optionally benjamini hochberg correction and log transform, 
    # then predict cell type, and assign not confident predictions with "unassigned".
    
    # scipy hypergeom probability mass function: p(k, M, n, N) 
    # n = W_i (number of genes in gene set/cell type after subsetting it from P)
    # k = w_i,n (number of genes in intersection of W_i and cell n's gene signature)
    # M = P (number of genes in input data after any gene filtering steps done before hyper geo testing)
    # N = n_genes (number of genes in cell's gene signature)
    
    n = list(W_len.values())

    # k is a pandas data frame where rows are the unknown cells and columns are the cell types of the gene sets. 
    # IMPORTANT: Each value in k is a w, the number of genes overlapping between a cell's gene signature (top n_genes closest genes to the cell) and a gene set W
    # (the number of genes in each intersection of (cell_n's gene signature and W_i gene set)
    k = pd.DataFrame(np.dot(np.transpose(gene_signatures), W_binary_matrix), index=gene_signatures.columns, columns=W_binary_matrix.columns)
    k = k.astype("int")

    M = len(genes)
    
    N = n_genes

    # Hypergeometric distribution stored in A: probability values (between 0 and 1) that unknown cell x is cell type y 
    # rows are genes in P, columns are cell types/gene sets
    A = pd.DataFrame()
    
    # Probability mass function:
    print("Performing hypergeometric test...\n")    
    i = 0
    cells = list(k.columns)
    # iterate the cell types/gene sets
    for i in range(0, len(cells)):
        prb = hypergeom.pmf(k=k.iloc[:, i], M=M, n=n[i], N=N) # generate probability for all unknown cells with the ith cell type/gene set
        A[cells[i]] = prb # prb is a vector
    A.index = k.index
    print("Probabilities:")
    print(A.T)

    """
    PREDICTING CELL TYPE: 
    - A cell is considered as ENRICHED in those gene sets for which the hypergeometric test p-value is < 1e-02 (0.01) -> 
    only if performing benjamini hochberg correction on probability values or performing neither BH correction nor log transformation,
    OR if performing log transformation on the probability values / after BH correction, -log10 corrected p-value > 2.
    - Simplified: enrichment ranges
        - case 1: p-value < 0.01 : (neither BH correction nor log transform) or (just BH correction on p-values)
        - case 2: p-value > 2 : (only log transform the probability values) or (BH correction then log transform p-values)
    - When a disjointed classification is required (for a cell, there are multiple values/gene sets that are in the enriched range), 
        - case 1: a cell will be assigned to the gene set with the LOWEST significant corrected p-value. 
        - case 2: a cell will be assigned to the gene set with the LARGEST significant corrected p-value. 
    - If no significant hits are found (cell is not enriched in any of the gene sets), a cell will remain unassigned.
    """

    # Benjamini Hochberg multiple testing correction:
    # adjust p-values to control the false discovery rate.
    if p_adjust:
        # perform benjamini correction by rows, each individual unknown cell and its p-values among the gene sets.
        A = pd.DataFrame(false_discovery_control(A, method='bh', axis=1), index=A.index, columns=A.columns)
        print("Benjamini Hochberg Correction:\n", A.T)         

    # CELL TYPE PREDICTING
    # prediction will be a series where the key is the unknown cell name and value is the cell type prediction
    if (p_adjust and not log_trans) or (not p_adjust and not log_trans): # just benjamini correction or neither
        prediction = A.idxmin(axis=1)  # take smallest p value for each unknown cell
        for cell in prediction.index:
            if np.min(A.loc[cell]) >= 0.01: # assign cell's whose most significant p value is not in enriched range as unassigned
                prediction[cell] = "unassigned"
    else: # BH correction then log transform, or no BH corr and log transform
        # -log base10 transform p-values in A
        A = -np.log10(A)
        print("-Log base10 transformation:\n", A.T)   
        prediction = A.idxmax(axis=1) # take largest p value for each unknown cell
        for cell in prediction.index:
            if np.max(A.loc[cell]) <= 2: # assign cell's whose most significant p value is not in enriched range as unassigned
                prediction[cell] = "unassigned"
    print("\nCell type/Gene set predictions:\n", prediction)

    return prediction