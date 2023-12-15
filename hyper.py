import numpy as np
import pandas as pd
import dist, sys, csv
from scipy.stats import hypergeom, false_discovery_control

# X is an object containing 2 pandas data frames of the coordinates for genes and cells (results from running RunMCA() )
# pathways is the GENE MARKER LIST
# p_adjust if TRUE, apply Benjamini Hochberg correction to p-value (values in matrix A from hypergeo distribution)
def RunCellHGT(X: object, pathways, n_features = 200, features = None, dims = range(50), minSize = 10, log_trans = True, p_adjust = True):
    # Find DISTANCES between the coordinates of features (genes) and cells
    # DT will be a pandas data frame, 2d matrix containing distances between genes and cells (rows are genes, cells are columns)
    DT = dist.GetDistances(X)
    print("Distance matrix:" , DT.shape)
    print(DT)
    # Optionally write distance results to CSV file
    #DT.to_csv('Results_csv/CellGeneDistances.csv')
    
    print("ranking genes")
    features = DT.index
    cells = DT.columns
    
    # Target ------------------------------------------------------------------
    # i is a matrix containing the INDICES of the genes with the smallest distance for each cell (column), sorted for each column -> (GENE RANKINGS FOR EACH CELL by row INDICE)
    i = np.argsort(DT, axis=0)[range(n_features), :] 
    # j is a vector of n_features repeated values in a range of the number of columns in the distance matrix. Example: 1111 2222 3333  -> the range is 1 to 3 inclusive but each value is repeated 4 times. Purpose is to be the column indice for the sparse matrix stored in TargetMatrix. i holds row indices and j holds column indices for the locations in the matrix to have value 1. 
    j = np.repeat(np.array(range(len(DT.columns))), n_features) 
    # TargetMatrix_df contains values 0 and 1 representing if the gene is within the cell's gene signature (top 200 genes with closest distance to that cell)
    TargetMatrix = np.zeros(shape=(len(features), len(cells))) # create matrix of all zeros
    j_index = 0
    
    for column_num in range(i.shape[1]):
        for row_num in range(i.shape[0]):
            TargetMatrix[ i[row_num, column_num], j[j_index] ] = 1
            j_index += 1    
    print("TargetMatrix:")
    TargetMatrix = pd.DataFrame(TargetMatrix, index=features, columns=cells)
    #TargetMatrix.to_csv('Results_csv/targetMatrix_py.csv') 
    print(TargetMatrix.shape)
    print(TargetMatrix)
    # ------------------------------------------------------------------------
    
    # Geneset (Gene marker dataset) -----------------------------------------------------------------
    nPathInit = len(pathways) # number of cell types in gene marker dataset
    cell_type_list = pathways.index
    # Row by row, filter the genes in the gene marker dataset to be only the genes inside of the features list (genes from the input gene expression data)
    new_gene_list= []
    for i in range(len(pathways)): # iterate cell type by cell type
        for gene in pathways[i]: # iterate gene by gene for ith cell type
            if gene in features: # check if gene is a gene in the list of unknown cells
                new_gene_list.append(gene)
        pathways[i] = new_gene_list
        new_gene_list = []
    print("Genes that are in the features list: ")
    print(pathways)
    # IMPORTANT GENE MARKER FILTERING STEP: Only keep the cell types (rows) in the gene reference matrix, pathways, where the number of matched genes for a row is >= minSize
    # pathways now stores data in a dictionary: {"cell type name 1": [list of genes], "cell type name 2": [list of genes], ...}
    dict_pathways = {}
    for i in range(len(pathways)):
        if len(pathways[i]) >= minSize: # CHANGE 1 TO minSize
            dict_pathways.update({cell_type_list[i]:pathways[i]}) 
    pathways = dict_pathways
    print("\n", pathways, "\n")
    nPathEnd = len(pathways)
    nFiltPath = nPathInit - nPathEnd
    # Stop Hypergeometric testing if there are no cell types in the gene reference data that have more than minSize (default val is 10) 
    # genes that are present in the input gene expression data containing the unknown cells
    if nPathEnd == 0:
        sys.exit(f"All pathways (possible cell types using the user-chosen gene marker dataset) have less than {minSize} features in common with the data")
    print(f"{nPathEnd} pathways kept for hypergeometric test out of {nPathInit}, {nFiltPath} filtered as less than {minSize} features was present in the data")
    print("\ncalculating features overlap\n")
    # PathwayMat contains row INDICES of where the genes left in the gene marker data (pathways) are located in the distance matrix (DT)
    # PathwayLen contains dictionary with length of gene list for each cell type
    new_gene_list= []
    PathwayLen = {} 
    for key in pathways:
        for gene in pathways.get(key):
            new_gene_list.append(features.to_list().index(gene))
        PathwayLen.update({key : len(new_gene_list)})
        pathways[key] = new_gene_list 
        new_gene_list = []
    PathwayMat = pathways
    print(PathwayMat)
    print(PathwayLen)
    j = [] # j will contain column indices
    length_list = list(PathwayLen.values())
    for i in range(len(PathwayMat)):
        j.extend(np.repeat(i, length_list[i]))    
    # IMPORTANT: PathwayMatrix is a sparse matrix for the GENE SET representing the genes in the distance matrix that are present in the reference gene marker dataset for every cell type.
    # Value is 1 if the gene is present for that cell type. Rows are genes, columns are cell types 
    PathwayMatrix = np.zeros(shape=(len(features), len(PathwayLen)))
    j_index = 0
    for key in PathwayMat:
        for index in pathways.get(key):
            PathwayMatrix[ index, j[j_index] ] = 1
            j_index += 1    
    PathwayMatrix =  pd.DataFrame(PathwayMatrix, index=features, columns=list(PathwayLen.keys()))
    print("PathwayMatrix:" , PathwayMatrix.shape)
    print(PathwayMatrix)
    #PathwayMatrix.to_csv("Results_csv/pathway_matrix.csv")
    # ------------------------------------------------------------------------
    
    # TargetMatrix contains 0s and 1s saying which of the genes are the top 200 closest to cell Y (represented by 1)
    # PathwayMatrix contains 0s and 1s saying which of the genes are present (represented by 1) in the cell type (cell types determined from gene marker dataset used)
     
     
    # Hypergeo ----------------------------------------------------------------
    
    # q: vector of quantiles representing the number of white balls drawn without replacement from an urn which contains both black and white balls.
        # (number of genes that are both present in that cell type and are one of the n_features (200) closest genes to that unknown cell - 1)
    # m: the number of white balls in the urn. (num of genes, for every cell type in gene reference dataset, that ARE in the distance Matrix)
    # n: the number of black balls in the urn. (num of genes in the distance matrix that are NOT in that cell type)
    # k: the number of balls drawn from the urn (default value is 200, representing the number of genes for every unknown cell, from distance matrix, that
        # are being considered for that unknown cell (determined when gene ranking by closest distance of gene x to unknown cell y ))
    
    # q is a pandas data frame where rows are the unknown cells and columns are the cell types from the PathwayMatrix. The values in the frame represent the number of genes that are both present in that cell type and are one of the n_features (200) closest genes to that unknown cell - 1. So if a value in q is 36, that means there were 37 genes that were both in that cell type and were one of the top 200 closest genes. The -1 may be done to help filter unknown cells that have 0 genes that meet the criteria above. 
    # Important understanding: the column of the max value in every row is the cell type that the unknown cell is most closely related to
    
    q = pd.DataFrame(np.dot(np.transpose(TargetMatrix), PathwayMatrix) - 1 , index=TargetMatrix.columns, columns=PathwayMatrix.columns) # GET RID OF THE -1, messing up distribution results
    
    # q.to_csv("hyperdistrib_csv/q.csv", index=False)
    # 112 OF Q's VALUES DIFFER BY a value of 1 when comparing python to R output (when using baron with panglao_pancreas)
    # (most likely because of the varying distance values by decimals between R and python resulting in the ranking of genes to be off)
    print("q:\n", q)
    
    m = PathwayLen # number of genes for each cell type after the filtering
    print("m:\n", m)
    # with open('hyperdistrib_csv/m.csv', 'w', newline='') as file:
    #     writer = csv.writer(file)
    #     writer.writerow(list(m.values())) 

    
    
    n = {}  
    for key in m:
        n.update({key : len(features) - m[key]})
    print("n:\n", n)
    # with open('hyperdistrib_csv/n.csv', 'w', newline='') as file:
    #     writer = csv.writer(file)
    #     writer.writerow(list(n.values())) 
    
    k = n_features
    
    print("performing hypergeometric test\n")
    # Hypergeometric distribution stored in variable A: results in probability values (between 0 and 1) that unknown cell x is cell type y (based off of similarities of genes in gene ranking of the unknown cell and the genes in that cell type)    
    
    m = list(m.values())
    n_list = list(n.values())
    #Python version variables for hyper distrib func
    N = 200
    n = m
    M = [m[i] + n_list[i] for i in range(len(m))] 
    # q is the equivalent of x parameter in hypergeom.pmf()
    q = q.astype("int")

    print("M", M)

    print("performing hypergeometric test\n")
    """
    Hypergeometric distribution stored in variable A: results in probability values (between 0 and 1) that unknown cell x is cell type y (based off of similarities of genes in gene ranking of the unknown cell and the genes in that cell type)
    pmf(k, M, n, N, loc=0) -> Probability mass function.

    PYTHON HYPER DISTRIB FUNCTION:
    M is the total number of objects, 
    n is total number of Type I objects (equivalent to Rs WHITE BALLS). 
    The random variate represents the number of Type I objects in N drawn without replacement from the total population.

    R vs Python
    q = k
    m+n = M
    m = n
    """

    # -> hypergeom distribution performed on every column (cell types) so iterate the columns for q, m, and n  
    A = pd.DataFrame()
    i = 0

    """
    Ex: We have a collection of (m+n) -> (14804 for Baron test input) total genes and 46 of these genes are in this cell type.
    If we choose 200 genes at random (number of genes being considered in gene ranking for an unknown cell), what is probability that 36 of these genes are both found in that cell type and unknown cell's gene ranking?

    M: total num of genes
    n: num of genes in this cell type (that are also found in the list of total genes)
    N: 200 (n_features)
    x: num of genes that are both found in that cell type and unknown cell's gene ranking
    """

    # loop to create A which contains results of hyper distribution for every unknown cell with every cell type
    r = range(0, q.shape[0])
    cells = list(q.columns)
    # iterates cell type by cell type
    for i in range(0, len(cells)):
        # shift distribution of + 1 -> FALSE
        # The + 1 is literally adding 1  to every value in the q column vector (undoing the -1 done when creating q)
        # SO python vs R distribution is different b/c the q values are off by 1 but
        # in python having +1 produces much more accurate results
        prb = hypergeom.pmf(q.iloc[r, i] + 1, M=M[i], n=n[i], N=N)
        A[cells[i]] = prb
    A.index = q.index
    print("A:")
    print(A)
    #A.to_csv("hyperdistrib_csv/A.csv")

    """
    RESULT (IMPORTANT) of RunCellHGT(): 
    -A cell is considered as enriched in those gene sets for which the hypergeometric test p-value is < 1e-02 (0.01)  (-log10 corrected p-value > 2), after Benjamini Hochberg multiple testing correction (p_adjust).
    -The RunCellHGT function will provide the -log10 corrected p-value for each cell and each signature evaluated, so a multi-class evaluation is enabled. When a disjointed classification is required, a cell will be assigned to the gene set with the lowest significant corrected p-value. If no significant hits are found, a cell will remain unassigned.
    """

    A = A.T

    p_adjust = True
    log_trans = True

    # Benjamini Hochberg multiple testing correction:
    # method of Benjamini, Hochberg, and Yekutieli control the false discovery rate, the expected proportion of false discoveries amongst the rejected hypotheses.
    if p_adjust:
        # column by column in A (columns are individual unknown cells),
        # perform Benjamini Hochberg correction
        for col in A.columns:
            A[col] = false_discovery_control(A[col], method='bh')
        print("Benjamini Hochberg Correction:\n", A)
        #A.to_csv("hyperdistrib_csv/A_Hochberg.csv")

    # -log base10 correct p-values in A, for each cell and each signature evaluated
    if log_trans:
        for col in A.columns:
            A[col] = -np.log10(A[col])
        print("log transformation:\n", A)
        #A.to_csv("hyperdistrib_csv/using_my_q_A_minus_log10Correction.csv")

    # A log bas 10 differences: values are different except for the significant value for each unkown cell. 
    # For cells where it seems to not have a cell type, they both have insignificant values so should still result in unnamed cell
    

    # CELL TYPE PREDICTING
    # result will be a dictionary where the key is the unknown cell name and value is the cell type prediction
    # {unknown_cell : cell_type}

    # pancreas_gs_prediction = {}
    # for unknown_cell in A.columns:
    #   pancreas_gs_prediction.append({unknown_cell : nA[unk]})

    HGT_pancreas_gs = A

    # For each cell, assess the signature (cell type) with the lowest corrected p-value (the max -log10 corrected p-value)
    pancreas_gs_prediction = HGT_pancreas_gs.idxmax()

    # For each cell, evaluate if the lowest p-value is significant
    # if the max value in the column in HGT_pancreas_gs is > 2, then cell type prediction is kept. Otherwise, 
    # cell type prediction is unassigned
    for cell in pancreas_gs_prediction.index:
        if np.max(HGT_pancreas_gs[cell]) <= 2:
            pancreas_gs_prediction[cell] = "unassigned"
    
    return pancreas_gs_prediction
    