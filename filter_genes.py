import numpy as np
import pandas as pd
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
        
    
"""
Genes expressing at least one count in less than n cells are removed.

input:
    - df: raw gene count data frame with genes as rows and cells as columns
    - n: min number of cells that the gene has to be expressed in
output:
    - list of genes kept
"""
def filter_genes(df: pd.DataFrame, n: int):
    print(f"Removing genes that are expressed in less than {n} cells...")
    genes = []
    df = df.map(lambda x: 1 if x > 0 else 0) # turn all gene count values > 0 to 1
    for gene in df.index:
        if np.sum(df.loc[gene]) >= n: genes.append(gene) # if sum of a row >= n (more than n cells express that gene)
    return genes