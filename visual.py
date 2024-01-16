import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
from umap import UMAP
from sklearn.manifold import TSNE
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

# https://pypi.org/project/umap-learn/
# umap-learn has requirement packages, make sure to pip install those as well (look at above link)

"""
UMAP or TSNE plot each cell by condensing the euclidean distances to each gene (allowing each point in the plot to represent a cell). Then, every cell/point 
is labelled with either the predicted cell type or actual cell type.
    input:
        - df: pandas data frame of euclidean distances between genes and cells where rows are cells and genes are columns
        - cell_types: list of the labelled cell types for each cell whether it be the actual or predicted labels
        - color_scale: dictionary where key is the cell type/gene set name and the value is the color that will be used for that label (color names are from plotly)
        - umap: if True, plots using UMAP. Otherwise, TSNE
        - n_neighbors: (only applies for UMAP) should be a value between 2 and 100. Low values will push UMAP to focus more on local structure by constraining the number 
        of neighboring points considered when analyzing the data in high dimensions, while high values will push UMAP towards representing the big-picture structure 
        while losing fine detail.
        - min_dist: (only applies for UMAP) should be a value between 0.0 and 0.99. Controls how tightly UMAP clumps points together, with low values leading to more 
        tightly packed embeddings. Larger values of min_dist will make UMAP pack points together more loosely, focusing instead on the preservation of the broad topological structure
    output:
        - plot opened in a web browser
"""
def Visualize_Distance_GeneExp(df: pd.DataFrame, cell_types: list, color_scale: dict, umap: True, n_neighbors=10, min_dist=0.5):
    if cell_types is None:
        print('ERROR: list of cell types for individual cells is not given')
        return
    if color_scale is None:
        print('ERROR: dictionary of color map for coloring of cell types is not given')
        return  
    if umap is True: # perform UMAP
        reducer = UMAP(n_neighbors=n_neighbors, min_dist=min_dist, n_components=2, init='random', random_state=0)
        embedding = reducer.fit_transform(X=df)
    else: # perform TSNE
        embedding = TSNE(n_components=2, learning_rate='auto', init='random', perplexity=3).fit_transform(X=df)
        
    # add cell type as a column so can color code the plot by cell types
    plot_df = pd.DataFrame(embedding, columns=["x", "y"])
    plot_df["cell_type"] = cell_types
    
    # scatter plot the UMAP using plotly
    fig = px.scatter(plot_df, x="x", y="y", color="cell_type", color_discrete_map=color_scale)
    fig.show()

"""
Visualize j dimensions of cell and gene coordinates in a 3D space (j dimensions condensed to 3). Utilizes plotly for plotting and 
umap-learn's UMAP function to condense the dimensions.
    input:
        - cellCoordinates: pandas data frame with cells as rows and dimensions as columns
        - geneCoordinates: pandas data frame with genes as rows and dimensions as columns
        - n_neighbors: should be a value between 2 and 100. Low values will push UMAP to focus more on local structure by constraining the number of neighboring 
        points considered when analyzing the data in high dimensions, while high values will push UMAP towards representing the big-picture structure while losing fine detail.
        - min_dist: should be a value between 0.0 and 0.99. Controls how tightly UMAP clumps points together, with low values leading to more tightly packed embeddings. 
        Larger values of min_dist will make UMAP pack points together more loosely, focusing instead on the preservation of the broad topological structure.
    output:
        - opens a link to the umap plot in a web browser. .html file of the plot is saved locally
"""
def Visualize_Coordinates(cellCoordinates: pd.DataFrame, geneCoordinates: pd.DataFrame, n_neighbors=15, min_dist=0.1):
    print("\nCreating UMAP of cell and gene coordinates in 3 dimensional space...\n")
    # construct umap object with parameters you want to be performed (n_components is number of dimensions you want the data to be reduced to) 
    umap_3d = UMAP(n_neighbors=n_neighbors, min_dist=min_dist, n_components=3, init='random', random_state=0)
    
    # Train our reducer, letting it learn about the manifold.
    # Fit X into an embedded space and return that transformed output.
    # the proj variables are n rows by 3 columns NDarray (each column is a dimension in 3D condensed representation from UMAP)
    # the values in the outputted ndarray change b/c the result is a representation of clusters/neighboring points meaning
    # the values outputted from umap are different than the magnitude of the orginal data b/c values used to graph visual clusters to show denisty of data
    
    # perform specifiied UMAP on data
    proj_3d_gene = umap_3d.fit_transform(X=geneCoordinates) 
    proj_3d_cell = umap_3d.fit_transform(X=cellCoordinates)
    
    # plot the results column by column (first column is x, 2nd column is y, and if have 3 dimensions, 3rd column is z)
    # each row in proj_3d_gene corresponds to a gene (gene list based on how was ordered in original data before performing umap),
    # likewise for proj_3d_cell where every row in it corresponds to a cell
    # THEREFORE, every point in the 3D graph is either an individual gene or cell
    fig = go.Figure(data=[
        go.Scatter3d(
            x=proj_3d_gene[:, 0], y=proj_3d_gene[:, 1], z=proj_3d_gene[:, 2],
            mode='markers',
            text=geneCoordinates.index.to_list(),  # text attribute labels every row where every row is an x,y,z coordinate (label every coordinate with gene name)
            name="Gene Coordinate",
            marker = dict(
                size=6,
                opacity=.5,
                colorscale=px.colors.sequential.Bluyl)
        ),
        go.Scatter3d(
            x=proj_3d_cell[:, 0], y=proj_3d_cell[:, 1], z=proj_3d_cell[:, 2],
            mode='markers', 
            text=cellCoordinates.index.to_list(), 
            name="Cell Coordinate",
            marker = dict(
                size=6,
                opacity=.5,
                colorscale=px.colors.sequential.Hot)
        )     
    ])
    
    fig.update_layout(scene = dict(
                    xaxis_title='MCA_1',
                    yaxis_title='MCA_2',
                    zaxis_title='MCA_3'),
                    )
    
    fig.show()
    fig.write_html("coordinate.html", include_plotlyjs=True, full_html=True)