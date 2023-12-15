import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
from umap import UMAP
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import seaborn as sns
# https://pypi.org/project/umap-learn/
# umap-learn has requirement packages, make sure to pip install those as well (look at above link)

def Visualize_Distance_GeneExp(df: pd.DataFrame, cell_types:pd.DataFrame, color_scale:dict, umap: True):
    if cell_types is None:
        print('ERROR: No pandas data frame given containing list of cell types for individual cells')
        return
    if color_scale is None:
        print('ERROR: No dictionary of color map given for coloring of cell types')
        return  
    if umap is True: # perform umap
        reducer = UMAP(n_neighbors=10, min_dist=0.5, n_components=2, init='random', random_state=0)
        embedding = reducer.fit_transform(X=df)
    else: # perform TSNE
        embedding = TSNE(n_components=2, learning_rate='auto',
                   init='random', perplexity=3).fit_transform(X=df)
    
    # add cell type as a column so can color code the plot by cell type
    plot_df = pd.DataFrame(embedding, columns=["x", "y"])
    plot_df["cell_type"] = cell_types.to_list()
    
    #using plotly
    fig = px.scatter(plot_df, x="x", y="y", color="cell_type", color_discrete_map=color_scale)
    fig.show()


# this method directly uses umap-learn's umap function 
# UMAP is going to help show density of the data
def Visualize_Coordinates(mca_result):
    # Genes is df1, Cells is df2
    df1 : pd.DataFrame = mca_result.featuresCoordinates
    df2 : pd.DataFrame = mca_result.cellsCoordinates
    
    # construct umap object with parameters you want to be performed (n_components is number of dimensions you want the data to be reduced to) 
    # important umap parameters to adjust: n_neighbors, min_dist
    # default values: n_neighbors=15, min_dist=0.1
    # for the coordinate representation, we want a close local view with points spread out not so tightly packed so,
    # go for lower n_neighbors and higher min_dist value
    umap_3d = UMAP(n_neighbors=10, min_dist=0.5, n_components=3, init='random', random_state=0)
    
    
    # Train our reducer, letting it learn about the manifold.
    # Fit X into an embedded space and return that transformed output.
    # the proj variables are n rows by 3 columns NDarray (each column is a dimension in 3D condensed representation from UMAP)
    # But why do the values from df1 and df2 become larger values after fit_transform? 
    # -> IMPORTANT: the values in the outputted ndarray change b/c the result is a representation of clusters/neighboring points meaning
    # the values outputted from umap are different than the magnitude of the orginal data b/c values used to graph visual clusters to show denisty of data
    # perform specifiied UMAP on data
    proj_3d_gene = umap_3d.fit_transform(X=df1) 
    proj_3d_cell = umap_3d.fit_transform(X=df2)
    
    # plot the results column by column (first column is x, 2nd column is y, and if have 3 dimensions, 3rd column is z)
    # color code the graph by what you visually want to see (grouping color has nothing to do with how umap is performed)
    # TO DO SO: must understand result. Here, each row in proj_3d_gene corresponds to a gene (gene list based on how was ordered in original data before performing umap),
    # likewise for proj_3d_cell where every row in it corresponds to a cell
    # THEREFORE, every point in the 3D graph is either an individual gene or cell
    fig = go.Figure(data=[
        go.Scatter3d(
            x=proj_3d_gene[:, 0], y=proj_3d_gene[:, 1], z=proj_3d_gene[:, 2],
            mode='markers',
            text=df1.index.to_list(),  #text attribute labels every row where every row is an x,y,z coordinate (label every coordinate with gene name)
            name="Gene Coordinate",
            marker = dict(
                size=6,
                opacity=.5,
                colorscale=px.colors.sequential.Bluyl)
        ),
        go.Scatter3d(
            x=proj_3d_cell[:, 0], y=proj_3d_cell[:, 1], z=proj_3d_cell[:, 2],
            mode='markers', 
            text=df2.index.to_list(), 
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
    