# CelliD-Python

## Description
Python package for cell identity recognition at individual cell level from single-cell RNA-seq data.<br/><br/>
Cell types for unknown cells are predicted, through a statistical approach, using raw gene expresssion <br/>
data among a list of cells, and a gene marker dataset consisting of a list of genes found in known cell types.<br/><br/>
CelliD-Python is based off of the
- R software: https://github.com/RausellLab/CelliD
- CelliD statistical method presented in the article: <a href="https://www.nature.com/articles/s41587-021-00896-6.epdf?sharing_token=cb8TdGrz0o3PXjXn_wZGCdRgN0jAjWel9jnR3ZoTv0Oa3WzvJLtg4J6wv_eRGblv7pCmV-VB-3abW6uWDvAeOER7rbNPidd1IsRjFITIK8SJ_d0RrfACjtlZFkN4l3DDZLXWnaDW2XZDF1uZ-2DWCHQNkva9vqKjz708F5zU2FU%3D">Gene signature extraction and cell identity recognition at the single-cell level with Cell-ID, Nature Biotechnology 2021</a>



## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [API](#api)
- [Motivation](#motivation)
- [Credits](#credits)
- [License](#license)

## Installation
```
pip install CelliD-Python
```

## Usage
Any data used in the examples can be found in the data folder.
### Step 1.1 - Store raw gene count data in a pandas data frame. 
Genes are rows and columns are cells. Row and column names MUST be included. Values will represent gene expression of gene X in cell Y. Here, the Baron pancreas single-cell RNA-seq data set provided in <a href="https://www.sciencedirect.com/science/article/pii/S2405471216302666?via%3Dihub">Baron et al. 2016</a> is being utilized.<br>
<br/>
```
import pandas as pd
from mca import RunMCA
from hyper import RunCellHGT
from visual import Visualize_Coordinates
from preprocessData import preprocess
from diff_exp import top_diff_exp_genes

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
```
Here, the data was stored in an R specific object (within a .rds file), so the data was read from RStudio then written to .csv files to be read from python. The raw gene count matrix was too large to be written to a single .csv file in one write so had to do 2 writes along with 2 reads in python and concatenate the data read to one pandas data frame. NOTE, some of the lines of code above are specific to this situation.

```
print(baron)
```

gene count data frame: <br/>
![alt text](assets/images/raw_gene_count.png) <br/>

### Optional Step 1.2 - Restricting the analysis to protein-coding genes
Using only genes found in the HgProteinCodingGenes list obtained from BioMart Ensembl release 100, version April 2020 (GrCH38.p13 for human, and GRCm38.p6 for mouse). Here the HgProteinCodingGenes was originally read from RStudio, after loading the CelliD library, using: data("HgProteinCodingGenes"). Then, written to a .csv file to be read from python.
```
# Optional Step - Restricting genes to protein-coding genes
HgProteinCodingGenes = pd.read_csv('HgProteinCodingGenes.csv')
baron = baron.loc[list(set(HgProteinCodingGenes['x'].tolist()) & set(baron.index))] # inner join of gene lists
```
```
print(baron)
```
gene count data frame: <br/>
![alt text](assets/images/raw_gene_count_protein.png) <br/>

### Optional Step 1.3 - Store meta data
You can store meta data for the cells in a pandas data frame which can be used to group the cells in the following step for gene filtering. In the meta data frame, rows are cells and columns are different meta data. Here, the data was stored in an R specific object (within a .rds file), so the data was read from RStudio then written to .csv files to be read from python.
```
# meta data (optional)
baron_meta = pd.read_csv('BaronMeta.csv')
baron_meta.index = baron_meta.iloc[:, 0].tolist() # row names
baron_meta = baron_meta.iloc[:, 1:]  # excluding the column of indexes from the .csv that were added after reading from .csv
```
meta data frame: <br/>
![alt text](assets/images/meta_data.png)<br/>



### Step 2 (CRUCIAL STEP for determining accuracy of cell type prediction)- Preprocess the raw gene count data through natural log transformation, normalization, and scaling.
&nbsp;&nbsp;&nbsp;&nbsp; Using the preprocess function which utilizes <a href="https://numpy.org/doc/stable/reference/generated/numpy.log.html">numpy's natural log function</a> along with <a href="https://scikit-learn.org/stable/modules/preprocessing.html">sklearn's normalizing and scaling functions</a>, optionally choose which preprocessing steps to be done and the order of them. Can either be normalizing, natural log transformation, then scaling, OR natural log transformation, normalizing, then scaling. From there, filter which of the steps you want to perform by specifying False for the preprocess() function attributes. Can optionally add 1 as well for the log, norm, scale order after performing log transformation (may improve results).

<br>&nbsp;&nbsp;&nbsp;&nbsp; Log transforming helps reduce the range of the values to avoid computer arithemtic overflow errors.  Normalizing makes values between 0 and 1, and scaling helps ensure that all features contribute equally to the model and avoid the domination of features with larger values.
<br>&nbsp;&nbsp;&nbsp;&nbsp;"There is no guarantee that the log-transformation will reduce skewness and make the data a better approximation of the normal distribution" and "log transformation can often increase – not reduce – the variability of data whether or not there are outliers." (<a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4120293/">Feng C, Wang H, Lu N, Chen T, He H, Lu Y, Tu XM. Log-transformation and its implications for data analysis. Shanghai Arch Psychiatry. 2014</a>). Therefore, best to try different combinations of preprocessing steps to see what works best for your data. <a href="#results">Here are the results from different combinations of preprocessing of the Baron single cell pancreas data.</a>

<br>&nbsp;&nbsp;&nbsp;&nbsp; If performing log transformation first, you can avoid problems with the log transformation turning zeros into -infinity values by adding 1 to the gene count data frame (this problem occurs if taking the log of 0 as this is undefined). Would recommend experimenting with adding 1 if doing normalization first as well. The gene count data frame is transposed since the preprocess function expects the input data frame to consist of cells as rows (samples) and genes as columns (features), allowing for the values for genes to be normalized across every cell.
```
# AnnData reads features/genes in columns and cells in rows. ALSO, the preprocessing 
# function expects cells to be rows and columns to be genes so transpose
baron = baron.transpose()
baron = baron + 1
assay_df = preprocess(baron, order=1, )
assay_df = pd.DataFrame(assay_df, index=baron.index, columns=baron.columns)
```

### Optional step - Filter genes based off of meta data grouping of cells.
The top x most differentially expressed genes (features) per group will be kept. 
```
wefweffrre
```

### Step 3 - Multiple Correspondence Analysis (MCA) dimensionality reduction method
assay_df = assay_df.transpose() # must put data in orientation: genes as rows, cells as columns for MCA method

### Optional step - Plot gene and cell coordinates 
Out of the default 50 dimensions, they are condensed into 3 dimensions using UMAP method for visualization

### Step 4 - Gene Marker Database

### Step 5 - Hypergeometric Testing

### Optional step - Plot cell type predictions and optionally actual cell type predictions
Plot the Umap of the distance matrix where each point represents a cell. They are labelled with the predicted cell type and another umap with the actual expected cell type (this is optional depending on meta data you have available) 

## API
- preprocess(pandas data frame)
- top_diff_exp_genes()

## Motivation
- Translating R package to python helped do...
- Added feature to take top differentially expressed genes based off of grouping of cells by some kind of meta data

## Results
The following are number of cells whose cell type was incorrectly predicted for the Baron single cell pancreas data for different approaches,
- preprocessing: +1, norm, scale -> 830 out of 8569 cells
- preprocessing: +1, norm, log, scale -> 886 out of 8569 cells
- preprocessing: +1, norm, scale -> 886 out of 8569 cells
- preprocessing: +1, log, +1, norm, scale -> 962 out of 8569 cells
- preprocessing: +1, norm, scale, using top 1100 differentially expressed genes per cell type group -> 978 out of 8569 cells
- preprocessing: +1, log, norm, scale, using top 1100 differentially expressed genes per cell type group  -> 1221 out of 8569 cells
- preprocessing: +1, log, norm, scale -> 2001 out of 8569 cells
- preprocessing: norm, scale (no +1 before hand) -> 4059 out of 8569 cells

## Credits

List your collaborators, if any, with links to their GitHub profiles.

If you used any third-party assets that require attribution, list the creators with links to their primary web presence in this section.


anndata: Annotated data
Isaac Virshup, Sergei Rybakov, Fabian J. Theis, Philipp Angerer, F. Alexander Wolf
bioRxiv 2021 Dec 19. doi: 10.1101/2021.12.16.473007.

## License

The last section of a high-quality README file is the license. This lets other developers know what they can and cannot do with your project. If you need help choosing a license, refer to [https://choosealicense.com/](https://choosealicense.com/).






