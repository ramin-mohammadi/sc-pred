import numpy as np
import pandas as pd
from sklearn import preprocessing

# uses numpy and sklearn to perform natural log transformation, normalize, and scale the data
# IMPORTANT: assay_df must be a pandas data frame where rows are cells (samples) and columns are genes (features)
def preprocess(assay_df: pd.DataFrame, order=1, log=True, normalize=True, scale=True, add1=False):
    print("Beginning to preprocess data...\n")
    print("Before preprocessing:\n", assay_df)
    
    # default preprocess order: normalize, natural log transform, scale
    if order == 1: 
        assay_df = norm(assay_df, normalize)
        assay_df = logTransform(assay_df, log)
        assay_df = scaleData(assay_df, scale)
  
    # natural log transform, normalize, scale
    else:
        assay_df = logTransform(assay_df, log)
        if add1:
            assay_df += 1 # may improve results
        assay_df = norm(assay_df, normalize)
        assay_df = scaleData(assay_df, scale)
  
    print("After preprocessing:\n", assay_df)
    return assay_df

# Normalize (represent values between 0 and 1). 
# Normalizes across axis 1 (by defualt), so across values in a single row (a row here should represent a cell so a cell's gene counts are normalized among all of the genes)    
# The default norm for normalize() is L2, also known as the Euclidean norm. The L2 norm formula is the square root of the sum of the squares of each value.
def norm(assay_df: pd.DataFrame, normalize: bool):
    if normalize:
        assay_df = preprocessing.normalize(assay_df)
        print("Normalized:\n", assay_df)
    return assay_df

# Natural log (base e) transformation. Helps reduce range of values
def logTransform(assay_df: pd.DataFrame, log: bool):
    if log:
        assay_df = np.log(assay_df)
        print("Natural log transformation:\n", assay_df)
    return assay_df

# Scaling Data (The purpose is to ensure that all features contribute equally to the model and to avoid the domination of features with larger values.)
def scaleData(assay_df: pd.DataFrame, scale: bool):
    if scale:
        scaler = preprocessing.StandardScaler().fit(assay_df)
        assay_df = scaler.transform(assay_df)
        print("Scaled:\n", assay_df)
    return assay_df