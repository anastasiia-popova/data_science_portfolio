# helper functions for the project 
# author: Anastasia Popova ppva.nastya@proton.me
import os

import numpy as np
import pandas as pd
import anndata as ad
import pyroe
import scanpy as sc

from scipy.stats import median_abs_deviation
from scipy.stats import spearmanr

import matplotlib.pyplot as plt 

# import mudata as md
# from mudata import MuData
# import muon as mu

def return_path(filename='Analysis_of_CITE-seq_data.ipynb'):
    """Uses os to return the correct path of a directory.
    Default is a directory containing Analysis_of_CITE-seq_data.ipynb file. 
    """
    absolute_path = os.path.abspath(filename)
    directory_name = os.path.dirname(absolute_path).replace("\\" , "/" )
    return directory_name 

def gene_names(gene_id_path, adata): 
    
    """
    Map gene IDs to gene names and return a DataFrame.

    Parameters:
    ----------
    gene_id_path : str
        Path to a directory containing the 'gene_id_to_name.tsv' file with gene ID to gene name mappings.
        
    adata : AnnData
        Anndata object containing the gene expression data.

    Returns:
    -------
    pd.DataFrame
        A DataFrame with gene names as the index and relevant gene information.
        
    Description:
    -----------
    This function takes a path to a directory containing a 'gene_id_to_name.tsv' file, which should have two columns:
    'gene_id' and 'gene_name', with tab as the delimiter. It also takes an AnnData object ('adata') containing gene 
    expression data. The function maps gene IDs in 'adata' to gene names using the provided mapping file and returns a 
    DataFrame with gene names as the index. The returned DataFrame contains relevant gene information from 'adata'.
    """
    
    gene_id_to_name = np.genfromtxt(f"{gene_id_path}/gene_id_to_name.tsv", delimiter='\t', dtype=str )

    gene_id_to_name = [tuple(g) for g in gene_id_to_name]
    gene_dict  = dict(set(gene_id_to_name))

    df_genes = adata.var.copy()
    df_genes  = df_genes.reset_index()
    df_genes["gene_name"] = df_genes["gene_ids"].map(gene_dict)
    #df_genes = df_genes.set_index("gene_name")
    
    return df_genes 

def gene_names_coding_proteins(gene_id_path, adata): 
    
    """
    Retrieve gene names and information for protein-coding genes from an AnnData object.

    This function reads a gene ID file (gene_types.txt) to extract
    information about protein-coding genes and their associated gene names. It then filters
    the provided AnnData object to include only protein-coding genes and returns a DataFrame
    containing their information.

    Parameters:
    - gene_id_path (str): The path to the directory containing the gene ID file (gene_types.txt).
    - adata (AnnData): An AnnData object representing gene expression data.

    Returns:
    - pd.DataFrame: A DataFrame containing ids and names for protein-coding genes, indexed by gene name.

    Note:
    - The 'gene_types.txt' file should have columns 'Gene stable ID version', 'Gene name', and 'Gene type'
      separated by tabs. This function assumes that 'Gene type' is used to identify protein-coding genes.
    - If 'Gene name' is missing for a gene, it will be replaced with the corresponding 'Gene stable ID'.
    - The function filters genes in the AnnData object to match the protein-coding genes found in the gene ID file.
    
    """
    
    gene_id_and_type = pd.read_csv(f"{gene_id_path}/gene_types.txt", sep='\t')

    gene_id_protein = gene_id_and_type[gene_id_and_type["Gene type"]== "protein_coding"][["Gene stable ID version", "Gene name", "Gene type" ]]
    unstable_genes = list(gene_id_protein[gene_id_protein["Gene name"].isna()].index)

    for index in unstable_genes:
        gene_id_protein.loc[index, "Gene name"] = gene_id_and_type.loc[index, "Gene stable ID"]

    gene_id_protein.set_index("Gene stable ID version", inplace=True)
    gene_id_protein.index.name = "gene_ids"

    gene_id_protein.rename(columns={"Gene name": "gene_name", "Gene type":"gene_type"}, inplace=True)
    gene_id_protein_dict = gene_id_protein.drop(columns=["gene_type"]).to_dict()

    
    df_genes = adata.var.copy()
    #df_genes  = df_genes.reset_index()    
    df_genes = df_genes[df_genes["gene_ids"].isin(list(gene_id_protein.index)) ]
    df_genes["gene_name"] = df_genes["gene_ids"].map(gene_id_protein_dict["gene_name"])
    
    #df_genes.set_index("gene_name", inplace=True)
    
    return df_genes


def is_outlier(adata, metric: str, nmads: int):
    """
    Detect outliers in a single metric using the median and median absolute deviation (MAD).

    Parameters:
    ----------
    adata : AnnData
        Anndata object containing the data to be analyzed.

    metric : str
        Name of the metric (e.g., gene expression) to analyze for outliers.

    nmads : int
        Number of median absolute deviations (MADs) to consider as a threshold for outliers.

    Returns:
    -------
    pd.Series
        A boolean Series indicating whether each data point in the specified metric is an outlier.

    Description:
    -----------
    This function detects outliers in a given metric of an AnnData object ('adata') using the median and median 
    absolute deviation (MAD). Outliers are determined based on whether the metric values fall outside a certain 
    range defined by the median ± (nmads * MAD). The function returns a boolean Series where True indicates 
    data points that are considered outliers.

    The median and MAD are robust statistical measures, making this method suitable for identifying outliers 
    in metrics that may have non-normal distributions or contain noise.
    
    For more details on this method, refer to the original source:
    https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html
    """
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier

## This method is from doubletdetection Python package
def convergence(clf, show=False, save=None, p_thresh=1e-7, voter_thresh=0.9):
    """Produce a plot showing number of cells called doublet per iter

    Args:
        clf (BoostClassifier object): Fitted classifier
        show (bool, optional): If True, runs plt.show()
        save (str, optional): filename for saved figure,
            figure not saved by default
        p_thresh (float, optional): hypergeometric test p-value threshold
            that determines per iteration doublet calls
        voter_thresh (float, optional): fraction of iterations a cell must
            be called a doublet

    Returns:
        matplotlib figure
    """
    log_p_thresh = np.log(p_thresh)
    doubs_per_run = []
    # Ignore numpy complaining about np.nan comparisons
    with np.errstate(invalid="ignore"):
        for i in range(clf.n_iters):
            cum_log_p_values = clf.all_log_p_values_[: i + 1]
            cum_vote_average = np.mean(
                np.ma.masked_invalid(cum_log_p_values) <= log_p_thresh, axis=0
            )
            cum_doublets = np.ma.filled((cum_vote_average >= voter_thresh).astype(float), np.nan)
            doubs_per_run.append(np.nansum(cum_doublets))

    # Ignore warning for convergence plot
    # with warnings.catch_warnings():
    #     warnings.filterwarnings(action="ignore", module="matplotlib", message="^tight_layout")

    f, ax = plt.subplots(1, 1, figsize=(4, 4), dpi=150)
    ax.plot(np.arange(len(doubs_per_run)), doubs_per_run)
    ax.set_xlabel("Number of Iterations")
    ax.set_ylabel("Number of Predicted Doublets")
    ax.set_title("Predicted Doublets per Iteration")

    if show is True:
        plt.show()
    if isinstance(save, str):
        f.savefig(save, format="png", bbox_inches="tight")

    return f

def intersection(lst1, lst2):
    """
    Find the intersection of two lists.

    Parameters:
    ----------
    lst1 : list
        The first list to be compared.

    lst2 : list
        The second list to be compared.

    Returns:
    -------
    list
        A new list containing elements that are common to both input lists (intersection).

    Description:
    -----------
    This function takes two input lists, 'lst1' and 'lst2', and returns a new list containing elements that are present 
    in both 'lst1' and 'lst2'. It computes the intersection of the two lists, effectively finding elements that are 
    common to both lists.

    Example:
    --------
    >>> list1 = [1, 2, 3, 4, 5]
    >>> list2 = [3, 4, 5, 6, 7]
    >>> result = intersection(list1, list2)
    >>> print(result)
    [3, 4, 5]
    """
    
    lst3 = [value for value in lst1 if value in lst2]
    return lst3

def spearman_corr(gene_and_ant_expr):
    
    """
    Compute Spearman rank correlation coefficients between gene expression and antibody measurements.

    Parameters:
    ----------
    gene_and_ant_expr : pd.DataFrame
        A DataFrame containing gene expression and antibody measurement data.

    Returns:
    -------
    pd.DataFrame
        A DataFrame where rows represent genes, columns represent antibodies, and values are Spearman 
        correlation coefficients.

    Description:
    -----------
    This function calculates Spearman correlation coefficients between gene expression levels and antibody 
    measurements for a given dataset. The input 'gene_and_ant_expr' should be a DataFrame where columns 
    correspond to different genes and antibodies. The function computes the correlation coefficient for each gene 
    with respect to four specific antibodies: IgM, IgG-Fc, Ig-kappa, and Ig-lambda.

    The resulting DataFrame provides a correlation matrix where rows represent genes, columns represent antibodies, 
    and the cell values indicate the Spearman correlation coefficients between the corresponding gene and antibody 
    expression profiles.

    Example:
    --------
    >>> gene_antibody_data = pd.DataFrame(...)
    >>> result = spearman_corr(gene_antibody_data)
    >>> print(result)
                 IgM   IgG-Fc  Ig-kappa  Ig-lambda
    Gene1    0.789   0.654    0.543     0.721
    Gene2    0.632   0.543    0.812     0.456
    Gene3    0.721   0.654    0.789     0.543

    Note:
    -----
    This function uses Spearman's rank correlation, which is a non-parametric measure of correlation that assesses 
    the strength and direction of monotonic relationships between variables.
    """
    
    all_names = list(gene_and_ant_expr.columns)
    ant_names = ["IgM","IgG-Fc","Ig-kappa","Ig-lambda"]
    gene_names = [x for x in all_names if x not in ant_names]
    spearman_corr = dict()


    for g in gene_names:
        spearman_corr[g] = {}
        for a in ant_names:
            spearman_corr_coef, _ = spearmanr(gene_and_ant_expr[g], gene_and_ant_expr[a])
            spearman_corr[g][a] = spearman_corr_coef
            
    return pd.DataFrame(spearman_corr)


def unique(list_):
    """
    Return a list containing unique elements from the input list, preserving their order.

    Parameters:
    ----------
    list_ : list
        The input list from which unique elements will be extracted.

    Returns:
    -------
    list
        A new list containing unique elements from the input list, in the order of their first occurrence.

    Description:
    -----------
    This function takes an input list 'list_' and returns a new list containing only the unique elements found in 
    'list_'. The order of elements in the output list is preserved, meaning they appear in the same order as 
    their first occurrence in the input list.

    Example:
    --------
    >>> input_list = [3, 1, 2, 3, 2, 4, 5, 1]
    >>> result = unique(input_list)
    >>> print(result)
    [3, 1, 2, 4, 5]

    Note:
    -----
    This function is useful when you need to extract and retain only the unique elements from a list while 
    maintaining the original order.
    """
    
    unique_list = []
    for x in list_:
        if x not in unique_list:
            unique_list.append(x)
            
    return unique_list
        
def top_correlated_genes(df_corr, antibody, n):
    """
    Retrieve the top 'n' genes most positively correlated with a specific antibody measurement.

    Parameters:
    ----------
    df_corr : pd.DataFrame
        A DataFrame containing correlation values between antibodies and genes.

    antibody : str
        The specific antibody for which you want to find correlated genes.

    n : int
        The number of top positively correlated genes to retrieve.

    Returns:
    -------
    pd.DataFrame
        A DataFrame with the 'n' genes most strongly correlated with the specified antibody.

    Description:
    -----------
    This function takes a DataFrame 'df_corr' containing correlation values between antibodies and genes and 
    identifies the top 'n' genes that are most strongly correlated with a specified 'antibody'. The function 
    sorts the genes in descending order based on their correlation with the antibody and returns a DataFrame 
    containing the top correlated genes along with their correlation values.
    """
    df = pd.DataFrame(df_corr.loc[antibody].sort_values(ascending=False)[:n]) #.sort_values(key=abs, ascending=False)[:n])
    return df

def antibody_corr_matrix(
    antibody_name,
    spearman_corr_1, spearman_corr_2, spearman_corr_3, spearman_corr_4,
    special_donor="Donor-2",
    n_top=10):
    
    """
    Generate a correlation matrix of top correlated genes with an antibody for multiple donors.

    Parameters:
    ----------
    antibody_name : str
        The name of the antibody for which correlation is computed.

    spearman_corr_1, spearman_corr_2, spearman_corr_3, spearman_corr_4 : pd.DataFrame
        DataFrames containing Spearman correlation values between genes and antibodies for four different donors.

    special_donor : str, optional
        The specific donor to focus on (e.g., "Donor-1", "Donor-2", "Donor-3", "Donor-4"). Default is "Donor-2".

    n_top : int, optional
        The number of top correlated genes to include in the correlation matrix. Default is 10.

    Returns:
    -------
    pd.DataFrame
        A DataFrame representing the correlation matrix for the specified antibody and donor(s).

    Description:
    -----------
    This function calculates the correlation between a given antibody and genes for multiple donors. It extracts the top 
    'n_top' genes that are most correlated with the specified 'antibody_name' for each donor separately. The resulting 
    correlation matrix includes columns for each donor and rows for the top correlated genes.

    You can specify a 'special_donor' to focus on, and the function will provide correlation information only for that 
    donor. If 'special_donor' is not specified, correlations for all donors will be included in the matrix.

    Example:
    --------
    >>> antibody_name = "IgM"
    >>> donor1_corr = pd.DataFrame(...)
    >>> donor2_corr = pd.DataFrame(...)
    >>> donor3_corr = pd.DataFrame(...)
    >>> donor4_corr = pd.DataFrame(...)
    >>> result = antibody_corr_matrix(antibody_name, donor1_corr, donor2_corr, donor3_corr, donor4_corr)
    >>> print(result)
                     Donor-1  Donor-2  Donor-3  Donor-4
    Gene1_correlation    0.789    0.654    0.543    0.721
    Gene2_correlation    0.632    0.543    0.812    0.456
    ...

    Note:
    -----
    The function extracts the top correlated genes based on the specified 'n_top' value and returns a correlation matrix 
    for the selected donor(s) and gene-antibody pairs.
    """
    
    if "Donor-1" in spearman_corr_1.columns: 
        filtered_corr_1 = spearman_corr_1.drop(columns=["Donor-1"])
        filtered_corr_2 = spearman_corr_2.drop(columns=["Donor-2"])
        filtered_corr_3 = spearman_corr_3.drop(columns=["Donor-3"])
        filtered_corr_4 = spearman_corr_4.drop(columns=["Donor-4"])
        
    else:
        filtered_corr_1 = spearman_corr_1.copy()
        filtered_corr_2 = spearman_corr_2.copy()
        filtered_corr_3 = spearman_corr_3.copy()
        filtered_corr_4 = spearman_corr_4.copy()
        
    
    Ig_top_1 = top_correlated_genes(filtered_corr_1 , antibody_name, n = n_top)
    Ig_top_2 = top_correlated_genes(filtered_corr_2 , antibody_name, n = n_top)
    Ig_top_3 = top_correlated_genes(filtered_corr_3 , antibody_name, n = n_top)
    Ig_top_4 = top_correlated_genes(filtered_corr_4 , antibody_name, n = n_top)
    
    #top_genes_Ig_all_donors =  Ig_top_2.index.to_list().copy() + Ig_top_1.index.to_list()  + Ig_top_3.index.to_list() + Ig_top_4.index.to_list()
    #top_genes_Ig_all_donors = unique(Ig_top_1.index.to_list() + Ig_top_2.index.to_list() + Ig_top_3.index.to_list() + Ig_top_4.index.to_list())
        
    if special_donor=="Donor-1":
        top_genes_Ig_all_donors =  Ig_top_1.index.to_list().copy() 
    elif special_donor=="Donor-2":
        top_genes_Ig_all_donors =  Ig_top_2.index.to_list().copy() 
    elif special_donor=="Donor-3":
        top_genes_Ig_all_donors =  Ig_top_3.index.to_list().copy() 
    elif special_donor=="Donor-4":
        top_genes_Ig_all_donors =  Ig_top_3.index.to_list().copy() 
    elif special_donor==None:
        top_genes_Ig_all_donors =  Ig_top_2.index.to_list().copy() + Ig_top_1.index.to_list()  + Ig_top_3.index.to_list() + Ig_top_4.index.to_list()
    else:
        print("The `special_donor` you entered is invalid. Please try again")
        raise KeyError('Donor name must be Donor-1, Donor-2, Donor-3, Donor-4 or None.')

    Ig_d1 = filtered_corr_1.loc[antibody_name][top_genes_Ig_all_donors].to_list()
    Ig_d2 = filtered_corr_2.loc[antibody_name][top_genes_Ig_all_donors].to_list()
    Ig_d3 = filtered_corr_3.loc[antibody_name][top_genes_Ig_all_donors].to_list()
    Ig_d4 = filtered_corr_4.loc[antibody_name][top_genes_Ig_all_donors].to_list()

    top_Ig_all = pd.DataFrame({"Donor-1": Ig_d1, "Donor-2": Ig_d2, "Donor-3": Ig_d3, "Donor-4": Ig_d4}, index=top_genes_Ig_all_donors)
    top_Ig_all.index.name = f"{antibody_name} top genes" 
    
    return top_Ig_all