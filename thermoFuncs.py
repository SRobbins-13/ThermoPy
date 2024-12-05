# -*- coding: utf-8 -*-
"""
Created on September 3, 2023

@authors: Samuel Robbins, Chelsea Mackaman-Lofland

"""

###############################################################
# Import required modules and notebook setup
###############################################################

# Data Handling/Visualization Modules
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.font_manager as fm
import matplotlib.colors as mcolors
import matplotlib.lines as mlines
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.axes import Axes

import seaborn as sns
import random

# Statistics Modules
import numpy as np
import statsmodels.api as sm
from scipy.stats import t

# Path Modules
import pathlib
from pathlib import Path

# Typing Modules
from typing import List, Tuple, Optional, Union

sns.set(style='white')

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# Shades of Grey
GREY10 = "#1a1a1a"
GREY30 = "#4d4d4d"
GREY40 = "#666666"
GREY50 = "#7f7f7f"
GREY60 = "#999999"
GREY75 = "#bfbfbf"
GREY91 = "#e8e8e8"
GREY98 = "#fafafa"

#######################################################################
# Functions for loading and saving data from/to Excel 
#######################################################################
def loadDataExcel(file_path: str, samples: str = 'Samples', aliquots: str = 'Aliquots', ID_col:str = 'Sample'
    )-> Tuple[pd.DataFrame, List[Tuple[str, str]], List[str], pd.DataFrame]:
    """
    Loads an Excel file containing thermochronologic data.

    Parameters
    ----------
    file_path : str
        Full file path of the data sheet or name of data sheet if in the same folder.
    samples : str, optional
        The name of the Excel worksheet that contains sample information. The default name is 'Samples'.
    aliquots : str, optional
        The name of the Excel worksheet that contains grain analysis information. The default name is 'Aliquots'.
    ID_col : str, optional
        The name of the column that contains unique sample identifiers. Default is 'Sample'.

    Returns
    -------
    Tuple[pd.DataFrame, List[Tuple[str, str]], List[str], pd.DataFrame]
        samples : pd.DataFrame
            Database of sample information indexed by sample name.
        sample_list : List[Tuple[str, str]]
            List of all samples in the database.
        transect_list : List[str]
            List of all transects for later plotting.
        aliquots : pd.DataFrame
            Database of aliquot data (many aliquots: 1 sample).
    """
    samples = pd.read_excel(file_path, sheet_name = samples)
    samples.set_index(ID_col, inplace = True, drop = False)
    sample_list = list(zip(samples.Sample, samples.Mineral))

    transect_list = list(set([item for item in samples.Transect.to_list() if not(pd.isnull(item)) == True]))

    aliquots = pd.read_excel(file_path, sheet_name = aliquots)
    aliquots.set_index('Aliquot', inplace = True, drop = False)

    return samples, sample_list, transect_list, aliquots

def writeToExcel(sample_stats: pd.DataFrame, aliquots: pd.DataFrame, filename: str, folder: str = 'Data_Sheets') -> None:
    """
    Writes the sample statistics and aliquot classification dataframes to an Excel file. This will link the keep/reject 
    designation with the relevant sample means, std. dev, etc. 

    Parameters
    ----------
    sample_stats : pd.DataFrame
        Database of sample statistics indexed by sample name.
    aliquots : pd.DataFrame
        Database of aliquot data with outlier classification indexed by aliquot name (many aliquots: 1 sample).
    filename : str
        File name to write to.
    folder : str, optional
        Subfolder to write to (default is 'Data_Sheets').

    Returns
    -------
    None
        Outputs dataframes to an Excel file with two sheets ('Samples' and 'Aliquots').
    """
    pathlib.Path(folder).mkdir(parents=True, exist_ok=True)
    filepath = '{}/{}'.format(folder, filename)
    
    
    # Create a Pandas Excel writer using XlsxWriter as the engine.
    writer = pd.ExcelWriter(filepath, engine='xlsxwriter')

    # Write each dataframe to a different worksheet.
    sample_stats.to_excel(writer, sheet_name='Samples',index = False)
    aliquots.to_excel(writer, sheet_name='Aliquots',index = False)

    # Close the Pandas Excel writer and output the Excel file.
    writer.close()

#######################################################################
# Functions for summary statistics of thermochronologic data 
#######################################################################
def determineOutliersIQR(samples: pd.DataFrame, aliquots: pd.DataFrame, sample_list: List[str])-> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Determines the quantiles and outlier bounds for each sample in the sample_list and returns updated sample and
    aliquot databases.

    Parameters
    ----------
    samples : pd.DataFrame
        Database of sample information indexed by sample name.
    aliquots : pd.DataFrame
        Database of aliquot data (many aliquots: 1 sample).
    sample_list: List[Tuple[str, str]]
        List of all samples in the database, each entry containing a sample identifier and a sample type.

    Returns
    -------
    samples_iqr : pd.DataFrame
        Database of sample information indexed by sample name with quartile information specified.
    aliquots_iqr : pd.DataFrame
        Database of aliquot information indexed by aliquot name with outlier specification ('keep' or 'reject').
    """
    pd.options.mode.chained_assignment = None

    samples_iqr = samples.copy()
    aliquots_iqr = aliquots.copy()

    for sample in sample_list:
        # fission track samples only have 1 "aliquot" per sample
        if sample[1].upper() == 'AFT':
            # Get all aliquots associated with a single sample
            sample_df = aliquots[aliquots.Sample == sample[0]]

            N = int(sample_df.shape[0])
            aliquot_list = sample_df.index.to_list()
            
            for grain in aliquot_list:
                aliquots_iqr.loc[[grain], 'outlier'] = 'keep'

        else: 
            # Get all aliquots associated with a single sample
            sample_df = aliquots[aliquots.Sample == sample[0]]

            N = int(sample_df.shape[0])
            aliquot_list = sample_df.index.to_list()

            # Calculate quantiles
            Q1 = sample_df.Corrected_Date_Ma.quantile(0.25)
            Q3 = sample_df.Corrected_Date_Ma.quantile(0.75)
            IQR = Q3-Q1

            # Determine Lower/Upper Bounds
            lowerB = Q1 - 1.5 * IQR
            upperB = Q3 + 1.5 * IQR

            # Add values to sample database
            samples_iqr.loc[[sample[0]],'N_grains'] = N
            samples_iqr.loc[[sample[0]],'Q1'] = Q1
            samples_iqr.loc[[sample[0]],'Q3'] = Q3
            samples_iqr.loc[[sample[0]],'IQR'] = IQR
            samples_iqr.loc[[sample[0]],'Lower_Bound'] = lowerB
            samples_iqr.loc[[sample[0]],'Upper_Bound'] = upperB

            # Determine outliers based on bounds
            for grain in aliquot_list:
                grain_date = sample_df.loc[[grain],'Corrected_Date_Ma'].iloc[0] 

                if grain_date < lowerB or grain_date > upperB:
                    aliquots_iqr.loc[[grain], 'outlier'] = 'reject'
                else:
                    aliquots_iqr.loc[[grain], 'outlier'] = 'keep'

    samples_iqr['N_grains'] = samples_iqr['N_grains'].astype('Int64')

    return samples_iqr, aliquots_iqr

def chauvenetsCriterionOutliers(samples: pd.DataFrame, aliquots: pd.DataFrame, sample_list: List[str])-> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Identifies outliers based on Chauvenet's Criterion. Method creates an acceptable band of data around the mean for each sample.
    Any individual aliquot age that falls outside its sample bounds is marked as an outlier.

    Parameters
    ----------
    samples : pd.DataFrame
        Database of sample information indexed by sample name.
    aliquots : pd.DataFrame
        Database of aliquot data (many aliquots: 1 sample).
    sample_list: List[Tuple[str, str]]
        List of all samples in the database, each entry containing a sample identifier and a sample type.


    Returns
    -------
    sample_chauv : pd.DataFrame
        Database of sample information indexed by sample name with Chauvenet Criterion specified.
    aliquots_chauv : pd.DataFrame
        Database of aliquot information indexed by aliquot name with outlier specification ('keep' or 'reject').
    """

    pd.options.mode.chained_assignment = None
    
    # Define the Chauvenet's Criterion Table
    n = [2,3,4,5,6,7,8,9,10,11,12,13,14,15]
    T_critical = [1.150,1.383,1.534,1.645,1.732,1.803,1.863,1.915,1.960,2.000,2.037,2.070,2.100,2.128]

    criterion_table = pd.DataFrame(list(zip(n,T_critical)),columns = ['n','T_critical'], index = n)

    # Create copies of samples and aliquots DFs so as to not alter original dataframes
    sample_chauv = samples.copy()
    aliquots_chauv = aliquots.copy()

    for sample in sample_list:

        sample_df = aliquots[aliquots.Sample == sample[0]]
        aliquot_list = sample_df.index.to_list()

        N = int(sample_df.shape[0])

        if N==0:
            print('{} has no associated aliquots in the data sheet'.format(sample))
            pass

        elif N == 1:
            sample_chauv.loc[[sample[0]],'N_grains'] = N
            for grain in aliquot_list:
                aliquots_chauv.loc[[grain],'outlier'] = 'keep'

        else: 
            T_critical = criterion_table.T_critical.loc[[N]].iloc[0]

            sample_chauv.loc[[sample[0]],'N_grains'] = N
            sample_chauv.loc[[sample[0]],'T_critical'] = T_critical

            mean = sample_df.Corrected_Date_Ma.mean()
            std = sample_df.Corrected_Date_Ma.std()

            for grain in aliquot_list:
                grain_date = sample_df.loc[[grain],'Corrected_Date_Ma'].iloc[0] 

                T = abs(grain_date - mean) / std

                aliquots_chauv.loc[[grain],'T'] = T
                aliquots_chauv.loc[[grain],'T_critical'] = T_critical

                if T > T_critical:
                    aliquots_chauv.loc[[grain],'outlier'] = 'reject'
                else:
                    aliquots_chauv.loc[[grain],'outlier'] = 'keep'



    sample_chauv['N_grains'] = sample_chauv['N_grains'].astype('Int64')

    return sample_chauv, aliquots_chauv

def overwriteOutlierClassification(outliers_to_keep: List[str], outliers_to_reject: List[str], aliqouts_db: pd.DataFrame)-> pd.DataFrame:
    """
    Overwrites the outlier classification of specific grains based on user specifications (i.e., changing from 'reject' 
    to 'keep' or vice versa).

    Parameters
    ----------
    outliers_to_keep : List[str]
        List of aliquots to change from 'reject' to 'keep'.
    outliers_to_reject : List[str]
        List of aliquots to change from 'keep' to 'reject'.
    aliqouts_db : pd.DataFrame
        Database of aliquot data (many aliquots: 1 sample) that have already been assigned an outlier classification.

    Returns
    -------
    pd.DataFrame
        Updated aliquot database with user-defined outlier specifications.
    """
    
    # Create copies of samples and aliquots DFs so as to not alter original dataframes
    aliquots_outliers = aliqouts_db.copy()
    
    # changing from 'reject' to 'keep'
    for grain in outliers_to_keep:
        aliquots_outliers.loc[grain,'outlier'] = 'keep'
    
    # changing from 'keep' to 'reject'
    for grain in outliers_to_reject:
        aliquots_outliers.loc[grain,'outlier'] = 'reject'
    
    return aliquots_outliers

def viewOutliers(samples: pd.DataFrame, aliquots: pd.DataFrame, sample_list: List[str], stat_scheme: str)-> pd.DataFrame:
    """
    Views outliers in the aliquot data based on the specified statistical scheme ('IQR' or 'Chauvenet').

    Parameters
    ----------
    samples : pd.DataFrame
        Database of sample information indexed by sample name.
    aliquots : pd.DataFrame
        Database of aliquot data (many aliquots: 1 sample).
    sample_list : List[str]
        List of all samples in the database.
    stat_scheme : str
        Statistical scheme to use for outlier detection ('IQR' or 'Chauvenet').

    Returns
    -------
    pd.DataFrame
        DataFrame of aliquots classified as outliers ('reject') according to the chosen scheme.
    """

    if stat_scheme == 'IQR':
        samples_iqr, aliquots_iqr = determineOutliersIQR(samples, aliquots, sample_list)

        return aliquots_iqr[aliquots_iqr.outlier=='reject']

    if stat_scheme == 'Chauvenet':
        samples_chauv, aliquots_chauv = chauvenetsCriterionOutliers(samples, aliquots, sample_list)

        return aliquots_chauv[aliquots_chauv.outlier=='reject']

def calculateFullSummaryStats(samples: pd.DataFrame, aliquots: pd.DataFrame, sample_list: List[str], stat_scheme: str, 
    keep_all_outliers: bool = False, aliquots_to_keep: List[str] = None, aliquots_to_reject: List[str] = None, 
    save_data_excel: bool = False, filename: str ='summary_statistics_export.xlsx') -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Calculates summary statistics for samples and aliquots, allowing for outlier detection via IQR or Chauvenet's Criterion.
    
    Parameters
    ----------
    samples : pd.DataFrame
        Database of sample information indexed by sample name.
    aliquots : pd.DataFrame
        Database of aliquot data (many aliquots: 1 sample).
    sample_list : List[str]
        List of all samples in the database.
    stat_scheme : str
        Statistical scheme to use for outlier detection ('IQR' or 'Chauvenet').
    keep_all_outliers : bool, optional
        If True, all outliers will be retained (default is False).
    aliquots_to_keep : List[str], optional
        List of aliquots to retain despite being marked as outliers.
    aliquots_to_reject : List[str], optional
        List of aliquots to reject despite not being marked as outliers.
    save_data_excel : bool, optional
        If True, saves the resulting statistics to an Excel file (default is False).
    filename : str, optional
        The filename for saving if `save_data_excel` is True (default is 'summary_statistics_export.xlsx').

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        The resulting samples and aliquots statistics dataframes.
    """
        
    def process_outliers(sample_df, aliquot_df, outliers_to_keep, outliers_to_reject):
        aliquots_stats = overwriteOutlierClassification(outliers_to_keep, outliers_to_reject, aliquot_df)
        samples_stats = calculateUnweightedSummaryStats(sample_df, aliquots_stats, sample_list)
        samples_stats = calculateWeightedInverseVarianceSummaryStats(samples_stats, aliquots_stats, sample_list)
        samples_stats = calculateWeightedRelativeDeviationSummaryStats(samples_stats, aliquots_stats, sample_list)
        return samples_stats, aliquots_stats

    # Define default method and handle each scheme
    if stat_scheme == 'IQR':
        sample_df, aliquot_df = determineOutliersIQR(samples, aliquots, sample_list)
    elif stat_scheme == 'Chauvenet':
        sample_df, aliquot_df = chauvenetsCriterionOutliers(samples, aliquots, sample_list)
    else:
        raise ValueError("stat_scheme must be either 'IQR' or 'Chauvenet'.")

    # Set up outlier handling
    if keep_all_outliers:
        outliers_to_keep = aliquot_df[aliquot_df.outlier == 'reject']['Aliquot'].tolist()
        outliers_to_reject = []
    else:
        outliers_to_keep = aliquots_to_keep if aliquots_to_keep else []
        outliers_to_reject = aliquots_to_reject if aliquots_to_reject else []

    # Calculate stats and save if requested
    samples_stats, aliquots_stats = process_outliers(sample_df, aliquot_df, outliers_to_keep, outliers_to_reject)
    if save_data_excel:
        writeToExcel(samples_stats, aliquots_stats, filename, folder='Summary_Statistics')

    return samples_stats, aliquots_stats

def calculateUnweightedSummaryStats(samples: pd.DataFrame, aliquots: pd.DataFrame, sample_list: List[str])-> pd.DataFrame:
    """
    Calculates unweighted summary statistics for each sample, including mean, standard deviation, standard error based
    on standard deviation and individual uncertainty, and MSWD.

    Parameters
    ----------
    samples : pd.DataFrame
        DataFrame containing sample information indexed by sample name.
    aliquots : pd.DataFrame
        DataFrame containing aliquot data indexed by aliquot name, where multiple aliquots may correspond to a single sample.
    sample_list : List[str]
        List of sample identifiers to be analyzed.

    Returns
    -------
    sample_update : pd.DataFrame
        Updated sample DataFrame with additional columns: 
        - 'Mean' : Unweighted mean of the sample ages.
        - 'StDev' : Unweighted standard deviation of the sample ages.
        - 'SEsd' : Standard error of the Unweighted standard deviation.
        - 'SEiu' : Standard error based on individual uncertainties.
        - 'MSWD' : Mean Square Weighted Deviation, representing goodness of fit.

    Notes
    -----
    For samples identified as 'AFT' or 'ZFT', the function uses the mean and standard deviation directly from the 
    'Corrected_Date_Ma' and 'Corrected_Uncertainty_1σ_Ma' columns without further calculation. 
    """
    sample_update = samples.copy()

    for sample in sample_list:
        # Check if the sample is of type 'AFT' or 'ZFT' (fission track)
        if sample[1].upper() in {'AFT', 'ZFT'}:
            mean = aliquots.loc[sample[0], 'Corrected_Date_Ma']
            stdev = aliquots.loc[sample[0], 'Corrected_Uncertainty_1σ_Ma']

            sample_update.loc[[sample[0]],'Mean'] = mean
            sample_update.loc[[sample[0]],'StDev'] = stdev
        else:
            # Filter aliquots for the sample and keep only non-outliers
            sample_df = aliquots[(aliquots.Sample == sample[0]) & (aliquots.outlier == 'keep')]
            ages = [x for x in sample_df.Corrected_Date_Ma]
            errors = [x for x in sample_df.Corrected_Uncertainty_1σ_Ma]

            n = len(sample_df)
            
            if n > 0:
                # Mean (x_bar)
                mean = sum(ages)/len(ages)

                # Standard Deviation (SD)
                stdev = errors[0] if n == 1 else np.sqrt(sum((x - mean) ** 2 for x in ages) / (n - 1))

                # Standard Error based on SD (SEsd)
                se_sd = stdev / np.sqrt(n)

                # Standard Error based on individual uncertainty (SEiu)
                se_iu = np.sqrt(sum(x**2 for x in errors)) / n

                # Mean Square Weighted Deviation (MSWD)
                mswd = None if n <= 1 else sum((x - mean) ** 2 / y ** 2 for x, y in zip(ages, errors)) / (n - 1)


                # Add values to sample database
                sample_update.loc[[sample[0]],'Mean'] = mean
                sample_update.loc[[sample[0]],'StDev'] = stdev
                sample_update.loc[[sample[0]],'SEsd'] = se_sd
                sample_update.loc[[sample[0]],'SEiu'] = se_iu
                sample_update.loc[[sample[0]],'MSWD'] = mswd

    return sample_update

def calculateWeightedInverseVarianceSummaryStats(samples: pd.DataFrame, aliquots: pd.DataFrame, sample_list: List[str])-> pd.DataFrame:
    """
    Calculates the summary statistics weighted by inverse variance for each sample, including weighted mean, standard deviation, and errors, and MSWD

    Parameters
    ----------
    samples : pd.DataFrame
        DataFrame containing sample information indexed by sample name.
    aliquots : pd.DataFrame
        DataFrame containing aliquot data indexed by aliquot name, where multiple aliquots may correspond to a single sample.
    sample_list : List[str]
        List of sample identifiers to be analyzed.

    Returns
    -------
    sample_update : pd.DataFrame
        Updated sample DataFrame with additional columns:
        - 'Mean_w' : Weighted mean of the sample ages.
        - 'StDev_w' : Weighted standard deviation of the sample ages.
        - 'SEsd_w' : Standard error of the weighted standard deviation.
        - 'SEiu_w' : Standard error based on individual uncertainties.
        - 'MSWD_w' : Mean Square Weighted Deviation, representing goodness of fit.
    
    Notes
    -----
    For samples identified as 'AFT' or 'ZFT', the function uses the mean and standard deviation directly from the 
    'Corrected_Date_Ma' and 'Corrected_Uncertainty_1σ_Ma' columns without further calculation. 
    """
    sample_update = samples.copy()

    for sample in sample_list:
        if sample[1].upper() in {'AFT', 'ZFT'}:
            mean = aliquots.loc[sample[0], 'Corrected_Date_Ma']
            stdev = aliquots.loc[sample[0], 'Corrected_Uncertainty_1σ_Ma']

            sample_update.loc[[sample[0]],'Mean_w'] = mean
            sample_update.loc[[sample[0]],'StDev_w'] = stdev
        else:
            # Filter aliquots for the sample and keep only non-outliers
            sample_df = aliquots[(aliquots.Sample == sample[0]) & (aliquots.outlier == 'keep')]
            ages = [x for x in sample_df.Corrected_Date_Ma]
            errors = [x for x in sample_df.Corrected_Uncertainty_1σ_Ma]

            n = len(sample_df)
            
            if n > 0:
                # Weighted mean (x̄_w)
                mean = (sum( x/y**2 for x,y in zip(ages,errors)) / sum(1/x**2 for x in errors))

                # Weighted standard deviation (SD_w)
                stdev = errors[0] if n <= 1 else np.sqrt(
                    sum((x - mean)**2 / y**2 for x, y in zip(ages, errors)) / ((n - 1) * sum(1 / y**2 for y in errors))
                )

                # Standard error based on SD (SEsd_w)
                se_sd = stdev / np.sqrt(n)

                # Standard error based on individual uncertainties (SEiu_w)
                se_iu = np.sqrt(1 / sum(1 / y**2 for y in errors))

                # Mean Square Weighted Deviation (MSWD_w)
                mswd = None if n <= 1 else sum((x - mean)**2 / y**2 for x, y in zip(ages, errors)) / (n - 1)

                # Add values to sample database
                sample_update.loc[[sample[0]],'Mean_w'] = mean
                sample_update.loc[[sample[0]],'StDev_w'] = stdev
                sample_update.loc[[sample[0]],'SEsd_w'] = se_sd
                sample_update.loc[[sample[0]],'SEiu_w'] = se_iu
                sample_update.loc[[sample[0]],'MSWD_w'] = mswd

    return sample_update

def calculateWeightedRelativeDeviationSummaryStats(samples: pd.DataFrame, aliquots: pd.DataFrame, sample_list: List[str])-> pd.DataFrame:
    """
    Calculates the summary statistics weighted based on relative deviation for each sample, including weighted mean, 
    standard deviation, errors, and MSWD.

    Parameters
    ----------
    samples : pd.DataFrame
        DataFrame containing sample information indexed by sample name.
    aliquots : pd.DataFrame
        DataFrame containing aliquot data indexed by aliquot name, where multiple aliquots may correspond to a single sample.
    sample_list : List[str]
        List of sample identifiers to be analyzed.

    Returns
    -------
    sample_update : pd.DataFrame
        Updated sample DataFrame with additional columns:
        - 'Mean_rw' : Weighted relative deviation mean of the sample ages.
        - 'StDev_rw' : Weighted standard deviation based on relative deviation.
        - 'SEsd_rw' : Standard error of the weighted standard deviation.
        - 'SEiu_rw' : Standard error based on individual uncertainties.
        - 'MSWD_rw' : Mean Square Weighted Deviation, representing goodness of fit.

    Notes
    -----
    For samples identified as 'AFT', the function directly assigns the mean and standard deviation values from the 
    'Corrected_Date_Ma' and 'Corrected_Uncertainty_1σ_Ma' columns without further calculation. 
    """
    sample_update = samples.copy()

    for sample in sample_list:
        if sample[1].upper() == 'AFT':
            mean = aliquots.loc[sample[0], 'Corrected_Date_Ma']
            stdev = aliquots.loc[sample[0], 'Corrected_Uncertainty_1σ_Ma']

            sample_update.loc[[sample[0]],'Mean_rw'] = mean
            sample_update.loc[[sample[0]],'StDev_rw'] = stdev
        else:
            # Filter aliquots for the sample and keep only non-outliers
            sample_df = aliquots[(aliquots.Sample == sample[0]) & (aliquots.outlier == 'keep')]
            ages = [x for x in sample_df.Corrected_Date_Ma]
            errors = [x for x in sample_df.Corrected_Uncertainty_1σ_Ma]
            
            # Calculate relative deviation (r_i = error_i / age_i)
            r = [err / age for age, err in zip(ages, errors)]

            n = len(sample_df)
            
            if n > 0:
                # Weighted relative mean
                mean = sum(x / r_val**2 for x, r_val in zip(ages, r)) / sum(1 / r_val**2 for r_val in r)

                # Weighted relative standard deviation
                stdev = errors[0] if n <= 1 else np.sqrt(
                    sum((x - mean)**2 / r_val**2 for x, r_val in zip(ages, r)) / ((n - 1) * sum(1 / r_val**2 for r_val in r))
                )

                # Standard error based on SD (SEsd_rw)
                se_sd = stdev / np.sqrt(n)

                # Standard error based on individual uncertainties (SEiu_rw)
                se_iu = mean * np.sqrt(1 / sum(1 / r_val**2 for r_val in r))

                # Mean Square Weighted Deviation (MSWD_rw)
                mswd = None if n <= 1 else sum((x - mean)**2 / err**2 for x, err in zip(ages, errors)) / (n - 1)

                # Add values to sample database
                sample_update.loc[[sample[0]],'Mean_rw'] = mean
                sample_update.loc[[sample[0]],'StDev_rw'] = stdev
                sample_update.loc[[sample[0]],'SEsd_rw'] = se_sd
                sample_update.loc[[sample[0]],'SEiu_rw'] = se_iu
                sample_update.loc[[sample[0]],'MSWD_rw'] = mswd

    return sample_update


#######################################################################
# Functions for plotting thermochronologic data 
#######################################################################
def plot_samples_eU_Rft(sample_list: List[str], aliquots: pd.DataFrame, radius: str, plot_histogram: bool = False, 
    bin_width: int = 10, kde_overlay: bool = False, savefig: bool = False) -> None:
    """
    Generates plots for each sample in `sample_list`, showing eU (effective Uranium) vs Corrected Age, 
    and Rft (or Rs) vs Corrected Age. Optionally includes a histogram of Corrected Age with configurable bin width 
    and optional KDE overlay.

    This function is based on Figure 4 from Flowers et al., 2022.

    Parameters
    ----------
    sample_list : List[str]
        List of sample identifiers to be plotted, which allows for grouping by transect, basin, or other categories.
    aliquots : pd.DataFrame
        DataFrame containing aliquot data, indexed by aliquot name, with multiple aliquots potentially corresponding 
        to a single sample.
    radius : str
        Specifies the radius column to use for plotting, either 'Rft' or 'Rs'.
    plot_histogram : bool, optional (default=False)
        If True, includes a histogram of Corrected Age for each sample.
    bin_width : int, optional (default=10)
        Sets the width of the bins for the histogram, applied if `plot_histogram` is True.
    kde_overlay : bool, optional (default=False)
        If True, overlays a Kernel Density Estimate (KDE) on the histogram to show the distribution more smoothly.
    savefig : bool, optional (default=False)
        If True, saves each plot as an image file in a new directory folder created for the output.

    Returns
    -------
    None
        This function outputs plots for each sample in `sample_list` and, if `savefig` is True, saves the plots in a
        designated folder. The plots illustrate the relationships between Corrected Age and either eU or radius values,
        with optional histogram and KDE overlay to further detail age distributions.
    
    Notes
    -----
    - `radius` must be specified as either 'Rft' or 'Rs' to indicate the appropriate radius column.
    - Plots generated follow the layout of Figure 4 from Flowers et al., 2022, enabling easy comparison to published data.
    """
    for sample in sample_list:
        grain_df = aliquots[aliquots.Sample == sample]
        (mineral,) = (set(grain_df.Mineral.to_list()))

        if mineral.upper() == 'ZHE':
            marker = 'D'
        elif mineral.upper() == 'AHE':
            marker = 'h'
        else:
            marker = 'o'

        ### Figure Set Up -----------------------------
        if plot_histogram:
            fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(22, 6), gridspec_kw={'width_ratios': [1, 1, 0.3]})
            plt.subplots_adjust(wspace=0.1)  # Adjust horizontal space between subplots
        else:
            fig, (ax1,ax2) = plt.subplots(1,2 , figsize = (17,4))

        r_label = f'{radius} (µm)'

        ### Data Plotting -----------------------------
        ## eU v Corrected Age
        eu = ax1.scatter(grain_df.eU, grain_df.Corrected_Date_Ma, c=grain_df[radius], cmap='OrRd',
                         marker=marker, edgecolors='black', linewidth=1, s=150)

        cbar = plt.colorbar(eu, ax=ax1, pad=0.05)
        cbar.ax.tick_params(labelsize=12, colors=GREY40)
        cbar.set_label(label=r_label, size=12, color=GREY40)

        # X-axis formatting for ax1
        if mineral.upper() == 'ZHE':
            x1_max = max(grain_df.eU) + 50
        else:
            x1_max = max(grain_df.eU) + 20
        x1_interval = 10 if x1_max < 120 else (20 if x1_max < 200 else (50 if x1_max < 600 else 100))
        ax1.set_xticks([x for x in np.arange(0, x1_max, x1_interval)])
        ax1.set_xlabel('eU (ppm)', fontname="Arial", fontsize=12, weight=500, color=GREY40)
        ax1.set_ylabel('Corrected Date (Ma)', fontname="Arial", fontsize=12, weight=500, color=GREY40)
        ax1.tick_params(axis="x", length=5, color=GREY91)
        ax1.tick_params(axis="y", length=8, color=GREY91)

        # Y-axis formatting for ax1 and ax2
        if mineral.upper() == 'ZHE':
            y_max = max(grain_df.Corrected_Date_Ma) + 75
        else:
            y_max = max(grain_df.Corrected_Date_Ma) + 20
        y_interval = 10 if y_max <= 100 else (20 if y_max <= 200 else 50)
        y_ticks = [x for x in np.arange(0, y_max, y_interval)]
        ax1.set_yticks(y_ticks)
        
        ## Rft v Corrected Age
        rft = ax2.scatter(grain_df[radius], grain_df.Corrected_Date_Ma, c=grain_df.eU, cmap='viridis',
                          marker=marker, edgecolors='black', linewidth=1, s=150)

        cbar2 = plt.colorbar(rft, ax=ax2, pad=0.05)
        cbar2.ax.tick_params(labelsize=12, colors=GREY40)
        cbar2.set_label(label='eU (ppm)', size=12, color=GREY40)

        x2_max = max(grain_df[radius]) + 20
        x2_interval = 10 if x2_max < 100 else 20
        ax2.set_xticks([x for x in np.arange(0, x2_max, x2_interval)])
        ax2.set_xlabel(r_label, fontname="Arial", fontsize=12, weight=500, color=GREY40)
        ax2.set_yticks([x for x in np.arange(0, y_max, y_interval)])
        ax2.set_ylabel('Corrected Date (Ma)', fontname="Arial", fontsize=12, weight=500, color=GREY40)
        ax2.tick_params(axis="x", length=5, color=GREY91)
        ax2.tick_params(axis="y", length=8, color=GREY91)

        if plot_histogram:
        ## Histogram of Corrected Date (y-axis) with optional KDE
            range_y = np.arange(0, y_max, bin_width)
            ax3.hist(grain_df.Corrected_Date_Ma, bins=range_y, orientation='horizontal', color='lightblue', edgecolor='black')
            
            # KDE overlay
            if kde_overlay:
                sns.kdeplot(y=grain_df.Corrected_Date_Ma, ax=ax3, color='darkblue', fill=True, bw_adjust=0.5)

            # Move y-axis label to the right and match y-ticks with ax1/ax2
            ax3.yaxis.set_label_position("right")
            ax3.yaxis.tick_right()
            ax3.set_ylim(ax1.get_ylim())  # Set y-limits to match ax1
            ax3.set_yticks(y_ticks)       # Set y-ticks to match ax1 and ax2
            ax3.set_yticklabels([str(int(tick)) for tick in y_ticks])  # Match y-tick labels

            ax3.set_ylabel('Corrected Date (Ma)', fontname="Arial", fontsize=12, weight=500, color=GREY40)
            ax3.set_xlabel('Frequency', fontname="Arial", fontsize=12, weight=500, color=GREY40)
            ax3.tick_params(axis="x", length=5, color=GREY91)
            ax3.tick_params(axis="y", length=8, color=GREY91)
        
        ### Final Formatting -----------------------------
        title = 'Sample {} - {}'.format(sample, mineral)
        fig.text(0.45, 0.93, title, color='black', fontsize=15, fontname="Arial", weight="bold")

        if savefig:
            pathlib.Path('Plots/eu_Rft').mkdir(parents=True, exist_ok=True)
            if plot_histogram:
                filepath = f'Plots/eu_Rft/sample_{sample}_eu_rft_plot_histogram.pdf'
            else:
                filepath = f'Plots/eu_Rft/sample_{sample}_eu_rft_plot.pdf'
            plt.savefig(filepath, dpi='figure', bbox_inches='tight', pad_inches=0.5)

        plt.show();

def confidence_intervals(transectData: pd.DataFrame, sample_data: pd.DataFrame, variable: str, chronometer: str, 
    weightedBy: str)-> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculates the regression line and confidence intervals for a given dataset.

    Parameters
    ----------
    transectData : pd.DataFrame
        DataFrame containing aliquot ages with columns 'Corrected_Date_Ma', 'Mean', 'StDev', and specified `variable`.
    sample_data : pd.DataFrame
        DataFrame containing sample data with columns 'Mineral', 'Mean', 'StDev', and specified `variable`.
    variable : str
        The independent variable for the regression, either 'Elevation_m' or 'Latitude'.
    chronometer : str
        The chronometer type to filter sample data, typically 'AHe', 'ZHe', 'AFT', or 'ZFT'.
    weightedBy : str
        Specifies the weighting method, options are 'unweighted', 'inverse_variance', or 'relative_deviation'.

    Returns
    -------
    p_x : np.ndarray
        Array of new test x-values being predicted for.
    p_y : np.ndarray
        Array of predicted y-values for the linear regression line.
    lower : np.ndarray
        Lower confidence interval boundary for the predicted y-values.
    upper : np.ndarray
        Upper confidence interval boundary for the predicted y-values.
    """
    # Define column mapping based on weightedBy
    column_mapping = {
        'unweighted': {'Mean': 'Mean', 'StDev': 'StDev', 'SEsd': 'SEsd', 'SEiu': 'SEiu'},
        'inverse_variance': {'Mean': 'Mean_w', 'StDev': 'StDev_w', 'SEsd': 'SEsd_w', 'SEiu': 'SEiu_w'},
        'relative_deviation': {'Mean': 'Mean_rw', 'StDev': 'StDev_rw', 'SEsd': 'SEsd_rw', 'SEiu': 'SEiu_rw'}
    }

    # Select columns based on weightedBy
    selected_columns = column_mapping.get(weightedBy, column_mapping['unweighted'])

    # Filter samples by the specified chronometer type
    transectSamples = sample_data[sample_data.Mineral == chronometer]

    # Convert relevant data columns to numpy arrays
    yi = transectData['Corrected_Date_Ma'].to_numpy()
    y = transectData[selected_columns['Mean']].to_numpy()
    x = transectData[variable].to_numpy()
    w = transectData[selected_columns['StDev']].to_numpy()
    y2 = transectSamples[selected_columns['Mean']].to_numpy()
    x2 = transectSamples[variable].to_numpy()
    w2 = transectSamples[selected_columns['StDev']].to_numpy()

    # Calculate weighted regression parameters based on chronometer type
    S, Sx, Sxx, Sy, Sxy = sum(1 / w2**2), sum(x2 / w2**2), sum((x2**2) / w2**2), sum(y2 / w2**2), sum((x2 * y2) / w2**2)
    if chronometer in ('AFT', 'ZFT'):
        Sx = sum(x / w2**2)
        Sxy = sum(x * y2 / w2**2)
    delta = S * Sxx - Sx**2

    # Calculate regression coefficients (slope, intercept) and their errors
    m = ((S * Sxy) - (Sx * Sy)) / delta
    m_err = S / delta
    b = ((Sxx * Sy) - (Sx * Sxy)) / delta
    b_err = Sxx / delta

    # Fit regression line and calculate residuals
    fit_y = m * x + b
    residuals = yi - fit_y

    # Define the variables needed to calculate the confidence interval
    mean_x = np.mean(x) # mean of x
    n = len(yi) # number of samples in original fit
    tstat = t.ppf(0.975, n-1) # appropriate t value (0.975 for 95% confidence interval)
    s_err = np.sum(residuals**2) # sum of the squares of the residuals

    # create series of new test x-values to predict for
    p_x = np.linspace(np.min(x),np.max(x),26)

    confs = tstat * np.sqrt((s_err / (n - 2)) * (1.0 / n + ((p_x - mean_x)**2) / ((np.sum(x**2)) - n * (mean_x**2))))

    # Predict y-values and calculate confidence intervals
    p_y = m * p_x + b
    lower, upper = p_y - confs, p_y + confs

    return p_x, p_y, lower, upper

def filter_regression_data(df: pd.DataFrame, mineral: str, exclude_outliers: bool, exclude_samples: Optional[List[str]], exclude_aliquots: Optional[List[str]]
    ) -> pd.DataFrame:
    """
    Filter data based on specific exclusions for outliers, samples, and aliquots.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing mineral data.
    mineral : str
        The mineral type to filter (e.g., 'AHe', 'ZHe').
    exclude_outliers : bool
        If True, exclude data points labeled as outliers.
    exclude_samples : Optional[List[str]]
        List of sample identifiers to exclude from the data.
    exclude_aliquots : Optional[List[str]]
        List of aliquot identifiers to exclude from the data.

    Returns
    -------
    filtered_df : pd.DataFrame
        Filtered DataFrame with specified outliers, samples, and aliquots removed.
    """
    filtered_df = df[df['Mineral'] == mineral]
    if exclude_outliers:
        filtered_df = filtered_df[filtered_df['outlier'] != 'reject']
    if exclude_samples:
        filtered_df = filtered_df[~filtered_df['Sample'].isin(exclude_samples)]
    if exclude_aliquots:
        filtered_df = filtered_df[~filtered_df['Aliquot'].isin(exclude_aliquots)]
        for val in exclude_aliquots:
            df.at[val, 'outlier'] = 'reject'
    return filtered_df

def lighten_color(color: str, amount: float=0.7)-> List[float]:
    """
    Lightens a given color by interpolating it with white.

    Parameters
    ----------
    color : str
        Color to lighten, specified as a valid matplotlib color string (e.g., 'blue', '#333333').
    amount : float
        Amount to lighten the color, where a higher value results in a lighter color. Defaults to 0.7.

    Returns
    -------
    List[float]
        RGBA values of the lightened color.
    """
    # Convert the color to RGB
    base_color = mcolors.to_rgba(color)

    # Interpolate between the color and white
    new_color = [1 - (1 - component) * (1 - amount) for component in base_color[:3]] + [base_color[3]]
    return new_color

def plot_regression(ax: Axes, df: pd.DataFrame, plot_df: pd.DataFrame, y_variable: str, color: str, 
    chronometer: str, weightedBy: str = 'unweighted') -> None:
    """
    Plots the regression line and confidence intervals on the given axes.

    Parameters
    ----------
    ax : Axes
        Matplotlib Axes object on which to plot.
    df : pd.DataFrame
        DataFrame containing data points to calculate the regression.
    plot_df : pd.DataFrame
        DataFrame with data for plotting the regression line and confidence intervals.
    y_variable : str
        The variable to plot on the y-axis (e.g., 'Latitude', 'Longitude', 'Elevation', 'Structural Level').
    color : str
        Color for the confidence interval fill.
    chronometer : str
        Chronometer type (e.g., 'AHe', 'ZHe') used in the regression calculation.
    weightedBy : str, default='unweighted'
        Specifies the weighting method, with options such as 'unweighted'.

    Returns
    -------
    None
    """
    if len(df['Sample'].unique()) > 1:
        p_x, p_y, lower, upper = confidence_intervals(df, plot_df, y_variable, chronometer=chronometer, weightedBy = weightedBy)
        ax.plot(p_y, p_x, color="grey", label="Regression Line")
        ax.fill_betweenx(p_x, lower, upper, color=color, alpha=0.3)

def plotElevationProfile(samples: pd.DataFrame, x_variable: str, transect: Optional[str] = None, colorBy: Optional[str] = None,
    x_bounds: Optional[Tuple[float, float]] = None, y_bounds: Optional[Tuple[float, float]] = None,
    label_samples: bool = True, label_offset: Tuple[float, float] = (0, 0),
    AHeColor: str = 'darkmagenta', AHeMarker: str = 'h', AHeMarkerSize: int = 12,
    ZHeColor: str = 'cadetblue', ZHeMarker: str = 'D', ZHeMarkerSize: int = 10,
    AFTColor: str = 'gray', AFTMarker: str = '^', AFTMarkerSize: int = 10,
    ZFTColor: str = 'forestgreen', ZFTMarker: str = 'o', ZFTMarkerSize: int = 10, 
    savefig: bool = False, savefigFileName: Optional[str] = None, saveFolder: str = 'Plots') -> None:
    """
    Plots the Elevation profile against Latitude or Longitude (based on user input) for the sample data, allowing customization by 
    chronometer, transect, and appearance settings.

    Parameters
    ----------
    samples : pd.DataFrame
        Dataframe containing sample data with columns for Longitude, Elevation_m, Transect, Mineral, Sample, and other relevant details.
    x_variable : str
        Specific X-variable to plot.
    transect : Optional[str]
        Specific transect to plot; defaults to None to plot all data.
    colorBy : Optional[str]
        Method for color-coding data points, 'chronometer' or 'transect'. Default is None for uniform color.
    x_bounds : Optional[Tuple[float, float]]
        Manually set x-axis bounds (lower_x, upper_x).
    y_bounds : Optional[Tuple[float, float]]
        Manually set y-axis bounds (lower_y, upper_y).
    label_samples : bool
        If True, labels each sample on the plot.
    label_offset : tuple, optional
        Offset for sample labels (x, y); default is (0, 0).
    AHeColor : str, optional
        Color for the AHe marker; default is 'darkmagenta'.
    AHeMarker : str, optional
        Style for the AHe marker; default is 'h'.
    AHeMarkerSize : int, optional
        Size for the AHe marker; default is 12.
    ZHeColor : str, optional
        Color for the ZHe marker; default is 'cadetblue'.
    ZHeMarker : str, optional
        Style for the ZHe marker; default is 'D'.
    ZHeMarkerSize : int, optional
        Size for the ZHe marker; default is 10.
    AFTColor : str, optional
        Color for the AFT marker; default is 'gray'.
    AFTMarker : str, optional
        Style for the AFT marker; default is '^'.
    AFTMarkerSize : int, optional
        Size for the AFT marker; default is 10.
    ZFTColor : str, optional
        Color for the ZFT marker; default is 'forestgreen'.
    ZFTMarker : str, optional
        Style for the ZFT marker; default is 'o'.
    ZFTMarkerSize : int, optional
        Size for the ZFT marker; default is 10.
    savefig : bool, optional
        If True, saves the plot.
    savefigFileName : str, optional
        Filename to save the figure under, if saving.
    saveFolder : str, optional
        Folder to save figures.

    Returns
    -------
    None
        Displays the plot; saves it if savefig is True.
    """
    # Set up transect colors for 'transect' color option
    transect_categories = samples['Transect'].unique()
    transect_colors = {category: '#' + ''.join(random.choice('0123456789ABCDEF') for _ in range(6)) for category in transect_categories}

    # Filter data based on transect if specified
    plot_data = samples[samples['Transect'] == transect] if transect else samples.copy()
    
    # Set up plot
    fig, ax = plt.subplots(figsize=(15, 8))

    for row in plot_data.itertuples():
        if colorBy == 'chronometer':
            if row.Mineral == 'AHe':
                marker, color, size = AHeMarker, AHeColor, AHeMarkerSize
            elif row.Mineral == 'ZHe':
                marker, color, size = ZHeMarker, ZHeColor, ZHeMarkerSize
            elif row.Mineral == 'AFT':
                marker, color, size = AFTMarker, AFTColor, AFTMarkerSize
            elif row.Mineral == 'ZFT':
                marker, color, size = ZFTMarker, ZFTColor, ZFTMarkerSize
        elif colorBy == 'transect':
            marker, color, size = 'o', transect_colors[row.Transect], 12
        else:  # Default color when colorBy is None
            marker, color, size = 'o', 'black', 10

        # Plot sample point
        ax.errorbar(x=getattr(row, x_variable), y=row.Elevation_m, ecolor=color, mfc=color, marker=marker,
                    ms=size, mec='black', mew=1)

        # Optional sample label, only if the sample is within x_bounds
        if label_samples and (
            (x_bounds is None or (x_bounds[0] <= row[x_variable] <= x_bounds[1])) and
            (y_bounds is None or (y_bounds[0] <= row.Elevation_m <= y_bounds[1]))
        ):
            ax.text(
                getattr(row, x_variable) + label_offset[0], 
                row.Elevation_m + label_offset[1], 
                str(row.Sample),
                fontsize=12, 
                fontname="Arial", 
                rotation=90, 
                va='center', 
                ha='left', 
                color='gray'
            )

    # Chronometer or transect legend
    if colorBy == 'chronometer':
        legend_handles = [
            mlines.Line2D([], [], color=AHeColor, marker=AHeMarker, linestyle='', markersize=AHeMarkerSize, label='AHe'),
            mlines.Line2D([], [], color=ZHeColor, marker=ZHeMarker, linestyle='', markersize=ZHeMarkerSize, label='ZHe'),
            mlines.Line2D([], [], color=AFTColor, marker=AFTMarker, linestyle='', markersize=AFTMarkerSize, label='AFT'),
            mlines.Line2D([], [], color=ZFTColor, marker=ZFTMarker, linestyle='', markersize=ZFTMarkerSize, label='ZFT')
        ]
        ax.legend(handles=legend_handles, title="Chronometers")
    elif colorBy == 'transect':
        legend_handles = [plt.Line2D([], [], color=transect_colors[cat], marker='o', linestyle='', label=cat)
                          for cat in transect_categories]
        ax.legend(handles=legend_handles, title="Transects")

    # Set plot bounds if specified
    if x_bounds:
        ax.set_xlim(x_bounds)
    if y_bounds:
        ax.set_ylim(y_bounds)

    # Axis Specifications
    ax.spines["left"].set_color('k')
    ax.spines["bottom"].set_color('k')
    ax.spines["right"].set_color('k')
    ax.spines["top"].set_color('k')

    ax.grid(True)

    ax.set_xlabel(x_variable,fontname= "Arial",fontsize=12,weight=500,color=GREY40)
    ax.set_ylabel('Elevation (m)',fontname= "Arial",fontsize=12,weight=500,color=GREY40)

    # Customize tick label font properties
    font_properties = fm.FontProperties(
        family="Arial",  # Font name
        size=10,         # Font size
        weight="bold"    # Font weight
    )

    # Set tick parameters for color and length
    ax.tick_params(axis='both', color=GREY91, length=5)

    # Apply font properties to x-axis and y-axis tick labels
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontproperties(font_properties)
        label.set_color("darkgray")  # Set font color

    # Save figure if specified
    if savefig:
        pathlib.Path(saveFolder).mkdir(parents=True, exist_ok=True)
        filepath = f"{saveFolder}/{savefigFileName or 'elevation_plot'}.pdf"
        plt.savefig(filepath, bbox_inches='tight', pad_inches=0.5)

    plt.show()

def plotAgeVersus(samples: pd.DataFrame, aliquots: pd.DataFrame, x_variable: str, y_variable: str,
    transect: Optional[str] = None,show_aliquots: bool = True, 
    weightedBy : str = 'unweighted', plot_SE: bool = False, SE_basedOn: str = 'max',
    x_bounds: Optional[Tuple[float, float]] = None, y_bounds: Optional[Tuple[float, float]] = None,
    label_samples: bool = True, label_offset: Tuple[float, float] = (0, 0), 
    AHeColor: str = 'darkmagenta', AHeMarker: str = 'h', AHeMarkerSize: int = 12,
    ZHeColor: str = 'cadetblue', ZHeMarker: str = 'D', ZHeMarkerSize: int = 10,
    AFTColor: str = 'gray', AFTMarker: str = '^', AFTMarkerSize: int = 10,
    ZFTColor: str = 'forestgreen', ZFTMarker: str = 'o', ZFTMarkerSize: int = 10,
    savefig: bool = False, savefigFileName: Optional[str] = None, saveFolder: str = 'Plots') -> None:
    """
    Plots the Cooling vs. a user-specified variable for the data, allowing customization of the x- and y-variables, transect, 
    and appearance settings.

    Allows customization of markers, colors, labels, and other details for different thermochronometers.

    Parameters
    ----------
    samples : pd.DataFrame
        Dataframe containing sample data with columns for Longitude, Elevation_m, Transect, Mineral, Sample, and other relevant details.
    aliquots : pd.DataFrame
        DataFrame containing aliquot data indexed by aliquot name, where multiple aliquots may correspond to a single sample.
    x_variable : str
        Specific X-variable to plot.
    y_variable : str
        Specific Y-variable to plot.
    transect : Optional[str]
        Specific transect to plot; defaults to None to plot all data.
    show_aliquots : bool
        If True, shows the individual aliquot ages for each sample.
    weightedBy : str, default='unweighted'
        Specifies the weighting method, with options such as 'unweighted'.
    plot_SE : bool, default=False
        Whether to plot standard error bars.
    SE_basedOn : str, default='max'
        Basis for SE calculation; 'max' uses the maximum standard error.
    x_bounds : Optional[Tuple[float, float]]
        Manually set x-axis bounds (lower_x, upper_x).
    y_bounds : Optional[Tuple[float, float]]
        Manually set y-axis bounds (lower_y, upper_y).
    label_samples : bool
        If True, labels each sample on the plot.
    label_offset : tuple, optional
        Offset for sample labels (x, y); default is (0, 0).
    AHeColor : str, optional
        Color for the AHe marker; default is 'darkmagenta'.
    AHeMarker : str, optional
        Style for the AHe marker; default is 'h'.
    AHeMarkerSize : int, optional
        Size for the AHe marker; default is 12.
    ZHeColor : str, optional
        Color for the ZHe marker; default is 'cadetblue'.
    ZHeMarker : str, optional
        Style for the ZHe marker; default is 'D'.
    ZHeMarkerSize : int, optional
        Size for the ZHe marker; default is 10.
    AFTColor : str, optional
        Color for the AFT marker; default is 'gray'.
    AFTMarker : str, optional
        Style for the AFT marker; default is '^'.
    AFTMarkerSize : int, optional
        Size for the AFT marker; default is 10.
    ZFTColor : str, optional
        Color for the ZFT marker; default is 'forestgreen'.
    ZFTMarker : str, optional
        Style for the ZFT marker; default is 'o'.
    ZFTMarkerSize : int, optional
        Size for the ZFT marker; default is 10.
    savefig : bool, optional
        If True, saves the plot.
    savefigFileName : str, optional
        Filename to save the figure under, if saving.
    saveFolder : str, optional
        Folder to save figures.

    Returns
    -------
    None
        Displays the plot; saves it if savefig is True.
    """
    # Set up data for plotting 
    full_aliquot_df = pd.merge(aliquots,samples[['Mean','StDev','Latitude','Longitude','Elevation_m','Structural_Level', 'Transect']],
                     on = 'Sample', how = 'left')

    full_aliquot_df.set_index('Aliquot', inplace = True, drop = False)

    plot_data = samples[samples['Transect'] == transect] if transect else samples.copy()
    full_aliquot_df = full_aliquot_df[full_aliquot_df['Transect'] == transect] if transect else full_aliquot_df

    # Figure Set up 
    fig, ax = plt.subplots(figsize=(15, 8))
    
    outlier_markers = {'keep': 'o', 'reject': 'X'}
    colors = {'keep': 'black', 'reject': 'gray'}
    mineral_markers = {'AHe': AHeMarker, 'ZHe': ZHeMarker, 'AFT': AFTMarker, 'ZFT': ZFTMarker}

    # Markers for sample means with error bars
    mineral_styles = {
        'AHe': (AHeMarker, AHeColor, AHeMarkerSize),
        'ZHe': (ZHeMarker, ZHeColor, ZHeMarkerSize),
        'AFT': (AFTMarker, AFTColor, AFTMarkerSize),
        'ZFT': (ZFTMarker, ZFTColor, ZFTMarkerSize),
    }

    # Define the column mapping based on the weightedBy option
    column_mapping = {
        'unweighted': {'Mean': 'Mean', 'StDev': 'StDev', 'SEsd': 'SEsd', 'SEiu': 'SEiu'},
        'inverse_variance': {'Mean': 'Mean_w', 'StDev': 'StDev_w', 'SEsd': 'SEsd_w', 'SEiu': 'SEiu_w'},
        'relative_deviation': {'Mean': 'Mean_rw', 'StDev': 'StDev_rw', 'SEsd': 'SEsd_rw', 'SEiu': 'SEiu_rw'}
    }

    # Get the correct columns based on weightedBy
    selected_columns = column_mapping.get(weightedBy, column_mapping['unweighted']) 

    # # Function to lighten a color (for standard error option)
    # def lighten_color(color, amount=0.7):
    #     # Convert the color to RGB
    #     base_color = mcolors.to_rgba(color)
    #     # Interpolate between the color and white
    #     new_color = [1 - (1 - component) * (1 - amount) for component in base_color[:3]] + [base_color[3]]
    #     return new_color

    if y_variable == 'Age':
        # Plotting All Data Points 
        if show_aliquots:
            # Create separate datasets for rejected and non-rejected data points
            rejected_df = full_aliquot_df[full_aliquot_df['outlier'] == 'reject']
            non_rejected_df = full_aliquot_df[full_aliquot_df['outlier'] == 'keep']

            # Plot rejected data points with X markers
            ax.scatter(
                x=rejected_df[x_variable],
                y=rejected_df.Corrected_Date_Ma,
                marker=outlier_markers['reject'],
                color=colors['reject']
            )

            # Plot non-rejected data points with mineral markers
            for mineral, marker in mineral_markers.items():
                df = non_rejected_df[non_rejected_df['Mineral'] == mineral]
                ax.scatter(
                    x=df[x_variable],
                    y=df.Corrected_Date_Ma,
                    marker=marker,
                    color=colors['keep']
                )
 
        # Plot Sample Means w/ Errorbars 
        for index, row in plot_data.iterrows():
            marker, color, size = mineral_styles.get(row['Mineral'], ('o', 'maroon', 10))
            label = row.Sample

            ax.errorbar(x=plot_data.loc[index][x_variable], y=row[selected_columns['Mean']],
                        yerr=row[selected_columns['StDev']], ecolor = color, elinewidth=3,
                        mfc = color, marker = marker, ms = size, mec = 'black', mew = 1
            )
            if plot_SE:
                # Apply the lighter color for plotting SE
                lighter_color = lighten_color(color, amount=0.7)  # adjust 'amount' for more brightness

                # Select the appropriate xerr column based on SE_basedOn
                if SE_basedOn == 'SD':
                    err_value = row[selected_columns['SEsd']]
                elif SE_basedOn == 'IU':
                    err_value = row[selected_columns['SEiu']]
                elif SE_basedOn == 'max':
                    err_value = max(row[selected_columns['SEsd']], row[selected_columns['SEiu']])
                else:
                    raise ValueError("Invalid option for SE_basedOn. Choose from 'SD', 'IU', or 'max'.")

                ax.errorbar(
                    plot_data.loc[index][x_variable], row[selected_columns['Mean']],
                    yerr=err_value, ecolor=lighter_color, fmt='none', mew=1, elinewidth=3
                )

            # Sample Labels
            if label_samples:
                label = row.Sample
                
                if x_bounds and y_bounds:
                    if row[selected_columns['Mean']] < y_bounds[1] and row[selected_columns['Mean']] > y_bounds[0] and plot_data.loc[index][x_variable] < x_bounds[1] and plot_data.loc[index][x_variable] > x_bounds[0]:
                                ax.text(x = plot_data.loc[index][x_variable] + label_offset[0],
                                y = row[selected_columns['Mean']] + label_offset[1],
                                s = label,
                                fontsize = 12,
                                fontname = "Arial",
                                rotation = 90,
                                va = 'center',
                                ha = 'left',
                                color = GREY60,
                                weight = 'book',
                                style = 'italic')
                
                elif y_bounds:
                    if row[selected_columns['Mean']] < y_bounds[1] and row[selected_columns['Mean']] > y_bounds[0]:
                        ax.text(x = plot_data.loc[index][x_variable] + label_offset[0],
                                y = row[selected_columns['Mean']] + label_offset[1],
                                s = label,
                                fontsize = 12,
                                fontname = "Arial",
                                rotation = 90,
                                va = 'center',
                                ha = 'left',
                                color = GREY60,
                                weight = 'book',
                                style = 'italic')
                elif x_bounds:
                    if plot_data.loc[index][x_variable] < x_bounds[1] and plot_data.loc[index][x_variable] > x_bounds[0]:
                        ax.text(x = plot_data.loc[index][x_variable] + label_offset[0],
                                y = row[selected_columns['Mean']] + label_offset[1],
                                s = label,
                                fontsize = 12,
                                fontname = "Arial",
                                rotation = 90,
                                va = 'center',
                                ha = 'left',
                                color = GREY60,
                                weight = 'book',
                                style = 'italic')
                else:
                    ax.text(x = plot_data.loc[index][x_variable] + label_offset[0],
                                y = row[selected_columns['Mean']] + label_offset[1],
                                s = label,
                                fontsize = 12,
                                fontname = "Arial",
                                rotation = 90,
                                va = 'center',
                                ha = 'left',
                                color = GREY60,
                                weight = 'book',
                                style = 'italic')

        ##########################
        ### Axes and Spine Customization -----------------------------
        ax.spines["left"].set_color('k')
        ax.spines["bottom"].set_color('k')
        ax.spines["right"].set_color('k')
        ax.spines["top"].set_color('k')

        ax.grid(True, color = GREY91, linestyle = '-', linewidth = 1)

        if y_bounds:
            ax.set_ylim(y_bounds)
        else: 
            plt.gca().set_ylim(bottom=0)

        if x_bounds:
            ax.set_xlim(x_bounds)

        ax.tick_params(axis="x", length=5, color=GREY91)
        ax.tick_params(axis="y", length=5, color=GREY91)
        
        # Labels
        if x_variable == 'Latitude':
            x_axis_label = 'Latitude'
        elif x_variable == 'Longitude':
            x_axis_label = 'Longitude'
        elif x_variable == 'Structural_Level':
            x_axis_label = 'Structural Level'
        elif x_variable == 'Elevation_m':
            x_axis_label = 'Elevation (m)'
        else:
            x_axis_label = 'Y variable'

        ax.set_xlabel(x_axis_label,
                fontname= "Arial",
                fontsize=12,
                weight=500,
                color=GREY40
            )

        ax.set_ylabel('Corrected Cooling Age (Ma)',
                        fontname= "Arial",
                        fontsize=12,
                        weight=500,
                        color=GREY40
                    )

        
    elif x_variable == 'Age':
        # Plotting All Data Points -----------------------------
        if show_aliquots:

            # Create separate datasets for rejected and non-rejected data points
            rejected_df = full_aliquot_df[full_aliquot_df['outlier'] == 'reject']
            non_rejected_df = full_aliquot_df[full_aliquot_df['outlier'] == 'keep']

            # Plot rejected data points with X markers
            ax.scatter(
                x=rejected_df.Corrected_Date_Ma,
                y=rejected_df[y_variable],
                marker=outlier_markers['reject'],
                color=colors['reject']
            )

            # Plot non-rejected data points with mineral markers
            for mineral, marker in mineral_markers.items():
                df = non_rejected_df[non_rejected_df['Mineral'] == mineral]
                ax.scatter(
                    x=df.Corrected_Date_Ma,
                    y=df[y_variable],
                    marker=marker,
                    color=colors['keep']
                )
        
        ##########################
        ### Sample Means w/ Errorbars -----------------------------
        for index, row in plot_data.iterrows():
            marker, color, size = mineral_styles.get(row['Mineral'], ('o', 'maroon', 10))
            label = row.Sample

            ax.errorbar(x=row[selected_columns['Mean']], y=plot_data.loc[index][y_variable],
                        xerr=row[selected_columns['StDev']], ecolor = color, elinewidth=3,
                        mfc = color, marker = marker, ms = size, mec = 'black', mew = 1
            )

            if plot_SE:
                # Function to lighten a color
                def lighten_color(color, amount=0.7):
                    # Convert the color to RGB
                    base_color = mcolors.to_rgba(color)
                    # Interpolate between the color and white
                    new_color = [1 - (1 - component) * (1 - amount) for component in base_color[:3]] + [base_color[3]]
                    return new_color

                # Apply the lighter color for plotting SE
                lighter_color = lighten_color(color, amount=0.7)  # adjust 'amount' for more brightness

                # Select the appropriate xerr column based on SE_basedOn
                if SE_basedOn == 'SD':
                    err_value = row[selected_columns['SEsd']]
                elif SE_basedOn == 'IU':
                    err_value = row[selected_columns['SEiu']]
                elif SE_basedOn == 'max':
                    err_value = max(row[selected_columns['SEsd']], row[selected_columns['SEiu']])
                else:
                    raise ValueError("Invalid option for SE_basedOn. Choose from 'SD', 'IU', or 'max'.")

                ax.errorbar(
                    row[selected_columns['Mean']], plot_data.loc[index][y_variable],
                    xerr=err_value, ecolor=lighter_color, fmt='none', mew=1, elinewidth=3
                )

            # Sample Labels
            if label_samples:
                label = row.Sample
                
                if x_bounds and y_bounds:
                    if plot_data.loc[index][y_variable] < y_bounds[1] and plot_data.loc[index][y_variable] > y_bounds[0] and row[selected_columns['Mean']] < x_bounds[1] and row[selected_columns['Mean']] > x_bounds[0]:
                                ax.text(x = row[selected_columns['Mean']] + label_offset[0],
                                y = plot_data.loc[index][y_variable] + label_offset[1],
                                s = label,
                                fontsize = 12,
                                fontname = "Arial",
                                rotation = 90,
                                va = 'center',
                                ha = 'left',
                                color = GREY60,
                                weight = 'book',
                                style = 'italic')
                
                elif y_bounds:
                    if plot_data.loc[index][y_variable] < y_bounds[1] and plot_data.loc[index][y_variable] > y_bounds[0]:
                        ax.text(x = row[selected_columns['Mean']] + label_offset[0],
                                y = plot_data.loc[index][y_variable] + label_offset[1],
                                s = label,
                                fontsize = 12,
                                fontname = "Arial",
                                rotation = 90,
                                va = 'center',
                                ha = 'left',
                                color = GREY60,
                                weight = 'book',
                                style = 'italic')
                elif x_bounds:
                    if row[selected_columns['Mean']] < x_bounds[1] and row[selected_columns['Mean']] > x_bounds[0]:
                        ax.text(x = row[selected_columns['Mean']] + label_offset[0],
                                y = plot_data.loc[index][y_variable] + label_offset[1],
                                s = label,
                                fontsize = 12,
                                fontname = "Arial",
                                rotation = 90,
                                va = 'center',
                                ha = 'left',
                                color = GREY60,
                                weight = 'book',
                                style = 'italic')
                else:
                    ax.text(x = row[selected_columns['Mean']] + label_offset[0],
                                y = plot_data.loc[index][y_variable] + label_offset[1],
                                s = label,
                                fontsize = 12,
                                fontname = "Arial",
                                rotation = 90,
                                va = 'center',
                                ha = 'left',
                                color = GREY60,
                                weight = 'book',
                                style = 'italic')

        ##########################
        ### Axes and Spine Customization -----------------------------
        ax.spines["left"].set_color('k')
        ax.spines["bottom"].set_color('k')
        ax.spines["right"].set_color('k')
        ax.spines["top"].set_color('k')

        ax.grid(True, color = GREY91, linestyle = '-', linewidth = 1)

        if y_bounds:
            ax.set_ylim(y_bounds)

        if x_bounds:
            ax.set_xlim(x_bounds)
        else:
            plt.gca().set_xlim(left=0)

        ax.tick_params(axis="x", length=5, color=GREY91)
        ax.tick_params(axis="y", length=5, color=GREY91)
        
        # Labels
        if y_variable == 'Latitude':
            y_axis_label = 'Latitude'
        elif y_variable == 'Longitude':
            y_axis_label = 'Longitude'
        elif y_variable == 'Structural_Level':
            y_axis_label = 'Structural Level'
        elif y_variable == 'Elevation_m':
            y_axis_label = 'Elevation (m)'
        else:
            y_axis_label = 'Y variable'

        ax.set_ylabel(y_axis_label,
                fontname= "Arial",
                fontsize=12,
                weight=500,
                color=GREY40
            )

        ax.set_xlabel('Corrected Cooling Age (Ma)',
                        fontname= "Arial",
                        fontsize=12,
                        weight=500,
                        color=GREY40
                    )

    ##########################
    ### Show and Save Figure -----------------------------
    if savefig:
        pathlib.Path(saveFolder).mkdir(parents=True, exist_ok=True)
        if savefigFileName:
            filepath = '{}/{}.pdf'.format(saveFolder, savefigFileName)
        else:
            filepath = '{}/age_longitude_by_analysis_type.pdf'.format(saveFolder)
        plt.savefig(filepath, dpi = 'figure', bbox_inches='tight', pad_inches = 0.5)

    fig.tight_layout()
    plt.show();   

def plotAgeVersus_wHistogram(samples: pd.DataFrame, aliquots: pd.DataFrame, x_variable: str, y_variable: str,
    transect: Optional[str] = None,show_aliquots: bool = True, 
    weightedBy : str = 'unweighted', plot_SE: bool = False, SE_basedOn: str = 'max',
    x_bounds: Optional[Tuple[float, float]] = None,y_bounds: Optional[Tuple[float, float]] = None,
    label_samples: bool = True, label_offset: Tuple[float, float] = (0, 0), bin_width: int = 10,
    AHeColor: str = 'darkmagenta', AHeMarker: str = 'h', AHeMarkerSize: int = 12,
    ZHeColor: str = 'cadetblue', ZHeMarker: str = 'D', ZHeMarkerSize: int = 10,
    AFTColor: str = 'gray', AFTMarker: str = '^', AFTMarkerSize: int = 10,
    ZFTColor: str = 'forestgreen', ZFTMarker: str = 'o', ZFTMarkerSize: int = 10,
    savefig: bool = False, savefigFileName: Optional[str] = None, saveFolder: str = 'Plots') -> None:
    """
    Plots Cooling Age against a specified variable (e.g., Latitude, Longitude, Elevation, Structural Level),
    with a histogram of aliquot cooling ages for each mineral type.

    Allows customization of markers, colors, labels, and regression details for different thermochronometers.

    Parameters
    ----------
    samples : pd.DataFrame
        Dataframe containing sample data with columns for Longitude, Elevation_m, Transect, Mineral, Sample, and other relevant details.
    aliquots : pd.DataFrame
        DataFrame containing aliquot data indexed by aliquot name, where multiple aliquots may correspond to a single sample.
    x_variable : str
        Specific X-variable to plot.
    y_variable : str
        Specific Y-variable to plot.
    transect : Optional[str]
        Specific transect to plot; defaults to None to plot all data.
    show_aliquots : bool
        If True, shows the individual aliquot ages for each sample.
    weightedBy : str, default='unweighted'
        Specifies the weighting method, with options such as 'unweighted'.
    plot_SE : bool, default=False
        Whether to plot standard error bars.
    SE_basedOn : str, default='max'
        Basis for SE calculation; 'max' uses the maximum standard error.
    x_bounds : Optional[Tuple[float, float]]
        Manually set x-axis bounds (lower_x, upper_x).
    y_bounds : Optional[Tuple[float, float]]
        Manually set y-axis bounds (lower_y, upper_y).
    label_samples : bool
        If True, labels each sample on the plot.
    label_offset : tuple, optional
        Offset for sample labels (x, y); default is (0, 0).
    bin_width : int, optional
        Bin width for the historgram; default is 10.
    AHeColor : str, optional
        Color for the AHe marker; default is 'darkmagenta'.
    AHeMarker : str, optional
        Style for the AHe marker; default is 'h'.
    AHeMarkerSize : int, optional
        Size for the AHe marker; default is 12.
    ZHeColor : str, optional
        Color for the ZHe marker; default is 'cadetblue'.
    ZHeMarker : str, optional
        Style for the ZHe marker; default is 'D'.
    ZHeMarkerSize : int, optional
        Size for the ZHe marker; default is 10.
    AFTColor : str, optional
        Color for the AFT marker; default is 'gray'.
    AFTMarker : str, optional
        Style for the AFT marker; default is '^'.
    AFTMarkerSize : int, optional
        Size for the AFT marker; default is 10.
    ZFTColor : str, optional
        Color for the ZFT marker; default is 'forestgreen'.
    ZFTMarker : str, optional
        Style for the ZFT marker; default is 'o'.
    ZFTMarkerSize : int, optional
        Size for the ZFT marker; default is 10.
    savefig : bool, optional
        If True, saves the plot.
    savefigFileName : str, optional
        Filename to save the figure under, if saving.
    saveFolder : str, optional
        Folder to save figures.

    Returns
    -------
    None
        Displays the plot; saves it if savefig is True.
    """
    # Set up data for plotting -----------------------------
    full_aliquot_df = pd.merge(aliquots,samples[['Mean','StDev','Latitude','Longitude','Elevation_m','Transect',
                                                'Structural_Level']], on = 'Sample', how = 'left')

    full_aliquot_df.set_index('Aliquot', inplace = True, drop = False)
    
    plot_data = samples[samples['Transect'] == transect] if transect else samples.copy()
    full_aliquot_df = full_aliquot_df[full_aliquot_df['Transect'] == transect] if transect else full_aliquot_df

    ##########################
    # Figure Set up -----------------------------
    hist_colors = {'ZHe': ZHeColor,'AHe': AHeColor, 'AFT': AFTColor, 'ZFT': ZFTColor}

    outlier_markers = {'keep': 'o', 'reject': 'X'}
    colors = {'keep': 'black', 'reject': 'gray'}
    mineral_markers = {'AHe': AHeMarker, 'ZHe': ZHeMarker, 'AFT': AFTMarker, 'ZFT': ZFTMarker}

    mineral_styles = {
    'AHe': (AHeMarker, AHeColor, AHeMarkerSize),
    'ZHe': (ZHeMarker, ZHeColor, ZHeMarkerSize),
    'AFT': (AFTMarker, AFTColor, AFTMarkerSize),
    'ZFT': (ZFTMarker, ZFTColor, ZFTMarkerSize),
    }

    # Define the column mapping based on the weightedBy option
    column_mapping = {
        'unweighted': {'Mean': 'Mean', 'StDev': 'StDev', 'SEsd': 'SEsd', 'SEiu': 'SEiu'},
        'inverse_variance': {'Mean': 'Mean_w', 'StDev': 'StDev_w', 'SEsd': 'SEsd_w', 'SEiu': 'SEiu_w'},
        'relative_deviation': {'Mean': 'Mean_rw', 'StDev': 'StDev_rw', 'SEsd': 'SEsd_rw', 'SEiu': 'SEiu_rw'}
    }

    # Get the correct columns based on weightedBy
    selected_columns = column_mapping.get(weightedBy, column_mapping['unweighted']) 

    # # Function to lighten a color (for standard error option)
    # def lighten_color(color, amount=0.7):
    #     # Convert the color to RGB
    #     base_color = mcolors.to_rgba(color)
    #     # Interpolate between the color and white
    #     new_color = [1 - (1 - component) * (1 - amount) for component in base_color[:3]] + [base_color[3]]
    #     return new_color
    
    ## orient the histogram plot correctly to correspond with 'Age' axis
    # Age on the x-axis
    if x_variable == 'Age':
        fig, (ax1,ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [2, 8]}, figsize = (15,8), sharex = True)

        ## Plotting Age Histogram (Axis 1) -----------------------------
        if y_bounds:
            filtered_aliquot_df = full_aliquot_df[(full_aliquot_df[y_variable] >= y_bounds[0]) & 
                                                (full_aliquot_df[y_variable] <= y_bounds[1])]

            sns.histplot(ax = ax1, x = filtered_aliquot_df.Corrected_Date_Ma, hue = filtered_aliquot_df.Mineral, 
                    palette=hist_colors, alpha = 0.8,
                    stat = 'count', binwidth = bin_width, kde = True, legend = True)

        else:
            sns.histplot(ax = ax1, x = full_aliquot_df.Corrected_Date_Ma, hue = full_aliquot_df.Mineral, 
                    palette=hist_colors, alpha = 0.8,
                    stat = 'count', binwidth = bin_width, kde = True, legend = True)

        ## Plotting Sample Means and Aliquot Data (Axis 2) -----------------------------
        if show_aliquots:
            # Create separate datasets for rejected and non-rejected data points
            rejected_df = full_aliquot_df[full_aliquot_df['outlier'] == 'reject']
            non_rejected_df = full_aliquot_df[full_aliquot_df['outlier'] == 'keep']

            # Plot rejected data points with X markers
            ax2.scatter(
                x=rejected_df.Corrected_Date_Ma,
                y=rejected_df[y_variable],
                marker=outlier_markers['reject'],
                color=colors['reject']
            )

            # Plot non-rejected data points with mineral markers
            for mineral, marker in mineral_markers.items():
                df = non_rejected_df[non_rejected_df['Mineral'] == mineral]
                ax2.scatter(
                    x=df.Corrected_Date_Ma,
                    y=df[y_variable],
                    marker=marker,
                    color=colors['keep']
                )

        ### Sample Means w/ Errorbars (Axis 2) 
        for index, row in plot_data.iterrows():
            marker, color, size = mineral_styles.get(row['Mineral'], ('o', 'maroon', 10))
            label = row.Sample

            # Plot Standard Deviation
            ax2.errorbar(x=row[selected_columns['Mean']], y=plot_data.loc[index][y_variable],
                    xerr=row[selected_columns['StDev']], ecolor = color, elinewidth=3,
                    mfc = color, marker = marker, ms = size, mec = 'black', mew = 1
            )

            # Plot Standard Error (if specified)
            if plot_SE:
                # Apply the lighter color for plotting SE
                lighter_color = lighten_color(color, amount=0.7)  # adjust 'amount' for more brightness

                # Select the appropriate xerr column based on SE_basedOn
                if SE_basedOn == 'SD':
                    err_value = row[selected_columns['SEsd']]
                elif SE_basedOn == 'IU':
                    err_value = row[selected_columns['SEiu']]
                elif SE_basedOn == 'max':
                    err_value = max(row[selected_columns['SEsd']], row[selected_columns['SEiu']])
                else:
                    raise ValueError("Invalid option for SE_basedOn. Choose from 'SD', 'IU', or 'max'.")

                ax2.errorbar(
                    row[selected_columns['Mean']], plot_data.loc[index][y_variable],
                    xerr=err_value, ecolor=lighter_color, fmt='none', mew=1, elinewidth=3
                )

            # Sample Labels
            if label_samples:
                if x_bounds and y_bounds:
                    if row[selected_columns['Mean']] < x_bounds[1] and row[selected_columns['Mean']] > x_bounds[0] and plot_data.loc[index][y_variable] < y_bounds[1] and plot_data.loc[index][y_variable] > y_bounds[0]:
                        ax2.text(x = row[selected_columns['Mean']] + label_offset[0],
                                y = plot_data.loc[index][y_variable] + label_offset[1],
                                s = label,
                                fontsize = 12,
                                fontname = "Arial",
                                rotation = 0,
                                va = 'center',
                                ha = 'left',
                                color = GREY60,
                                weight = 'book',
                                style = 'italic')

                elif y_bounds:
                    if row[y_variable] < y_bounds[1] and row[y_variable] > y_bounds[0]:
                        ax2.text(x = row[selected_columns['Mean']] + label_offset[0],
                                y = plot_data.loc[index][y_variable] + label_offset[1],
                                s = label,
                                fontsize = 12,
                                fontname = "Arial",
                                rotation = 0,
                                va = 'center',
                                ha = 'left',
                                color = GREY60,
                                weight = 'book',
                                style = 'italic')
                elif x_bounds:
                    if row[selected_columns['Mean']] < x_bounds[1] and row[selected_columns['Mean']] > x_bounds[0]:
                        ax2.text(x = row[selected_columns['Mean']] + label_offset[0],
                                y = plot_data.loc[index][y_variable] + label_offset[1],
                                s = label,
                                fontsize = 12,
                                fontname = "Arial",
                                rotation = 0,
                                va = 'center',
                                ha = 'left',
                                color = GREY60,
                                weight = 'book',
                                style = 'italic')
                else:
                    ax2.text(x = row[selected_columns['Mean']] + label_offset[0],
                                y = plot_data.loc[index][y_variable] + label_offset[1],
                                s = label,
                                fontsize = 12,
                                fontname = "Arial",
                                rotation = 0,
                                va = 'center',
                                ha = 'left',
                                color = GREY60,
                                weight = 'book',
                                style = 'italic')
        ##########################
        ### Axes and Spine Customization ----------------------------- 
        ax2.spines["left"].set_color('k')
        ax2.spines["bottom"].set_color('k')
        ax2.spines["right"].set_color('k')
        ax2.spines["top"].set_color('k')

        ax2.grid(True, color = GREY91, linestyle = '-', linewidth = 1)

        if y_bounds:
            ax2.set_ylim(y_bounds)

        if x_bounds:
            ax2.set_xlim(x_bounds)
        else: 
            plt.gca().set_xlim(left=0)

        ax2.set_xlabel('Corrected Cooling Age (Ma)',
                        fontname= "Arial",
                        fontsize=12,
                        weight=500,
                        color=GREY40
                    )

        if y_variable == 'Latitude':
            y_axis_label = 'Latitude'
        elif y_variable == 'Longitude':
            y_axis_label = 'Longitude'
        elif y_variable == 'Structural_Level':
            y_axis_label = 'Structural Level'
        elif y_variable == 'Elevation_m':
            y_axis_label = 'Elevation (m)'
        else:
            y_axis_label = 'Y variable'

        ax2.set_ylabel(y_axis_label,
                        fontname= "Arial",
                        fontsize=12,
                        weight=500,
                        color=GREY40
                    )

        ax2.tick_params(axis="x", length=5, color=GREY91)
        ax2.tick_params(axis="y", length=5, color=GREY91)

    ## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Age on the y-axis
    elif y_variable == 'Age':
        fig, (ax1,ax2) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [1, 8]}, figsize = (12,6), sharey = True)

        ## Plotting Age Histogram (Axis 1) 
        if x_bounds:
            filtered_aliquot_df = full_aliquot_df[(full_aliquot_df[x_variable] >= x_bounds[0]) & (full_aliquot_df[x_variable] <= x_bounds[1])]
            
            sns.histplot(ax = ax1, y = filtered_aliquot_df.Corrected_Date_Ma, hue = filtered_aliquot_df.Mineral, 
                    palette=hist_colors, alpha = 0.8,
                    stat = 'count', binwidth = bin_width, kde = True, legend = True)
        
        else:
            sns.histplot(ax = ax1, y = full_aliquot_df.Corrected_Date_Ma, hue = full_aliquot_df.Mineral, 
                    palette=hist_colors, alpha = 0.8,
                    stat = 'count', binwidth = bin_width, kde = True, legend = True)
            
        ## Plotting Sample Means and Aliquot Data (Axis 2) 
        if show_aliquots:
            # Create separate datasets for rejected and non-rejected data points
            rejected_df = full_aliquot_df[full_aliquot_df['outlier'] == 'reject']
            non_rejected_df = full_aliquot_df[full_aliquot_df['outlier'] == 'keep']

            # Plot rejected data points with X markers
            ax2.scatter(
                x=rejected_df[x_variable],
                y=rejected_df.Corrected_Date_Ma,
                marker=outlier_markers['reject'],
                color=colors['reject']
            )

            # Plot non-rejected data points with mineral markers
            for mineral, marker in mineral_markers.items():
                df = non_rejected_df[non_rejected_df['Mineral'] == mineral]
                ax2.scatter(
                    x=df[x_variable],
                    y=df.Corrected_Date_Ma,
                    marker=marker,
                    color=colors['keep']
                )
        
        ### Sample Means w/ Errorbars (Axis 2) 
        for index, row in plot_data.iterrows():
            marker, color, size = mineral_styles.get(row['Mineral'], ('o', 'maroon', 10))
            label = row.Sample
        
            ax2.errorbar(x=plot_data.loc[index][x_variable], y=row[selected_columns['Mean']],
                    yerr=row[selected_columns['StDev']], ecolor = color, elinewidth=3,
                    mfc = color, marker = marker, ms = size, mec = 'black', mew = 1
            )

            # Plot Standard Error (if specified)
            if plot_SE:
                # Apply the lighter color for plotting SE
                lighter_color = lighten_color(color, amount=0.7)  # adjust 'amount' for more brightness
        
                # Select the appropriate xerr column based on SE_basedOn
                if SE_basedOn == 'SD':
                    err_value = row[selected_columns['SEsd']]
                elif SE_basedOn == 'IU':
                    err_value = row[selected_columns['SEiu']]
                elif SE_basedOn == 'max':
                    err_value = max(row[selected_columns['SEsd']], row[selected_columns['SEiu']])
                else:
                    raise ValueError("Invalid option for SE_basedOn. Choose from 'SD', 'IU', or 'max'.")
        
                ax2.errorbar(
                    plot_data.loc[index][x_variable], row[selected_columns['Mean']],
                    yerr=err_value, ecolor=lighter_color, fmt='none', mew=1, elinewidth=3
                )

            # Sample Labels
            if label_samples:
                if x_bounds and y_bounds:
                    if row[selected_columns['Mean']] < y_bounds[1] and row[selected_columns['Mean']] > y_bounds[0] and plot_data.loc[index][x_variable] < x_bounds[1] and plot_data.loc[index][x_variable] > x_bounds[0]:
                        ax2.text(x = plot_data.loc[index][x_variable] + label_offset[0],
                                y = plot_data.loc[index].Mean + label_offset[1],
                                s = label,
                                fontsize = 12,
                                fontname = "Arial",
                                rotation = 90,
                                va = 'center',
                                ha = 'left',
                                color = GREY60,
                                weight = 'book',
                                style = 'italic')
                
                elif y_bounds:
                    if row[selected_columns['Mean']] < y_bounds[1] and row[selected_columns['Mean']] > y_bounds[0]:
                        ax2.text(x = plot_data.loc[index][x_variable] + label_offset[0],
                                y = row[selected_columns['Mean']] + label_offset[1],
                                s = label,
                                fontsize = 12,
                                fontname = "Arial",
                                rotation = 90,
                                va = 'center',
                                ha = 'left',
                                color = GREY60,
                                weight = 'book',
                                style = 'italic')
                elif x_bounds:
                    if row[x_variable] < x_bounds[1] and row[x_variable] > x_bounds[0]:
                        ax2.text(x = plot_data.loc[index][x_variable] + label_offset[0],
                                y = row[selected_columns['Mean']] + label_offset[1],
                                s = label,
                                fontsize = 12,
                                fontname = "Arial",
                                rotation = 90,
                                va = 'center',
                                ha = 'left',
                                color = GREY60,
                                weight = 'book',
                                style = 'italic')
                else:
                    ax2.text(x = plot_data.loc[index][x_variable] + label_offset[0],
                                y = row[selected_columns['Mean']] + label_offset[1],
                                s = label,
                                fontsize = 12,
                                fontname = "Arial",
                                rotation = 90,
                                va = 'center',
                                ha = 'left',
                                color = GREY60,
                                weight = 'book',
                                style = 'italic')
        ##########################
        ### Axes and Spine Customization -----------------------------
        ax1.invert_xaxis()
        
        ax2.spines["left"].set_color('k')
        ax2.spines["bottom"].set_color('k')
        ax2.spines["right"].set_color('k')
        ax2.spines["top"].set_color('k')

        ax2.grid(True, color = GREY91, linestyle = '-', linewidth = 1)

        if y_bounds:
            ax2.set_ylim(y_bounds)
        else: 
            plt.gca().set_ylim(bottom=0)

        if x_bounds:
            ax2.set_xlim(x_bounds)
            
        ax1.set_ylabel('Corrected Cooling Age (Ma)',
                        fontname= "Arial",
                        fontsize=12,
                        weight=500,
                        color=GREY40
                    )
        
        if x_variable == 'Latitude':
            x_axis_label = 'Latitude'
        elif x_variable == 'Longitude':
            x_axis_label = 'Longitude'
        elif x_variable == 'Structural_Level':
            x_axis_label = 'Structural Level'
        elif x_variable == 'Elevation_m':
            x_axis_label = 'Elevation (m)'
        else:
            x_axis_label = 'X variable'
        
        ax2.set_xlabel(x_axis_label,
                        fontname= "Arial",
                        fontsize=12,
                        weight=500,
                        color=GREY40
                    )

        ax2.tick_params(axis="x", length=5, color=GREY91)
        ax2.tick_params(axis="y", length=5, color=GREY91)


    ##########################
    ### Show and Save Figure -----------------------------
    if savefig:
        pathlib.Path(saveFolder).mkdir(parents=True, exist_ok=True)
        if savefigFileName:
            filepath = '{}/{}.pdf'.format(saveFolder, savefigFileName)
        else:
            if transect:
                filepath = '{}/{}_v_Age_{}_transect.pdf'.format(saveFolder, x_variable, transect)
            else:
                filepath = '{}/{}_v_Age.pdf'.format(saveFolder, x_variable)
        plt.savefig(filepath, dpi = 'figure', bbox_inches='tight', pad_inches = 0.5)

    fig.tight_layout()
    plt.show();

def plot_AgeVersus_wZoomRegression(samples: pd.DataFrame, aliquots: pd.DataFrame, transect: Optional[str], y_variable: str,
    inset_xlim: Tuple[float, float], full_x_bounds: Optional[Tuple[float, float]] = None, y_bounds: Optional[Tuple[float, float]] = None,
    inset_x_interval: int = 10, full_x_interval: int = 50,
    weightedBy : str = 'unweighted', plot_SE: bool = False, SE_basedOn: str = 'max',
    label_samples: bool = True, label_offset: Tuple[int, int] = (5, 5),
    AHeColor: str = 'cadetblue', AHeMarker: str = 'h',AHeMarkerSize: int = 12,
    ZHeColor: str = 'darkmagenta',ZHeMarker: str = 'D',ZHeMarkerSize: int = 10,
    AFTColor: str = 'gray',AFTMarker: str = '^',AFTMarkerSize: int = 12,
    ZFTColor: str = 'forestgreen',ZFTMarker: str = 'o',ZFTMarkerSize: int = 12,
    AHe_regression: bool = False, AHeRegressionColor: str = 'cadetblue', excludeAHeOutliers: Optional[List[str]] = None,excludeAHeSamplesRegression: Optional[List[str]] = None, excludeAHeAliquotsRegression: Optional[List[str]] = None,
    ZHe_regression: bool = False, ZHeRegressionColor: str = 'darkmagenta', excludeZHeOutliers: Optional[List[str]] = None,excludeZHeSamplesRegression: Optional[List[str]] = None, excludeZHeAliquotsRegression: Optional[List[str]] = None,
    AFTRegression: bool = False, AFTRegressionColor: str = 'gainsboro', excludeAFTOutliers: Optional[List[str]] = None,excludeAFTSamplesRegression: Optional[List[str]] = None, excludeAFTAliquotsRegression: Optional[List[str]] = None,
    ZFTRegression: bool = False, ZFTRegressionColor: str = 'forestgreen', excludeZFTOutliers: Optional[List[str]] = None,excludeZFTSamplesRegression: Optional[List[str]] = None, excludeZFTAliquotsRegression: Optional[List[str]] = None,
    savefig: bool = False, savefigFileName: Optional[str] = None, saveFolder: str = 'Plots'
    ) -> None:
    """
    Plots Cooling Age against a specified y-axis variable (e.g., Latitude, Longitude, Elevation, Structural Level),
    with an inset zoom and regression options for each mineral type.

    Allows customization of markers, colors, labels, and regression details for different thermochronometers.

    Parameters
    ----------
    samples : pd.DataFrame
        DataFrame containing sample data and summary statistics, indexed by sample name.
    aliquots : pd.DataFrame
        DataFrame containing aliquot data, indexed by aliquot name. Each sample may have multiple aliquots.
    transect : str, optional
        Specific transect to plot. If None, the full dataset is plotted.
    y_variable : str
        The variable to plot on the y-axis (e.g., 'Latitude', 'Longitude', 'Elevation', 'Structural Level').
    inset_xlim : tuple of float
        X-axis limits for the inset plot, representing the zoomed-in view.
    full_x_bounds : tuple of float, optional
        X-axis limits for the full plot.
    y_bounds : tuple of float, optional
        Y-axis limits for the plot.
    inset_x_interval : int, default 10
        Interval for x-ticks in the inset plot.
    full_x_interval : int, default 50
        Interval for x-ticks in the full plot.
    weightedBy : str, default='unweighted'
        Specifies the weighting method, with options such as 'unweighted'.
    plot_SE : bool, default=False
        Whether to plot standard error bars.
    SE_basedOn : str, default='max'
        Basis for SE calculation; 'max' uses the maximum standard error.
    label_samples : bool, default True
        Whether to label sample points on the plot.
    label_offset : tuple of int, default (5, 5)
        Offset for sample labels to position them more clearly.
    AHeColor : str, default 'cadetblue'
        Color for Apatite He (AHe) markers.
    AHeMarker : str, default 'h'
        Marker style for AHe data points.
    AHeMarkerSize : int, default 12
        Marker size for AHe data points.
    ZHeColor : str, default 'darkmagenta'
        Color for Zircon He (ZHe) markers.
    ZHeMarker : str, default 'D'
        Marker style for ZHe data points.
    ZHeMarkerSize : int, default 10
        Marker size for ZHe data points.
    AFTColor : str, default 'gray'
        Color for Apatite Fission Track (AFT) markers.
    AFTMarker : str, default '^'
        Marker style for AFT data points.
    AFTMarkerSize : int, default 12
        Marker size for AFT data points.
    ZFTColor : str, default 'forestgreen'
        Color for Zircon Fission Track (ZFT) markers.
    ZFTMarker : str, default 'o'
        Marker style for ZFT data points.
    ZFTMarkerSize : int, default 12
        Marker size for ZFT data points.
    AHe_regression : bool, default False
        Whether to perform regression on AHe data points.
    AHeRegressionColor : str, default 'cadetblue'
        Color for the AHe regression line and confidence intervals.
    excludeAHeOutliers : list of str, optional
        List of outliers to exclude from AHe regression.
    excludeAHeSamplesRegression : list of str, optional
        List of samples to exclude from AHe regression.
    excludeAHeAliquotsRegression : list of str, optional
        List of aliquots to exclude from AHe regression.
    ZHe_regression : bool, default False
        Whether to perform regression on ZHe data points.
    ZHeRegressionColor : str, default 'darkmagenta'
        Color for the ZHe regression line and confidence intervals.
    excludeZHeOutliers : list of str, optional
        List of outliers to exclude from ZHe regression.
    excludeZHeSamplesRegression : list of str, optional
        List of samples to exclude from ZHe regression.
    excludeZHeAliquotsRegression : list of str, optional
        List of aliquots to exclude from ZHe regression.
    AFTRegression : bool, default False
        Whether to perform regression on AFT data points.
    AFTRegressionColor : str, default 'gainsboro'
        Color for the AFT regression line and confidence intervals.
    excludeAFTOutliers : list of str, optional
        List of outliers to exclude from AFT regression.
    excludeAFTSamplesRegression : list of str, optional
        List of samples to exclude from AFT regression.
    excludeAFTAliquotsRegression : list of str, optional
        List of aliquots to exclude from AFT regression.
    ZFTRegression : bool, default False
        Whether to perform regression on ZFT data points.
    ZFTRegressionColor : str, default 'forestgreen'
        Color for the ZFT regression line and confidence intervals.
    excludeZFTOutliers : list of str, optional
        List of outliers to exclude from ZFT regression.
    excludeZFTSamplesRegression : list of str, optional
        List of samples to exclude from ZFT regression.
    excludeZFTAliquotsRegression : list of str, optional
        List of aliquots to exclude from ZFT regression.
    savefig : bool, default False
        Whether to save the generated plot to a file.
    savefigFileName : str, optional
        Name of the file to save the plot to. Ignored if `savefig` is False.
    saveFolder : str, default 'Plots'
        Folder where the plot will be saved if `savefig` is True.

    Returns
    -------
    None
        Displays the plot; saves it if savefig is True.
    """
    # Merge sample and aliquot data for complete data set
    full_aliquot_df = pd.merge(
        aliquots,
        samples[['Mean', 'StDev', 'Mean_w', 'StDev_w','Mean_rw', 'StDev_rw','Latitude', 'Longitude', 'Elevation_m', 'Structural_Level', 'Transect']],
        on='Sample',
        how='left'
    )
    full_aliquot_df.set_index('Aliquot', inplace=True, drop=False)

    plot_data = samples[samples['Transect'] == transect] if transect else samples.copy()
    full_aliquot_df = full_aliquot_df[full_aliquot_df['Transect'] == transect] if transect else full_aliquot_df

    # Set up figure and axes
    fig, (ax1, ax2) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [1, 4]}, figsize=(12, 6), sharey=True)

    # Regression for each chronometer
    if AHe_regression:
        ahe_regression_df = filter_regression_data(full_aliquot_df, 'AHe', excludeAHeOutliers, excludeAHeSamplesRegression, excludeAHeAliquotsRegression)
        plot_regression(ax1, ahe_regression_df, plot_data, y_variable, AHeRegressionColor, 'AHe', weightedBy)
        plot_regression(ax2, ahe_regression_df, plot_data, y_variable, AHeRegressionColor, 'AHe', weightedBy)

    if ZHe_regression:
        zhe_regression_df = filter_regression_data(full_aliquot_df, 'ZHe', excludeZHeOutliers, excludeZHeSamplesRegression, excludeZHeAliquotsRegression)
        plot_regression(ax1, zhe_regression_df, plot_data, y_variable, ZHeRegressionColor, 'ZHe', weightedBy)
        plot_regression(ax2, zhe_regression_df, plot_data, y_variable, ZHeRegressionColor, 'ZHe', weightedBy)

    if AFTRegression:
        aft_regression_df = filter_regression_data(full_aliquot_df, 'AFT', excludeAFTOutliers, excludeAFTSamplesRegression, excludeAFTAliquotsRegression)
        plot_regression(ax1, aft_regression_df, plot_data, y_variable, AFTRegressionColor, 'AFT', weightedBy)
        plot_regression(ax2, aft_regression_df, plot_data, y_variable, AFTRegressionColor, 'AFT', weightedBy)

    if ZFTRegression:
        zft_regression_df = filter_regression_data(full_aliquot_df, 'ZFT', excludeZFTOutliers, excludeZFTSamplesRegression, excludeZFTAliquotsRegression)
        plot_regression(ax1, zft_regression_df, plot_data, y_variable, ZFTRegressionColor, 'ZFT', weightedBy)
        plot_regression(ax2, zft_regression_df, plot_data, y_variable, ZFTRegressionColor, 'ZFT', weightedBy)

    # Plot data points
    outlier_markers = {'keep': 'o', 'reject': 'X'}
    colors = {'keep': 'black', 'reject': 'gray'}
    mineral_markers = {'AHe': AHeMarker, 'ZHe': ZHeMarker, 'AFT': AFTMarker, 'ZFT': ZFTMarker}
    
    for ax in [ax1, ax2]:
        for outlier, marker in outlier_markers.items():
            df = full_aliquot_df[full_aliquot_df['outlier'] == outlier]
            ax.scatter(df['Corrected_Date_Ma'], df[y_variable], marker=marker, color=colors[outlier])

    # Plot sample means with error bars
    mineral_styles = {
        'AHe': (AHeMarker, AHeColor, AHeMarkerSize),
        'ZHe': (ZHeMarker, ZHeColor, ZHeMarkerSize),
        'AFT': (AFTMarker, AFTColor, AFTMarkerSize),
        'ZFT': (ZFTMarker, ZFTColor, ZFTMarkerSize),
    }

    # Define the column mapping based on the weightedBy option
    column_mapping = {
        'unweighted': {'Mean': 'Mean', 'StDev': 'StDev', 'SEsd': 'SEsd', 'SEiu': 'SEiu'},
        'inverse_variance': {'Mean': 'Mean_w', 'StDev': 'StDev_w', 'SEsd': 'SEsd_w', 'SEiu': 'SEiu_w'},
        'relative_deviation': {'Mean': 'Mean_rw', 'StDev': 'StDev_rw', 'SEsd': 'SEsd_rw', 'SEiu': 'SEiu_rw'}
    }

    # Get the correct columns based on weightedBy
    selected_columns = column_mapping.get(weightedBy, column_mapping['unweighted']) 

    for index, row in plot_data.iterrows():
        marker, color, size = mineral_styles.get(row['Mineral'], ('o', 'maroon', 10))
        label = row['Sample']
        for ax in [ax1, ax2]:
            # ax.errorbar(row['Mean'], row[y_variable], xerr=row['StDev'], ecolor=color, mfc=color, marker=marker, ms=size, mec='black', mew=1, elinewidth = 3)

            # Use the selected columns for x-value and xerr
            ax.errorbar(
                row[selected_columns['Mean']], row[y_variable],
                xerr=row[selected_columns['StDev']], ecolor=color,
                mfc=color, marker=marker, ms=size, mec='black', mew=1, elinewidth=3
            )
            
            if plot_SE:
                # Apply the lighter color for plotting SE
                lighter_color = lighten_color(color, amount=0.7)  # adjust 'amount' for more brightness

                # Select the appropriate xerr column based on SE_basedOn
                if SE_basedOn == 'SD':
                    xerr_value = row[selected_columns['SEsd']]
                elif SE_basedOn == 'IU':
                    xerr_value = row[selected_columns['SEiu']]
                elif SE_basedOn == 'max':
                    xerr_value = max(row[selected_columns['SEsd']], row[selected_columns['SEiu']])
                else:
                    raise ValueError("Invalid option for SE_basedOn. Choose from 'SD', 'IU', or 'max'.")

                
                #ax.errorbar(row['Mean'], row[y_variable], xerr=xerr_value, ecolor=lighter_color, fmt='none', mew=1, elinewidth=3)

                ax.errorbar(
                    row[selected_columns['Mean']], row[y_variable],
                    xerr=xerr_value, ecolor=lighter_color, fmt='none', mew=1, elinewidth=3
                )

        ## Sample Labels
        if label_samples:
            
            # setting alignment for labels in the inset
            if plot_data.loc[index]['Mean'] < (inset_xlim/2):
                ha = 'left'
            elif plot_data.loc[index]['Mean'] > (inset_xlim/2):
                ha = 'right'

            if plot_data.loc[index]['Mineral'] == 'ZHe':
                ha = 'left'

            # labels for inset plot
            if y_bounds:
                if row[y_variable] <= y_bounds[1] and row[y_variable] >= y_bounds[0] and plot_data.loc[index]['Mean'] < inset_xlim:
                    ax1.text(x=plot_data.loc[index]['Mean'] + label_offset[0],
                            y=plot_data.loc[index][y_variable] + label_offset[1],
                            s = label,
                            fontsize = 12,
                            fontname = "Arial",
                            va = 'center',
                            ha = ha,
                            color = GREY60,
                            weight = 'book',
                            style = 'italic')
            
            elif plot_data.loc[index]['Mean'] < inset_xlim:
                ax1.text(x=plot_data.loc[index]['Mean'] + label_offset[0],
                         y=plot_data.loc[index][y_variable] + label_offset[1],
                         s=label,
                         fontsize = 12,
                         fontname = "Arial",
                         va = 'center',
                         ha = ha,
                         color = GREY60,
                         weight = 'book',
                         style = 'italic')
            
            # labels for full plot
            if full_x_bounds and y_bounds:
                if row[y_variable] <= y_bounds[1] and row[y_variable] >= y_bounds[0] and row.Mean >= full_x_bounds[1] and row.Mean <= full_x_bounds[0]:
                    ax2.text(x=plot_data.loc[index]['Mean'] + label_offset[0],
                            y=plot_data.loc[index][y_variable] + label_offset[1],
                            s = label,
                            fontsize = 12,
                            fontname = "Arial",
                            va = 'center',
                            ha = 'left',
                            color = GREY60,
                            weight = 'book',
                            style = 'italic')
                        
            elif full_x_bounds:
                if row.Mean >= full_x_bounds[1] and row.Mean <= full_x_bounds[0]:
                    ax2.text(x=plot_data.loc[index]['Mean'] + label_offset[0],
                            y=plot_data.loc[index][y_variable] + label_offset[1],
                            s = label,
                            fontsize = 12,
                            fontname = "Arial",
                            va = 'center',
                            ha = 'left',
                            color = GREY60,
                            weight = 'book',
                            style = 'italic')
            elif y_bounds:
                if row[y_variable] <= y_bounds[1] and row[y_variable] >= y_bounds[0]:
                    ax2.text(x=plot_data.loc[index]['Mean'] + label_offset[0],
                            y=plot_data.loc[index][y_variable] + label_offset[1],
                            s = label,
                            fontsize = 12,
                            fontname = "Arial",
                            va = 'center',
                            ha = 'left',
                            color = GREY60,
                            weight = 'book',
                            style = 'italic')
            else:    
                ax2.text(x=plot_data.loc[index]['Mean'] + label_offset[0],
                         y=plot_data.loc[index][y_variable] + label_offset[1],
                         s=label,
                         fontsize = 12,
                         fontname = "Arial",
                         va = 'center',
                         ha = 'left',
                         color = GREY60,
                         weight = 'book',
                         style = 'italic')

    ### Axes and Spine Customization -----------------------------
    ## Inset
    ax1.spines["left"].set_color('k')
    ax1.spines["bottom"].set_color('k')
    ax1.spines["right"].set_color('k')
    ax1.spines["top"].set_color('k')

    ax1.grid(True, color = GREY91, linestyle = '-', linewidth = 1)
    ax1.set_xlim(0,inset_xlim) 

    if y_bounds:
        ax1.set_ylim(y_bounds)

    ax1.set_xticks([x for x in np.arange(0, inset_xlim + 10, inset_x_interval)])
    ax1.set_xticklabels(
        [int(x) for x in np.arange(0, inset_xlim + 10, inset_x_interval)],
        fontname= "Arial",
        fontsize=10,
        weight=500,
        color=GREY40
    )

    ax1.set_xlabel('Corrected Cooling Age (Ma)',
                    fontname= "Arial",
                    fontsize=12,
                    weight=500,
                    color=GREY40
                  )
        
    if y_variable == 'Latitude':
        y_axis_label = 'Latitude'
    elif y_variable == 'Longitude':
        y_axis_label = 'Longitude'
    elif y_variable == 'Structural_Level':
        y_axis_label = 'Structural Level'
    elif y_variable == 'Elevation_m':
        y_axis_label = 'Elevation (m)'
    else:
        y_axis_label = y_variable
    
    ax1.set_ylabel(y_axis_label,
                    fontname= "Arial",
                    fontsize=12,
                    weight=500,
                    color=GREY40
                  )

    ax1.tick_params(axis="x", length=5, color=GREY91)
    ax1.tick_params(axis="y", length=5, color=GREY91)

    ## All Data
    ax2.spines["left"].set_color('k')
    ax2.spines["bottom"].set_color('k')
    ax2.spines["right"].set_color('k')
    ax2.spines["top"].set_color('k')
    
    if full_x_bounds:
        ax2.set_xlim(full_x_bounds)
        
        ax2.grid(True, color = GREY91, linestyle = '-', linewidth = 1)
        ax2.set_xlim(full_x_bounds)
        ax2.set_xticks([x for x in np.arange(0, full_x_bounds[1], full_x_interval)])
        ax2.set_xticklabels(
            [int(x) for x in np.arange(0, full_x_bounds[1], full_x_interval)],
            fontname= "Arial",
            fontsize=10,
            weight=500,
            color=GREY40
        )
        
    else:
        max_age = max(full_aliquot_df.Corrected_Date_Ma) + 50

        ax2.grid(True, color = GREY91, linestyle = '-', linewidth = 1)
        ax2.set_xlim(0,max_age)
        ax2.set_xticks([x for x in np.arange(0, max_age, full_x_interval)])
        ax2.set_xticklabels(
            [int(x) for x in np.arange(0, max_age, full_x_interval)],
            fontname= "Arial",
            fontsize=10,
            weight=500,
            color=GREY40
        )

    ax2.set_xlabel('Corrected Cooling Age (Ma)',
                    fontname= "Arial",
                    fontsize=12,
                    weight=500,
                    color=GREY40
                  )

    ax2.set_ylabel(y_axis_label,
                    fontname= "Arial",
                    fontsize=12,
                    weight=500,
                    color=GREY40
                  )

    ax2.tick_params(axis="x", length=5, color=GREY91)
    ax2.tick_params(axis="y", length=0, color=GREY91)

    ## Transect Label
    props = dict(boxstyle='square', facecolor='white', alpha=0.8)
    
    if transect:
        ax2.text(0.92,0.92, transect,
                 transform=ax2.transAxes,
                 fontsize = 18,
                 fontname = "Arial",
                 color = 'k',
                 va = 'center',
                 ha = 'center',
                 bbox = props,
                 weight = 'bold')
    
    if savefig:
        pathlib.Path(saveFolder).mkdir(parents=True, exist_ok=True)
        filepath = '{}/{}.pdf'.format(saveFolder, savefigFileName)
        plt.savefig(filepath, dpi = 'figure', bbox_inches='tight', pad_inches = 0.5)

    plt.show()

def plot_AgeVersus_wZoomRegressionHistogram(samples: pd.DataFrame, aliquots: pd.DataFrame, transect: Optional[str], y_variable: str,
    inset_xlim: Tuple[float, float], x_interval_inset: int = 10, full_x_interval: int = 50, 
    full_x_bounds: Optional[Tuple[float, float]] = None, y_bounds: Optional[Tuple[float, float]] = None,
    weightedBy: str = 'unweighted', plot_SE: bool = False, SE_basedOn: str = 'max',
    insetBinWidth: int = 2, fullBinWidth: int = 5,
    stat: str = 'count', kde: bool = True, histLegend: bool = True,
    AHeColor: str = 'cornflowerblue', AHeMarker: str = 'h', AHeMarkerSize: int = 10,
    ZHeColor: str = 'firebrick', ZHeMarker: str = 'D', ZHeMarkerSize: int = 8,
    AFTColor: str = 'gray', AFTMarker: str = '^', AFTMarkerSize: int = 12,
    ZFTColor: str = 'forestgreen', ZFTMarker: str = 'o', ZFTMarkerSize: int = 12,
    AHe_regression: bool = False, AHeRegressionColor: str = 'lightsteelblue', excludeAHeSamples: Optional[List[str]] = None,
    ZHe_regression: bool = False, ZHeRegressionColor: str = 'thistle', excludeZHeSamples: Optional[List[str]] = None,
    AFTRegression: bool = False, AFTRegressionColor: str = 'gainsboro', excludeAFTSamples: Optional[List[str]] = None,
    ZFTRegression: bool = False, ZFTRegressionColor: str = 'lightgreen', excludeZFTSamples: Optional[List[str]] = None,
    label_transects: bool = False, separateZrLabels: bool = False, plotDepoAges: bool = False,
    savefig: bool = False, savefigFileName: Optional[str] = None, saveFolder: str = 'Plots'
    ) -> None:
    """
    Plots Cooling Age against a specified y-axis variable (e.g., Latitude, Longitude, Elevation, Structural Level),
    with an inset zoom, regression options for each mineral type, and a histogram of aliquot cooling ages for each mineral type

    Allows customization of markers, colors, labels, and regression details for different thermochronometers.

    Parameters
    ----------
    samples : pd.DataFrame
        DataFrame containing sample data and summary statistics, indexed by sample name.
    aliquots : pd.DataFrame
        DataFrame containing aliquot data, indexed by aliquot name. Each sample may have multiple aliquots.
    transect : str, optional
        Specific transect to plot. If None, the full dataset is plotted.
    y_variable : str
        The variable to plot on the y-axis (e.g., 'Latitude', 'Longitude', 'Elevation', 'Structural Level').
    inset_xlim : tuple of float
        X-axis limits for the inset plot, representing the zoomed-in view.
    x_interval_inset : int, default=10
        Tick interval for the x-axis of the inset plot.
    full_x_interval : int, default=50
        Tick interval for the x-axis of the full plot.
    full_x_bounds : tuple of float, optional
        x-axis limits for the full plot. If None, defaults to auto-scaling.
    y_bounds : tuple of float, optional
        y-axis limits for the plot. If None, defaults to auto-scaling.
    weightedBy : str, default='unweighted'
        Specifies the weighting method, with options such as 'unweighted'.
    plot_SE : bool, default=False
        Whether to plot standard error bars.
    SE_basedOn : str, default='max'
        Basis for SE calculation; 'max' uses the maximum standard error.
    insetBinWidth : int, default=2
        Bin width for the histogram in the inset plot.
    fullBinWidth : int, default=5
        Bin width for the histogram in the full plot.
    stat : str, default='count'
        Statistic to display in the histogram; options include 'count' or 'density'.
    kde : bool, default=True
        Whether to overlay a kernel density estimate on the histogram.
    histLegend : bool, default=True
        Whether to display a legend for the histogram.
    AHeColor : str, default='cornflowerblue'
        Color for AHe markers.
    AHeMarker : str, default='h'
        Marker style for AHe data points.
    AHeMarkerSize : int, default=10
        Size for AHe markers.
    ZHeColor : str, default='firebrick'
        Color for ZHe markers.
    ZHeMarker : str, default='D'
        Marker style for ZHe data points.
    ZHeMarkerSize : int, default=8
        Size for ZHe markers.
    AFTColor : str, default='gray'
        Color for AFT markers.
    AFTMarker : str, default='^'
        Marker style for AFT data points.
    AFTMarkerSize : int, default=12
        Size for AFT markers.
    ZFTColor : str, default='forestgreen'
        Color for ZFT markers.
    ZFTMarker : str, default='o'
        Marker style for ZFT data points.
    ZFTMarkerSize : int, default=12
        Size for ZFT markers.
    AHe_regression : bool, default=False
        Whether to plot regression for AHe data.
    AHeRegressionColor : str, default='lightsteelblue'
        Color for the AHe regression line.
    excludeAHeSamples : list of str, optional
        Samples to exclude from AHe regression.
    ZHe_regression : bool, default=False
        Whether to plot regression for ZHe data.
    ZHeRegressionColor : str, default='thistle'
        Color for the ZHe regression line.
    excludeZHeSamples : list of str, optional
        Samples to exclude from ZHe regression.
    AFTRegression : bool, default=False
        Whether to plot regression for AFT data.
    AFTRegressionColor : str, default='gainsboro'
        Color for the AFT regression line.
    excludeAFTSamples : list of str, optional
        Samples to exclude from AFT regression.
    ZFTRegression : bool, default=False
        Whether to plot regression for ZFT data.
    ZFTRegressionColor : str, default='lightgreen'
        Color for the ZFT regression line.
    excludeZFTSamples : list of str, optional
        Samples to exclude from ZFT regression.
    label_transects : bool, default=False
        Whether to label transects in the plot.
    separateZrLabels : bool, default=False
        Whether to separate labels for Zr.
    plotDepoAges : bool, default=False
        Whether to plot depositional ages.
    savefig : bool, default=False
        Whether to save the plot as a file.
    savefigFileName : str, optional
        Filename for saving the plot.
    saveFolder : str, default='Plots'
        Directory to save the plot file if savefig is True.

    Returns
    -------
    None
        Displays the plot; saves it if savefig is True.
    """
    # Merge sample and aliquot data for complete data set
    full_aliquot_df = pd.merge(
        aliquots,
        samples[['Mean', 'StDev', 'Mean_w', 'StDev_w','Mean_rw', 'StDev_rw','Latitude', 'Longitude', 'Elevation_m', 'Structural_Level', 'Transect']],
        on='Sample',
        how='left'
    )

    # Set up Dataset
    full_aliquot_df.set_index('Aliquot', inplace=True, drop=False)

    plot_data = samples[samples['Transect'] == transect] if transect else samples.copy()
    full_aliquot_df = full_aliquot_df[full_aliquot_df['Transect'] == transect] if transect else full_aliquot_df
    
    # Figure Set Up 
    fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2, 2, gridspec_kw={'width_ratios': [1, 4], 'height_ratios':[2,5]}, 
                                          figsize = (12,8), sharex = 'col', sharey = 'row')
    
    # Histogram of cooling ages
    hist_colors = {'ZHe': ZHeColor,'AHe': AHeColor, 'AFT': AFTColor, 'ZFT': ZFTColor}
    
    # Inset histogram
    sns.histplot(ax = ax1, x = full_aliquot_df.Corrected_Date_Ma, hue = full_aliquot_df.Mineral, 
                palette=hist_colors, alpha = 0.8,
                stat = stat, binwidth = insetBinWidth, kde = kde, legend = histLegend)

    # Full histogram
    sns.histplot(ax = ax2, x = full_aliquot_df.Corrected_Date_Ma, hue = full_aliquot_df.Mineral, 
                palette=hist_colors, alpha = 0.8,
                stat = stat, binwidth = fullBinWidth, kde = kde, legend = histLegend)
    
    # Regression for each chronometer
    if AHe_regression:
        ahe_regression_df = filter_regression_data(full_aliquot_df, 'AHe', True, excludeAHeSamples, None)
        plot_regression(ax3, ahe_regression_df, plot_data, y_variable, AHeRegressionColor, 'AHe', weightedBy)
        plot_regression(ax4, ahe_regression_df, plot_data, y_variable, AHeRegressionColor, 'AHe', weightedBy)

    if ZHe_regression:
        zhe_regression_df = filter_regression_data(full_aliquot_df, 'ZHe', True, excludeZHeSamples, None)
        plot_regression(ax3, zhe_regression_df, plot_data, y_variable, ZHeRegressionColor, 'ZHe', weightedBy)
        plot_regression(ax4, zhe_regression_df, plot_data, y_variable, ZHeRegressionColor, 'ZHe', weightedBy)

    if AFTRegression:
        aft_regression_df = filter_regression_data(full_aliquot_df, 'AFT', True, excludeAFTSamples, None)
        plot_regression(ax3, aft_regression_df, plot_data, y_variable, AFTRegressionColor, 'AFT', weightedBy)
        plot_regression(ax4, aft_regression_df, plot_data, y_variable, AFTRegressionColor, 'AFT', weightedBy)

    if ZFTRegression:
        zft_regression_df = filter_regression_data(full_aliquot_df, 'ZFT', True, excludeZFTSamples, None)
        plot_regression(ax3, zft_regression_df, plot_data, y_variable, ZFTRegressionColor, 'ZFT', weightedBy)
        plot_regression(ax4, zft_regression_df, plot_data, y_variable, ZFTRegressionColor, 'ZFT', weightedBy)

    # Define marker and color mappings
    outlier_markers = {'keep':'o','reject':'X'}
    colors = {'keep':'black', 'reject':'gray'}
    mineral_markers = {'AHe': AHeMarker, 'ZHe': ZHeMarker, 'AFT': AFTMarker, 'ZFT': ZFTMarker}
    
    # Separate data into rejected and non-rejected datasets
    rejected_df = full_aliquot_df[full_aliquot_df['outlier'] == 'reject']
    non_rejected_df = full_aliquot_df[full_aliquot_df['outlier'] == 'keep']
    
    # Plot data for each axis
    for axs in [ax3, ax4]:
         # Plot rejected data points
        axs.scatter(
            x=rejected_df.Corrected_Date_Ma,
            y=rejected_df[y_variable],
            marker=outlier_markers['reject'],
            color=colors['reject']
        )

        # Plot non-rejected data points by mineral type
        for mineral, marker in mineral_markers.items():
            mineral_df = non_rejected_df[non_rejected_df['Mineral'] == mineral]
            if not mineral_df.empty:
                axs.scatter(
                    x=mineral_df.Corrected_Date_Ma,
                    y=mineral_df[y_variable],
                    marker=marker,
                    color=colors['keep']
                )

    # Markers for sample means with error bars
    mineral_styles = {
        'AHe': (AHeMarker, AHeColor, AHeMarkerSize),
        'ZHe': (ZHeMarker, ZHeColor, ZHeMarkerSize),
        'AFT': (AFTMarker, AFTColor, AFTMarkerSize),
        'ZFT': (ZFTMarker, ZFTColor, ZFTMarkerSize),
    }

    # Define the column mapping based on the weightedBy option
    column_mapping = {
        'unweighted': {'Mean': 'Mean', 'StDev': 'StDev', 'SEsd': 'SEsd', 'SEiu': 'SEiu'},
        'inverse_variance': {'Mean': 'Mean_w', 'StDev': 'StDev_w', 'SEsd': 'SEsd_w', 'SEiu': 'SEiu_w'},
        'relative_deviation': {'Mean': 'Mean_rw', 'StDev': 'StDev_rw', 'SEsd': 'SEsd_rw', 'SEiu': 'SEiu_rw'}
    }

    # Get the correct columns based on weightedBy
    selected_columns = column_mapping.get(weightedBy, column_mapping['unweighted']) 

    # Plot Sample Means w/ Errorbars 
    for index, row in plot_data.iterrows():
        marker, color, size = mineral_styles.get(row['Mineral'], ('o', 'maroon', 10))
        label = row.Sample
        
        for ax in [ax3,ax4]:
            # Use the selected columns for x-value and xerr
            ax.errorbar(
                row[selected_columns['Mean']], row[y_variable],
                xerr=row[selected_columns['StDev']], ecolor=color,
                mfc=color, marker=marker, ms=size, mec='black', mew=1, elinewidth=3
            )

            if plot_SE:
                # Apply the lighter color for plotting SE
                lighter_color = lighten_color(color, amount=0.7)  # adjust 'amount' for more brightness

                # Select the appropriate xerr column based on SE_basedOn
                if SE_basedOn == 'SD':
                    xerr_value = row[selected_columns['SEsd']]
                elif SE_basedOn == 'IU':
                    xerr_value = row[selected_columns['SEiu']]
                elif SE_basedOn == 'max':
                    xerr_value = max(row[selected_columns['SEsd']], row[selected_columns['SEiu']])
                else:
                    raise ValueError("Invalid option for SE_basedOn. Choose from 'SD', 'IU', or 'max'.")

                ax.errorbar(
                    row[selected_columns['Mean']], row[y_variable],
                    xerr=xerr_value, ecolor=lighter_color, fmt='none', mew=1, elinewidth=3
                )

    # Labels
    if label_transects:
        
        label_df = plot_data.copy()
        position, place = 'Mean','median' # getting values to base label placement on

        # if labeling zircon transects separately than apatite transects
        if separateZrLabels:
            for index, row in label_df.iterrows():
                if row['Mineral'] == 'ZHe':
                    label_df.loc[index, 'Transect'] = 'z' + row['Transect']

            position, place = 'Q3','max'

        label_df = label_df.groupby('Transect').agg({y_variable:'mean',position:place})

        for index, row in label_df.iterrows():
            label = index

            if label_df.loc[index][position] < inset_xlim:
                ax3.text(x=label_df.loc[index][position] + 5,
                         y=label_df.loc[index][y_variable] + 0.1,
                         s=label,
                         fontsize = 12,
                         fontname = "Arial",
                         va = 'center',
                         ha = 'center',
                         color = 'k',
                         weight = 'bold',
                         style = 'italic')

            ax4.text(x=label_df.loc[index][position] + 20,
                     y=label_df.loc[index][y_variable] + 0.1,
                     s=label,
                     fontsize = 12,
                     fontname = "Arial",
                     va = 'center',
                     ha = 'center',
                     color = 'k',
                     weight = 'bold',
                     style = 'italic')
            
    ### Depositional Ages
    if plotDepoAges:
        
        depoAge_df = plot_data[['Sample','Transect','Mineral','Latitude', 'Longitude', 'Elevation_m', 'Structural_Level',
                              'Depositional_Age_LB','Depositional_Age_UB']]
        depoAge_df = depoAge_df.groupby(['Transect','Depositional_Age_LB',
                                         'Depositional_Age_UB']).agg({'Latitude':['min','max'],
                                                                     'Longitude':['min','max'],
                                                                     'Elevation_m':['min','max'],
                                                                     'Structural_Level':['min','max']}).reset_index()
        depoAge_df.columns = [f'{i}_{j}' for i, j in depoAge_df.columns]
        
        # CURRENTLY ONLY WORKS FOR LATITUDE - ADD ANOTHER OPTION FOR ELEVATION
        for index, row in depoAge_df.iterrows():
            width = row['Depositional_Age_UB_'] - row['Depositional_Age_LB_']

            if y_variable == 'Latitude':
                height = row['Latitude_max'] - row['Latitude_min']
                anchor_y = row['Latitude_min']
                
                # for transects that span small Latitude range
                if height < 0.15:
                    height = 0.15

            elif y_variable == 'Elevation_m':
                height = row['Elevation_m_max'] - row['Elevation_m_min']
                anchor_y = row['Elevation_m_min']
                
                # for transects that span small Elevation range
                if transect:
                    if height < 20:
                        height = 20
                else:
                    if height < 50:
                        height = 50
            
            elif y_variable == 'Structural_Level':
                height = row['Structural_Level_max'] - row['Structural_Level_min']
                anchor_y = row['Structural_Level_min']
                
                # for transects that span small Structural_Level range
                if transect:
                    if height < 20:
                        height = 20
                else:
                    if height < 50:
                        height = 50
            
        
            anchor_x = row['Depositional_Age_LB_']
            
            ax4.add_patch(Rectangle((anchor_x,anchor_y),width,height,
                                 edgecolor = 'k',
                                 facecolor = 'whitesmoke',
                                 fill=True,
                                 alpha = 0.8,
                                 hatch = '//',
                                 zorder = 1000))

    ### AXES AND SPINE CUSTOMIZATION -----------------------------
    ### Zoom In Histogram
    ax1.spines["left"].set_color('k')
    ax1.spines["bottom"].set_color('k')
    ax1.spines["right"].set_color('k')
    ax1.spines["top"].set_color('k')

    ax1.grid(True, axis = 'x', color = GREY91, linestyle = '-', linewidth = 1)

    ax1.set_xlim(0,inset_xlim)

    ax1.set_ylabel('Number of Grains',
                    fontname= "Arial",
                    fontsize=12,
                    weight=500,
                    color=GREY40
                  )

    ax1.tick_params(axis="y", length=5, color=GREY91)

    ### All Data Histogram
    ax2.spines["left"].set_color('k')
    ax2.spines["bottom"].set_color('k')
    ax2.spines["right"].set_color('k')
    ax2.spines["top"].set_color('k')

    ax2.grid(True, axis = 'x', color = GREY91, linestyle = '-', linewidth = 1)

    ### Zoom In Scatter Plot
    ax3.spines["left"].set_color('k')
    ax3.spines["bottom"].set_color('k')
    ax3.spines["right"].set_color('k')
    ax3.spines["top"].set_color('k')

    ax3.grid(True, color = GREY91, linestyle = '-', linewidth = 1)
    ax3.set_xlim(0,inset_xlim)

    ax3.set_xticks([x for x in np.arange(0, inset_xlim + 10, x_interval_inset)])
    ax3.set_xticklabels(
        [int(x) for x in np.arange(0, inset_xlim + 10, x_interval_inset)],
        fontname= "Arial",
        fontsize=10,
        weight=500,
        color=GREY40
    )

    ax3.set_xlabel('Corrected Cooling Age (Ma)',
                    fontname= "Arial",
                    fontsize=12,
                    weight=500,
                    color=GREY40
                  )
    
    y_labels = {
        'Latitude': 'Latitude (º)',
        'Longitude': 'Longitude (º)',
        'Elevation_m': 'Elevation (m)',
        'Structural_Level': 'Structural Level'
        }

    # Get the y_label based on y_variable, defaulting to y_variable if not found in y_labels
    y_label = y_labels.get(y_variable, y_variable)

    ax3.set_ylabel(
        y_label,
        fontname="Arial",
        fontsize=12,
        weight=500,
        color=GREY40
    )

    ax3.tick_params(axis="x", length=5, color=GREY91)
    ax3.tick_params(axis="y", length=5, color=GREY91)

    ### All Data Scatter Plot
    ax4.spines["left"].set_color('k')
    ax4.spines["bottom"].set_color('k')
    ax4.spines["right"].set_color('k')
    ax4.spines["top"].set_color('k')

    max_age = max(full_aliquot_df.Corrected_Date_Ma) + 50

    ax4.grid(True, color = GREY91, linestyle = '-', linewidth = 1)

    if y_bounds:
        ax3.set_ylim(y_bounds)
        ax4.set_ylim(y_bounds)

    if full_x_bounds:
        ax4.set_xlim(full_x_bounds)
        
        ax4.grid(True, color = GREY91, linestyle = '-', linewidth = 1)
        ax4.set_xlim(full_x_bounds)
        ax4.set_xticks([x for x in np.arange(0, full_x_bounds[1], full_x_interval)])
        ax4.set_xticklabels(
            [int(x) for x in np.arange(0, full_x_bounds[1], full_x_interval)],
            fontname= "Arial",
            fontsize=10,
            weight=500,
            color=GREY40
        )
    else:
        ax4.set_xlim(0,max_age)
        ax4.set_xticks([x for x in np.arange(0, max_age, full_x_interval)])
        ax4.set_xticklabels(
            [int(x) for x in np.arange(0, max_age, full_x_interval)],
            fontname= "Arial",
            fontsize=10,
            weight=500,
            color=GREY40
        )

    ax4.set_xlabel('Corrected Cooling Age (Ma)',
                    fontname= "Arial",
                    fontsize=12,
                    weight=500,
                    color=GREY40
                  )
                  
    ax4.tick_params(axis="x", length=5, color=GREY91)
    ax4.tick_params(axis="y", length=0, color=GREY91)
    
    if savefig:
        pathlib.Path(saveFolder).mkdir(parents=True, exist_ok=True)
        filepath = '{}/{}.pdf'.format(saveFolder, savefigFileName)
        plt.savefig(filepath, dpi = 'figure', bbox_inches='tight', pad_inches = 0.5)
            
    plt.show();
