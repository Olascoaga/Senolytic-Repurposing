# -*- coding: utf-8 -*-
"""
Created on Sat Aug  8 19:04:31 2020

@author: olask
"""

import glob
import pandas as pd
from rdkit.Chem import PandasTools
from rdkit.Chem import Descriptors
import os
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import shapiro
from statsmodels.stats.diagnostic import lilliefors
from scipy.stats import kruskal
from scikit_posthocs import posthoc_dunn as dunnett

sns.set_context("paper", font_scale=1.5)
os.makedirs('results', exist_ok=True)
os.makedirs('graphs', exist_ok=True)

csv_files = glob.glob('*.csv') # Especificamos un patron de busqueda de archivos
csv_names = [] # Variable en donde guardaremos el nombre de los fingerprints 
for filename in csv_files:
    csv_names.append(filename.split(".csv")[0]) # Separamos el nombre de la extensión .csv

dfs = [] # Lista vacia en donde guardaremos los df
for filename in csv_files: # Recorremos los archivos y guardamos cada uno en la lista dfs[]
    data = pd.read_csv(filename)
    dfs.append(data)
del csv_files
del data

def molecular_descriptors(dfs):
    # Calculo de los descriptores moleculares
    for df in dfs:
        PandasTools.AddMoleculeColumnToFrame(df, smilesCol='SMILES', molCol='mol', includeFingerprints=True)
        df['MW'] = df['mol'].map(Descriptors.MolWt) # molecular weight
        df['RBs'] = df['mol'].map(Descriptors.NumRotatableBonds) # number of rotatable bonds
        df['HBAs'] = df['mol'].map(Descriptors.NumHAcceptors) # hydrogen bond acceptors
        df['HBDs'] = df['mol'].map(Descriptors.NumHDonors) # hydrogen bond donors
        df['TPSA'] = df['mol'].map(Descriptors.TPSA) # topological polar surface area
        df['logP'] = df['mol'].map(Descriptors.MolLogP) # octanol/water partition coefficient
    return dfs


# Guardando los descriptores en CVS
molecular_descriptors(dfs)
for i in range(len(dfs)):
    dfs[i].loc[:, dfs[i].columns != 'mol'].to_csv( r'results\\' + csv_names[i] + '.csv' )
    
descriptors_labels = list(dfs[0].columns) # Obtain the name of the descriptor from dataframe
del descriptors_labels[0:3] # del first 3 columns (id, smiles and mol)
descriptors = []
for descriptor in descriptors_labels:
    descriptors.append([]) # create a list of list for each descriptor
    
for df in dfs: # Iterating over each df in dfs list 
    i = 0
    for column in df: #Iterating over each column in df 
        if column in descriptors_labels:
            descriptors[i].append(df[column]) # Append to a descriptor in each descriptor position
            i += 1
        else:
            pass
      
for i in range(len(descriptors)):
    descriptor = pd.DataFrame(descriptors[i]).transpose()
    descriptor.columns = csv_names
    descriptor = descriptor.melt(var_name='DB', value_name=descriptors_labels[i])
    # Violinplots
    plt.figure()
    b = sns.violinplot(x="DB", y=descriptors_labels[i], data=descriptor)
    b.set_ylabel(descriptors_labels[i], fontsize=20)
    b.set_xlabel('')
    sns.despine()
    fig = b.get_figure()
    fig.savefig( r'graphs\\' + 'violinplot_' + descriptors_labels[i] + '.png' )
    # Boxplots
    plt.figure()
    plt.boxplot(descriptors[i], labels=csv_names, showfliers=True)
    plt.ylabel(descriptors_labels[i], fontsize=20)
    sns.despine()
    plt.savefig(r'graphs\\' + 'boxplot_' + descriptors_labels[i] + '.png' )


# Histogram 
kwargs = dict(histtype='step', alpha=0.6, normed=False)
for descriptor in descriptors_labels:
    plt.figure()
    for df in dfs:
        plt.hist(df[descriptor], linewidth=2, fill=False, **kwargs,)
        #plt.legend(csv_names)
        plt.xlabel(descriptor)
        plt.ylabel('Frequency')
        sns.despine()
        plt.savefig(r'graphs\\' + 'hist_' + descriptor + '.png' )

# Calculate the max size of datasets
data_size = 0
for descriptor in descriptors:
    for data in descriptor:
        size = len(data)
        if size > data_size:
            data_size = len(data)
        else: pass
    break
        
index = csv_names
columns = descriptors_labels
descriptor_statistics = pd.DataFrame(index=index, columns=columns)
provisional = []
pose = 0
for descriptor in descriptors:
    for data in descriptor:
        if data_size <= 2000:
            # Shapiro-Wilk normality test
            stat, p = shapiro(data)
        else:
            # Kolmogorov Smirnov - lilliefors
            stat, p = lilliefors(data, dist='norm')
        # interpret
        alpha = 0.05
        if p > alpha:
            provisional.append('Gaussian')
        else:
            provisional.append('No Gaussian')             
    descriptor_statistics[descriptors_labels[pose]] = provisional
    provisional = []  
    pose += 1
    
# determinating if is normal or not
if (descriptor_statistics.RBs == 'No Gaussian').sum() >= (descriptor_statistics.RBs == 'Gaussian').sum():
    normal = False
else:
    normal = True

if normal:
    print("ANOVA")
else:
    for descriptor in descriptors:
        stat, p = kruskal(*descriptor)
        df = pd.DataFrame(descriptor).transpose()
        df.columns = csv_names
        df = df.melt(var_name='groups', value_name='values')
        posthoc = dunnett(df, val_col='values', group_col='groups', p_adjust = 'holm')
        print(p)
        print(posthoc)
    