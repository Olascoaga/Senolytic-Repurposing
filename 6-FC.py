#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 10:33:08 2023

@author: samael
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os
import glob

# Carga el archivo de genes de interés
with open('scaps_19.txt', 'r') as f:
    genes_interes = [line.strip() for line in f]

# Encuentra todos los archivos TSV en el directorio
archivos_tsv = glob.glob("*.tsv")

# Crea un diccionario para almacenar los datos de cada tejido
datos_tejidos = {}

for archivo in archivos_tsv:
    tejido, rango = os.path.splitext(archivo)[0].split('_')
    if tejido not in datos_tejidos:
        datos_tejidos[tejido] = {}
    datos_tejidos[tejido][rango] = pd.read_csv(archivo, sep='\t')

resultados = []

for tejido, datos in datos_tejidos.items():
    jovenes = datos.get('jovenes')
    viejos = datos.get('viejos')
    
    # En caso de que no haya datos para jóvenes o viejos en un tejido, saltar al siguiente tejido
    if jovenes is None or viejos is None:
        continue

    # Filtra los genes de interés
    jovenes = jovenes[jovenes['Description'].isin(genes_interes)]
    viejos = viejos[viejos['Description'].isin(genes_interes)]
    
    # Aplica log10 más 1 a las columnas numéricas
    jovenes[jovenes.select_dtypes(include=[np.number]).columns] = np.log10(jovenes[jovenes.select_dtypes(include=[np.number]).columns] + 1)
    viejos[viejos.select_dtypes(include=[np.number]).columns] = np.log10(viejos[viejos.select_dtypes(include=[np.number]).columns] + 1)
    
    for gene in genes_interes:
        expresion_jovenes = jovenes[jovenes['Description'] == gene].iloc[:,1:]
        expresion_viejos = viejos[viejos['Description'] == gene].iloc[:,1:]
        
        # Calcula fold change
        fold_change = np.log2(expresion_viejos.mean().mean()) / np.log2(expresion_jovenes.mean().mean())
        
        # Realiza la prueba de Wilcoxon-Mann-Whitney 
        _, p_value = stats.ranksums(expresion_jovenes.values.flatten(), expresion_viejos.values.flatten())
        
        resultados.append({
            'Tejido': tejido,
            'Gene': gene,
            'Fold Change': fold_change,
            'p-value': p_value
        })

# Convierte los resultados en un DataFrame para fácil visualización
df_resultados = pd.DataFrame(resultados)

# Colores para las diferentes condiciones
colores = {'jovenes': 'blue', 'viejos': 'red'}

for tejido, datos in datos_tejidos.items():
    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(polar=True))
    for rango, df in datos.items():
        # Filtra los genes de interés
        df = df[df['Description'].isin(genes_interes)]
        # Calcula la media de las columnas de las muestras y la convierte a log2
        medias = np.log2(df.iloc[:, 1:].mean(axis=1) + 1) # Se añade 0.001 para evitar el logaritmo de 0
        # Preparación de los datos para el gráfico de radar
        angles = np.linspace(0, 2 * np.pi, len(medias), endpoint=False).tolist()
        stats = medias.tolist()
        # Se cierra el gráfico
        stats.append(stats[0])
        angles.append(angles[0])
        # Añade al gráfico
        ax.plot(angles, stats, color=colores[rango], linewidth=1, label=rango)
        ax.scatter(angles, stats, color=colores[rango], alpha=0.75)
    # Añade leyendas y etiquetas
    labels = df['Description'].tolist()
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(labels, fontsize=12, color='black', rotation=0)
    ax.grid(True)
    plt.legend(loc='upper right', bbox_to_anchor=(0.1, 0.1))
    plt.show()
