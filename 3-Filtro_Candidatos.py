#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 18:40:50 2023

@author: samael
"""

import glob
import os
import pandas as pd
import numpy as np

# Obtén la lista de archivos CSV en la carpeta actual
archivos_csv = glob.glob('*.csv')

# Define el nombre del archivo quimioteca.csv y el nombre de la columna 'ID'
archivo_quimioteca = 'lib/quimioteca.csv'
archivo_farmacos = 'lib/experimental.csv'

# Cargar el DataFrame de quimioteca.csv
df_quimioteca = pd.read_csv(archivo_quimioteca)
df_farmacos = pd.read_csv(archivo_farmacos)
nombres = df_quimioteca['ID'].tolist()
farmacos = df_farmacos['ID'].tolist()

matrices = []
# Itera sobre los archivos CSV encontrados en la carpeta actual
for archivo_csv in archivos_csv:
    matriz = pd.read_csv(archivo_csv, index_col=0)

    matriz.columns = nombres 
    matriz.index = nombres
    matriz = matriz.loc[farmacos, farmacos]

    matriz = matriz.loc[~matriz.index.duplicated(keep='first')]
    matriz = matriz.loc[:, ~matriz.columns.duplicated(keep='first')]
    matrices.append(matriz)

matriz_promedio = np.mean(matrices, axis=0)
# Crear un DataFrame a partir de matriz_promedio
matriz_promedio = pd.DataFrame(matriz_promedio, index=matrices[0].index, columns=matrices[0].columns)

matriz_promedio = matriz_promedio.T.round(3)
matriz_promedio[matriz_promedio < 0.7] = 0

# Crear un diccionario para mapear los nombres originales a "Compound 1, 2, 3"
diccionario_nombres = {nombre: f"Compound {i+1}" for i, nombre in enumerate(matriz_promedio.index)}

# Reemplazar los nombres en la matriz_promedio
matriz_promedio = matriz_promedio.rename(index=diccionario_nombres, columns=diccionario_nombres)

# Guardar el diccionario en un archivo CSV
diccionario_df = pd.DataFrame(list(diccionario_nombres.items()), columns=['Nombre Original', 'Nombre Nuevo'])
diccionario_df.to_csv('results/diccionario_nombres.csv', index=False)

matriz_promedio = matriz_promedio.iloc[:70].T
"""
matriz_promedio = matriz_promedio.stack()
matriz_promedio = matriz_promedio[(matriz_promedio != 1) & (matriz_promedio != 0)]
# Convertir la Serie en un DataFrame
matriz_promedio = matriz_promedio.reset_index()

# Asignar nombres a las tres primeras columnas
matriz_promedio.columns = ['Nombre1', 'Nombre2', 'Nombre3'] + list(matriz_promedio.columns[3:])
"""
matriz_promedio.to_csv('results/consensus_matrix.csv')
