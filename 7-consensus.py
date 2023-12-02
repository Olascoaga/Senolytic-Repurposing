#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 22:46:39 2023

@author: samael
"""

import os
import pandas as pd

# Ruta de la carpeta donde se encuentran los archivos CSV
carpeta = "."

# Obtener la lista de archivos CSV en la carpeta
archivos_csv = [archivo for archivo in os.listdir(carpeta) if archivo.endswith(".csv")]

# Crear un diccionario para almacenar los datos de Z score de cada programa y el ranking original
datos_programas = {}

# Obtener la lista completa de moléculas para el consenso
moléculas_consenso = set()

# Iterar a través de los archivos CSV y procesarlos
for archivo_csv in archivos_csv:
    # Construir la ruta completa del archivo
    ruta_completa = os.path.join(carpeta, archivo_csv)
    
    # Leer el archivo CSV en un DataFrame
    df = pd.read_csv(ruta_completa)
    
    # Almacenar las moléculas en la lista completa de moléculas para el consenso
    moléculas_consenso.update(df['Ligand'])
    
    # Almacenar los datos de Z score en el diccionario usando el nombre del archivo como clave
    datos_programas[archivo_csv] = df

# Crear un DataFrame de consenso con todas las moléculas
consenso = pd.DataFrame({'Ligand': list(moléculas_consenso)})

# Calcular el Z score promedio para el consenso como el promedio de los Z scores de los dos programas
for programa, datos in datos_programas.items():
    datos = datos[['Ligand', 'Z score']].rename(columns={'Z score': programa})
    consenso = pd.merge(consenso, datos, on='Ligand', how='left')

# Calcular el Z score promedio
consenso['Z score promedio'] = consenso.iloc[:, 1:].mean(axis=1)

# Generar el ranking de consenso
ranking_consenso = consenso[['Ligand', 'Z score promedio']].sort_values(by='Z score promedio', ascending=False)
ranking_consenso['Ranking Consenso'] = range(1, len(ranking_consenso) + 1)

# Agregar la posición en el ranking de consenso al DataFrame de consenso
consenso = pd.merge(consenso, ranking_consenso[['Ligand', 'Ranking Consenso']], on='Ligand', how='left')

# Generar los rankings de cada programa y las afinidades
for programa, datos in datos_programas.items():
    ranking_programa = datos[['Ligand', 'Affinity', 'Z score']].sort_values(by='Z score', ascending=False)
    ranking_programa['Ranking ' + programa] = range(1, len(ranking_programa) + 1)
    ranking_programa = ranking_programa.rename(columns={'Affinity': 'Affinity ' + programa, 'Z score': 'Z score ' + programa})
    consenso = pd.merge(consenso, ranking_programa[['Ligand', 'Ranking ' + programa, 'Affinity ' + programa]], on='Ligand', how='left')

consenso.to_csv('consenso.csv', index=False)