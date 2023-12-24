#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 17 11:32:08 2023

@author: samael
"""

import pandas as pd
import os
from scipy import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multitest import multipletests

# Obtener la lista de todos los archivos CSV en el directorio actual
archivos = [f for f in os.listdir('.') if f.endswith('.csv')]

# Función para procesar archivos y obtener un DataFrame
def procesar_archivos(archivos, columnas):
    df = pd.DataFrame(columns=columnas)
    for archivo in archivos:
        datos_archivo = pd.read_csv(archivo, header=None)
        datos_archivo = datos_archivo.to_numpy().flatten().tolist()

        datos_temp = {col: [] for col in columnas}
        for col in columnas:
            if datos_archivo:
                datos_temp[col] = datos_archivo[:3]
                datos_archivo = datos_archivo[3:]

        df_temp = pd.DataFrame(datos_temp)
        df = pd.concat([df, df_temp], ignore_index=True)
    return df

# Analizar cada fármaco
resultados = []
for archivo in set(a.split('_')[0] for a in archivos):
    archivos_farmaco_control = [f for f in archivos if f.startswith(f'{archivo}_control')]
    archivos_farmaco_senescente = [f for f in archivos if f.startswith(f'{archivo}_senescente')]

    # Columnas del DataFrame para este fármaco
    columnas = ['Control', 'DMSO', 'D+Q', '10 uM', '20 uM', '30 uM']

    # Procesar archivos de control y senescentes
    df_control = procesar_archivos(archivos_farmaco_control, columnas)
    df_senescentes = procesar_archivos(archivos_farmaco_senescente, columnas)

    # Reestructurar para ANOVA
    df_senescentes_melt = df_senescentes.melt(var_name='group', value_name='value')

    # ANOVA
    fvalue, pvalue = stats.f_oneway(*[df_senescentes_melt[df_senescentes_melt['group'] == group]['value']
                                      for group in df_senescentes_melt['group'].unique()])

    # Resultados de ANOVA
    resultado_anova = {'Fármaco': archivo, 'ANOVA F-valor': fvalue, 'ANOVA P-valor': pvalue}

    # Si el p-valor es menor a 0.05, realizar prueba de Tukey
    if pvalue < 1:
        tukey_results = pairwise_tukeyhsd(df_senescentes_melt['value'], df_senescentes_melt['group'], alpha=0.05)
        tukey_df = pd.DataFrame(data=tukey_results._results_table.data[1:], columns=tukey_results._results_table.data[0])

        # Ajuste de p-valores usando Benjamini-Hochberg
        reject, pvals_corrected, _, _ = multipletests(tukey_df['p-adj'], method='fdr_bh')
        tukey_df['p-adj corrected'] = pvals_corrected
        tukey_df['reject H0'] = reject
        resultado_anova['Tukey'] = tukey_df[['group1', 'group2', 'p-adj', 'p-adj corrected', 'reject H0']]

    resultados.append(resultado_anova)

# Mostrar los resultados
for resultado in resultados:
    print(resultado)
