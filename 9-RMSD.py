#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 04:08:56 2023

@author: samael
"""

import pandas as pd
import os
from io import StringIO
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['font.family'] = 'Arial'

# Definir los colores para cada ligando
colores_ligandos = {
    "Apo": "black",
    "Sotrastaurin": "purple",
    "Bicuculina": "green",
    "Tolvaptan": "blue",
    "Cryptotanshinona": "orange",
    "Inhibidor": "red"
}

# Ruta de la carpeta principal
carpeta_principal = '.'

# Crear un gráfico con 4 subgráficos (paneles)
fig, axs = plt.subplots(2, 2, figsize=(15, 12))  # Ajusta el layout según tus necesidades
axs = axs.flatten()  # Convierte el arreglo 2D en 1D para facilitar el acceso
proteina_idx = 0  # Índice para llevar la cuenta de la proteína actual

# Iterar a través de cada carpeta de proteína
for carpeta_proteina in os.listdir(carpeta_principal):
    ruta_carpeta = os.path.join(carpeta_principal, carpeta_proteina)
    if os.path.isdir(ruta_carpeta):
        datos_proteina = pd.DataFrame()  # DataFrame para almacenar datos de esta proteína

        # Iterar a través de cada archivo .tab en la carpeta de la proteína
        for archivo in os.listdir(ruta_carpeta):
            if archivo.endswith('.tab'):
                ruta_completa = os.path.join(ruta_carpeta, archivo)

                # Leer y procesar el contenido del archivo .tab
                with open(ruta_completa, 'r') as file:
                    contenido = file.read()
                    contenido_modificado = contenido.replace('\t', '').replace(' ', ',')
                df = pd.read_csv(StringIO(contenido_modificado), sep=',+', engine='python', nrows=1000)  # Lee solo los primeros 1000 registros

                # Verificar si las columnas requeridas están presentes
                if 'RMSDCa' in df.columns and 'Time[ns]' in df.columns:
                    df['Ligando'] = archivo.split('.')[0]  # Asignar nombre de ligando
                    datos_proteina = pd.concat([datos_proteina, df], ignore_index=True)

        # Filtrar y graficar datos de esta proteína
        if not datos_proteina.empty:
            sns.lineplot(
                ax=axs[proteina_idx],
                data=datos_proteina,
                x='Time[ns]',
                y='RMSDCa',
                hue='Ligando',
                palette=colores_ligandos
            )
            axs[proteina_idx].set_title(carpeta_proteina)
            axs[proteina_idx].set_xlabel('Tiempo [Nanosegundos]')
            axs[proteina_idx].set_ylabel('RMSD (Å)')
            axs[proteina_idx].legend(title='Ligando')
            proteina_idx += 1

# Mostrar el gráfico
plt.tight_layout()
plt.savefig('RMSD.png', dpi=600)
plt.show()
