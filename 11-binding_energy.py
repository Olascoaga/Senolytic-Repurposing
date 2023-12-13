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

# Crear un gráfico con subgráficos (paneles), ajusta según el número de proteínas
fig, axs = plt.subplots(2, 2, figsize=(15, 12))  # Ajusta el layout según tus necesidades
axs = axs.flatten()  # Convierte el arreglo 2D en 1D para facilitar el acceso
proteina_idx = 0  # Índice para llevar la cuenta de la proteína actual

# Iterar a través de cada carpeta de proteína
for carpeta_proteina in os.listdir(carpeta_principal):
    ruta_carpeta = os.path.join(carpeta_principal, carpeta_proteina)
    if os.path.isdir(ruta_carpeta):
        datos_binding = pd.DataFrame()  # DataFrame para almacenar datos de esta proteína

        # Iterar a través de cada archivo .tab en la carpeta de la proteína
        for archivo in os.listdir(ruta_carpeta):
            if archivo.endswith('.tab'):
                ruta_completa = os.path.join(ruta_carpeta, archivo)

                # Leer y procesar el contenido del archivo .tab
                with open(ruta_completa, 'r') as file:
                    contenido = file.read()
                    contenido_modificado = contenido.replace('\t', ',').replace(' ', ',')
                df = pd.read_csv(StringIO(contenido_modificado), sep=',+', engine='python', header=None)

                # Ignorar la primera línea y seleccionar las columnas relevantes
                df = df.iloc[1:]
                df = df.iloc[:1000]
                df = df[df.columns[1:3]]

                # Asegurarse de que el DataFrame tiene exactamente 2 columnas
                if df.shape[1] == 2:
                    # Renombrar las columnas para claridad
                    df.columns = ['Time[ns]', 'Binding Energy']
                    df['Ligand'] = archivo.split('.')[0]  # Asignar nombre de ligando
                    datos_binding = pd.concat([datos_binding, df], ignore_index=True)

        # Procesar y graficar datos de esta proteína
        if not datos_binding.empty:
            # Convertir a numérico y ajustar el tiempo
            datos_binding['Time[ns]'] = pd.to_numeric(datos_binding['Time[ns]'], errors='coerce')
            datos_binding['Binding Energy'] = pd.to_numeric(datos_binding['Binding Energy'], errors='coerce')
            datos_binding = datos_binding.dropna(axis=0)
            datos_binding['Time[ns]'] = datos_binding['Time[ns]'] / 1000

            sns.lineplot(
                ax=axs[proteina_idx],
                data=datos_binding,
                x='Time[ns]',
                y='Binding Energy',
                hue='Ligand',
                palette=colores_ligandos
            )
            axs[proteina_idx].set_title(carpeta_proteina)
            axs[proteina_idx].set_xlabel('Tiempo [Nanosegundos]')
            axs[proteina_idx].set_ylabel('Energía de unión [kJ/mol]')
            axs[proteina_idx].legend(title='Ligando')
            proteina_idx += 1

# Ajustar el layout y mostrar el gráfico
plt.tight_layout()
plt.savefig('Binding_Energy.png', dpi=600)
plt.show()
