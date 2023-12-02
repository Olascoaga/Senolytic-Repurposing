import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import numpy as np

# Directorio donde se encuentran los archivos CSV
directorio = '.'

# Obtener la lista de archivos CSV en el directorio
archivos_csv = [archivo for archivo in os.listdir(directorio) if archivo.endswith('.csv')]

# Crear una figura con 2 filas y 2 columnas
fig, axs = plt.subplots(2, 2, figsize=(12, 10))
fig.subplots_adjust(hspace=0.3, wspace=0.55)  # Ajustar el espacio entre paneles

# Lista para almacenar los puntos de scatterplot
scatterplots = []

# Iterar a través de los archivos CSV y crear un gráfico para cada uno
for i, archivo in enumerate(archivos_csv):
    # Leer el archivo CSV
    df = pd.read_csv(os.path.join(directorio, archivo))
    
    # Extraer el nombre de la proteína a partir del nombre del archivo
    nombre_proteina = os.path.splitext(archivo)[0]
    
    # Ordenar los datos por "Consensus Ranking"
    df = df.sort_values(by='Consensus Ranking')
    
    # Tomar solo algunas etiquetas de las moléculas (por ejemplo, las 5 primeras)
    num_etiquetas_mostrar = 0
    nombres_moleculas = df['Ligand'].head(num_etiquetas_mostrar)
    
    # Calcular la media de "Affinity Autodock" y "Affinity Vina"
    df['Affinity (kcal/mol)'] = (df['Affinity Autodock'] + df['Affinity Vina']) / 2
    
    # Crear un gráfico de scatterplot en el panel correspondiente con ejes invertidos y gradiente de color
    scatterplot = sns.scatterplot(
        data=df, x='Ranking Vina', y='Ranking Autodock',
        hue='Affinity (kcal/mol)', palette='viridis_r', ax=axs[i//2, i%2
    ])
    
    axs[i//2, i%2].set_title(nombre_proteina)
    axs[i//2, i%2].set_xlabel('Ranking Vina')
    axs[i//2, i%2].set_ylabel('Ranking Autodock')
    
    # Invertir los ejes
    axs[i//2, i%2].invert_xaxis()
    axs[i//2, i%2].invert_yaxis()
    
    for j, nombre_molecula in enumerate(nombres_moleculas):
        axs[i//2, i%2].text(
            df.loc[j, 'Ranking Vina'] + 1,  # Ajusta el desplazamiento en el eje x
            df.loc[j, 'Ranking Autodock'] + 1,  # Ajusta el desplazamiento en el eje y
            nombre_molecula,
            fontsize=8
        )
    
    axs[i//2, i%2].legend(title='Afinidad (kcal/mol)', loc='upper left', bbox_to_anchor=(1, 1))

nombre_archivo = 'scatterplot.png'
plt.savefig(nombre_archivo, dpi=600, bbox_inches='tight')

plt.close()
