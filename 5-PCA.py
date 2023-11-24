# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 23:17:20 2020

@author: olask
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
 

archivo = open("descriptores.txt", "r")
descriptores = archivo.readlines()
descriptores = list(map(lambda s: s.strip(), descriptores))

df = pd.read_csv('todo.csv', usecols=descriptores)
X = df.iloc[:,0:11].values
y = df.iloc[:,11].values

X_std = StandardScaler().fit_transform(X)
cov_mat = np.cov(X_std.T)
eig_vals, eig_vecs = np.linalg.eig(cov_mat)
    
# A partir de los autovalores, calculamos la varianza explicada
tot = sum(eig_vals)
var_exp = [(i / tot)*100 for i in sorted(eig_vals, reverse=True)]
cum_var_exp = np.cumsum(var_exp)

#  Hacemos una lista de parejas (autovector, autovalor) 
eig_pairs = [(np.abs(eig_vals[i]), eig_vecs[:,i]) for i in range(len(eig_vals))]

# Ordenamos estas parejas den orden descendiente con la función sort
eig_pairs.sort(key=lambda x: x[0], reverse=True)

# Representamos en un diagrama de barras la varianza explicada por cada autovalor, y la acumulada
with plt.style.context('ggplot'):
    plt.figure(figsize=(6, 4))

    plt.bar(range(11), var_exp, alpha=0.6, align='center',
            label='Varianza explicada individual', color='b')
    plt.step(range(11), cum_var_exp, where='mid', linestyle='-', label='Varianza explicada acumulada')
    plt.ylabel('Ratio de Varianza Explicada')
    plt.xlabel('Componentes Principales')
    plt.legend(loc='best')
    plt.tight_layout()


#Generamos la matríz a partir de los pares autovalor-autovector
matrix_w = np.hstack((eig_pairs[0][1].reshape(11,1),
                      eig_pairs[1][1].reshape(11,1)))

Y = X_std.dot(matrix_w)

with plt.style.context('seaborn-white'):
    plt.figure(figsize=(9, 6))
    for lab, col in zip(('Senolytics', 'Exp', 'NP', 'Aproved', 'N015' ),
                        ('purple', 'orange', 'red', 'blue', 'green')):
        plt.scatter(Y[y==lab, 0],
                    Y[y==lab, 1],
                    label=lab,
                    c=col)
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.show()