import os
from PIL import Image
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def calculate_background_intensity(background_image):
    background_gray = background_image.convert('L')
    background_array = np.array(background_gray)
    background_intensity = np.mean(background_array)
    return background_intensity

def process_images(group_name, background_intensity):
    pixel_values = []

    for filename in os.listdir('.'):
        if filename.endswith('.png') and group_name in filename:
            img = Image.open(filename)
            img_gray = img.convert('L')
            img_array = np.array(img_gray)
            img_array = img_array.astype(np.float64) - background_intensity
            img_array[img_array < 0] = 0
            img_array[img_array > 255] = 255
            img_array = img_array.astype(np.uint8)
            if img_array.size > 0:
                pixel_values.append(np.mean(img_array))

    return pixel_values

background_image_control = Image.open('background_control.png')
background_image_senes = Image.open('background_senes.png')

background_intensity_control = calculate_background_intensity(background_image_control)
background_intensity_senes = calculate_background_intensity(background_image_senes)

pixel_values_ctrl = process_images('CTRL', background_intensity_control)
pixel_values_senes = process_images('SENES', background_intensity_senes)

# Realizar una prueba t independiente en lugar de Wilcoxon
t_stat, p_value = stats.ttest_ind(pixel_values_ctrl, pixel_values_senes)

print(f't-statistic: {t_stat}')
print(f'p-value: {p_value}')

sns.set(style="white")

labels = ['Control', 'Control_Space', 'Senescence', 'Senescence_Space']
means = [np.mean(pixel_values_ctrl), 0, np.mean(pixel_values_senes), 0]
errors = [np.std(pixel_values_ctrl) / np.sqrt(len(pixel_values_ctrl)),
          np.nan,
          np.std(pixel_values_senes) / np.sqrt(len(pixel_values_senes)),
          np.nan]

# Crear un DataFrame con los resultados
df = pd.DataFrame({
    'Group': labels,
    'Mean Pixel Intensity': means,
    'SEM': errors
})

# Crear el gráfico de barras con Seaborn
bar_plot = sns.barplot(x='Group', y='Mean Pixel Intensity', data=df, palette=['white', 'white', 'black', 'black'], edgecolor='gray', capsize=0.2)

# Añadir barras de error integradas en las barras
for i, bar in enumerate(bar_plot.patches):
    if i % 2 == 0:  # Solo para los grupos de interés
        height = bar.get_height()
        error = errors[i]
        bar_plot.errorbar(bar.get_x() + bar.get_width() / 2., height, yerr=error, color='black', capsize=5, capthick=2)

bar_plot.set_xticklabels(['Control', '', 'Senescentes', ''])
bar_plot.set_ylabel('Intensidad media', fontsize=14)
bar_plot.set_title('Caspasa-3 activa', fontsize=16)
sns.despine()
plt.show()
