import pandas as pd

# Cargar el archivo CSV
input_file_path = "merged_output_Edited.csv"
filtered_output_path = "filtered_output.csv"

# Leer el archivo CSV
merged_output_edited = pd.read_csv(input_file_path)

# Eliminar las filas donde todas las columnas después de 'Gene' son NaN
filtered_result = merged_output_edited.dropna(subset=merged_output_edited.columns[1:], how='all')

# Guardar el archivo filtrado
filtered_result.to_csv(filtered_output_path, index=False)

