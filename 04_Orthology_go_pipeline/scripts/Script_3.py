import pandas as pd

# Ruta del archivo filtrado
filtered_input_path = "filtered_output.csv"
deduplicated_output_path = "deduplicated_output.csv"

# Leer el archivo CSV
filtered_table = pd.read_csv(filtered_input_path)

# Eliminar filas duplicadas
deduplicated_table = filtered_table.drop_duplicates()

# Guardar el archivo deduplicado
deduplicated_table.to_csv(deduplicated_output_path, index=False)

