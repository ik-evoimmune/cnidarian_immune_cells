import pandas as pd

# Ruta del archivo deduplicado
deduplicated_input_path = "deduplicated_output.csv"
filtered_no_C_output_path = "filtered_no_C_output.csv"

# Leer el archivo CSV
deduplicated_table = pd.read_csv(deduplicated_input_path)

# Eliminar filas donde 'GO_ASPECT' es 'C'
filtered_no_C_table = deduplicated_table[deduplicated_table['GO_ASPECT'] != 'C']

# Guardar el archivo resultante
filtered_no_C_table.to_csv(filtered_no_C_output_path, index=False)

