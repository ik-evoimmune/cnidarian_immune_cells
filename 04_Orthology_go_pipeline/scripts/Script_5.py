import pandas as pd

# Ruta del archivo de entrada
goterms_input_path = "Goterms.csv"
deduplicated_goterms_output_path = "deduplicated_Goterms.csv"

# Leer el archivo CSV
goterms_table = pd.read_csv(goterms_input_path)

# Eliminar filas duplicadas
deduplicated_goterms_table = goterms_table.drop_duplicates()

# Guardar el archivo deduplicado
deduplicated_goterms_table.to_csv(deduplicated_goterms_output_path, index=False)

