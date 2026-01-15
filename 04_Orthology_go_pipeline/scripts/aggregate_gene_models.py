import pandas as pd

# Load the file with a comma as the delimiter
file_path = 'Filtered_blastp_results.csv'
df = pd.read_csv(file_path, delimiter=",")

# Group by 'Gene' and aggregate 'Gene_Model' values with semicolon-separated strings
result = df.groupby("Gene")["Gene_Model"].apply(lambda x: "; ".join(sorted(set(x)))).reset_index()

# Save the processed output
output_path = '/path/to/Processed_Gene_GeneModel.csv'
result.to_csv(output_path, index=False)

print(f"Processed file saved to: {output_path}")

