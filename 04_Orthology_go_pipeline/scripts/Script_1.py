import pandas as pd

# Load the uploaded CSV files
correspondence_path = "Correspondence_Genes.csv"
merged_table_path = "merged_table_edited_1.csv"

# Read the CSV files
correspondence_df = pd.read_csv(correspondence_path)
merged_table_df = pd.read_csv(merged_table_path)

# Explode the 'Gene_Model' column in the correspondence table
correspondence_df['Gene_Model'] = correspondence_df['Gene_Model'].str.split('; ')
correspondence_exploded = correspondence_df.explode('Gene_Model')

# Merge the exploded correspondence table with the merged table based on 'Gene_Model'
merged_result = correspondence_exploded.merge(merged_table_df, on='Gene_Model', how='left')

# Save the output
output_path = "merged_output.csv"
merged_result.to_csv(output_path, index=False)

