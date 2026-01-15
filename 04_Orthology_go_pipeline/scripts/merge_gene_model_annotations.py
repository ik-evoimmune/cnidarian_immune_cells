import pandas as pd

# Load the two files
file1_path = 'Arnau.csv'
file2_path = 'Processed_Gene_GeneModel.csv'

# Read the files
df1 = pd.read_csv(file1_path, delimiter=",", header=None, names=["Gene", "Gene_Model"])
df2 = pd.read_csv(file2_path, delimiter=",", header=0)

# Merge the two dataframes on the "Gene" column
merged = pd.merge(df1, df2, on="Gene", how="outer", suffixes=("_file1", "_file2"))

# Function to combine and deduplicate entries while cleaning up spaces and semicolons
def combine_models(row):
    models_file1 = str(row["Gene_Model_file1"]).split(";") if pd.notna(row["Gene_Model_file1"]) else []
    models_file2 = str(row["Gene_Model_file2"]).split(";") if pd.notna(row["Gene_Model_file2"]) else []
    # Remove duplicates, extra spaces, and trailing semicolons
    combined = sorted(set(model.strip().rstrip(";") for model in models_file1 + models_file2))
    return "; ".join(combined)

merged["Combined_Gene_Model"] = merged.apply(combine_models, axis=1)

# Select relevant columns for the output
result = merged[["Gene", "Combined_Gene_Model"]]

# Save the final output
output_path = 'Merged_Gene_GeneModel_Corrected.csv'
result.to_csv(output_path, index=False)

print(f"Merged and cleaned file saved to: {output_path}")

