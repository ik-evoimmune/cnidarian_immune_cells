import pandas as pd

# Load the data file
df = pd.read_csv('Edited_blastp_Arnau.vs.NVE_old.csv')

# Convert evalue column to numeric for comparison
df['evalue'] = pd.to_numeric(df['evalue'], errors='coerce')

# Sort values by Gene, evalue (ascending), and bitscore (descending)
df_sorted = df.sort_values(by=['Gene', 'evalue', 'bitscore'], ascending=[True, True, False])

# Group by Gene and select rows with evalue == 0, or the first row otherwise
filtered_df = df_sorted.groupby('Gene').apply(
    lambda group: group[group['evalue'] == 0] if (group['evalue'] == 0).any() else group.iloc[:1]
).reset_index(drop=True)

# Save the filtered results
filtered_df.to_csv('Filtered_blastp_results.csv', index=False)

