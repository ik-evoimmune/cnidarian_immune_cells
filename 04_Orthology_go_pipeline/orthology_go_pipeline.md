---
title: "Orthology and Gene Model Correspondence Pipeline"
author: "Adrian Jaimes-Becerra"
date: "2026-01-07"
---

# Orthology and Gene Model Correspondence Pipeline

## Author

**Adrian Jaimes-Becerra** — Development of this pipeline.

This document describes the workflow used to infer orthologous gene groups, construct gene–gene model correspondence tables, and run Gene Ontology (GO) annotations across different gene models for downstream comparative analyses.

## **OrthoFinder**

How to run:

``` bash
sbatch run_orthofinder.slurm.sh
```

``` bash

#!/bin/bash
#SBATCH -c 12
#SBATCH --mem-per-cpu 2000
#SBATCH --time=96:00:00
#SBATCH -J orthofinder
#SBATCH -o orthofinder.%A.out
#SBATCH -e orthofinder.%A.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<ajjb80@gmail.com>

eval "$(/sci/labs/yehum79/adrianjjb/miniconda3/bin/conda shell.zsh hook)"
conda activate orthofinder

orthofinder -f  /sci/labs/yehum79/adrianjjb/Itamar_Orthology/
```

This SLURM batch script runs **OrthoFinder** to infer orthologous gene groups across multiple species using protein FASTA files as input. The analysis is executed in a dedicated Conda environment on an HPC cluster, allocating 12 CPUs and 96 hours of runtime. OrthoFinder performs performs all-versus-all sequence similarity searches, orthogroup clustering, and gene tree inference for downstream comparative genomic analyseperforms all-versus-all sequence similarity searches, orthogroup clustering, and gene tree inference for downstream comparative genomic analysess.

## Correspondence Between Gene Models

A set of custom Python scripts was used to define high-confidence correspondences between genes and alternative gene models based on BLASTp results.

How to run:

``` bash

python filter_blastp_best_hits.py
```

``` python

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
```

This script processes BLASTp results by ranking hits for each gene based on e-value (ascending) and bitscore (descending). For each gene, it preferentially retains hits with an e-value of zero; if none are present, the best-ranked hit is selected. The output is a filtered BLASTp table used for defining high-confidence gene–gene model correspondences.

``` bash

python aggregate_gene_models.py
```

``` python

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

```

This script groups filtered BLASTp results by gene and aggregates all associated gene models into a single semicolon-separated list. Duplicate gene models are removed and entries are sorted to generate a clean gene → gene model correspondence table, which is later used for integration with orthology results.

``` bash

python merge_gene_model_annotations.py
```

``` python

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

```

This script merges two independent gene–gene model annotation tables using an outer join on gene identifiers. Gene models from both sources are combined, deduplicated, and cleaned of formatting artifacts (extra spaces and trailing delimiters). The resulting unified table ensures maximal recovery of gene–gene model correspondences for downstream comparative analyses.

A series of custom Python scripts were developed to construct correspondence tables between Gene Ontology (GO) terms and alternative gene models, enabling consistent functional annotation across datasets.

**Script 1 – Gene model expansion and GO table merging**

**Description:**

Expands multi-gene model entries and merges the gene correspondence table with GO annotation data to generate a unified gene–GO mapping table.

``` python

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

```

**Script 2 – Removal of empty GO annotations**

**Description:**

Filters the merged table by removing rows lacking GO annotations, retaining only genes with at least one associated GO term.

``` python

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
```

**Script 3 – Deduplication of gene–GO entries**

**Description:**

Removes duplicate gene–GO associations to generate a non-redundant correspondence table.

``` python

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
```

**Script 4 – GO ontology filtering**

**Description:**

Filters GO annotations by ontology category, excluding Cellular Component (C) terms and retaining Biological Process and Molecular Function annotations.

``` python

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
```

**Script 5 – Deduplication of GO term reference table**

**Description:**

Removes duplicate entries from the GO term reference table to produce a curated, non-redundant GO term dataset.

``` python

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
```
