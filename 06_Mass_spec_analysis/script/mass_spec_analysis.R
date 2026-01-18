# Load libraries
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)

mass_spec<- readxl::read_excel("~/immune_cells/cnidarian_immune_cells/06_Mass_spec_analysis/input_files/MaxQuant-Perseus-stats-Rubi-20210201.xlsx", sheet = 2)

glimpse(mass_spec)

#Read dictionary file
dictionary<- read.csv("~/immune_cells/cnidarian_immune_cells/06_Mass_spec_analysis/input_files/Dictionary_Kozlovski_et.al.csv")

# Assign gene names to NVE gene models 

library(stringr)

# 1) Build a mapping: NVE gene -> protein_homolog(s)
gene2homolog <- dictionary %>%
  transmute(
    NVE_gene_model = as.character(NVE_gene_model),
    protein_homolog = na_if(trimws(as.character(protein_homolog)), "")
  ) %>%
  separate_rows(NVE_gene_model, sep = ";") %>%          # handle multi-IDs
  filter(!is.na(NVE_gene_model)) %>%
  group_by(NVE_gene_model) %>%
  summarise(
    protein_homolog_match = paste(unique(na.omit(protein_homolog)), collapse = ";"),
    .groups = "drop"
  )


# 2) Add homolog info to mass_spec WITHOUT duplicating final rows
mass_spec_with_homolog <- mass_spec %>%
  mutate(
    row_id = row_number(),
    `Fasta headers` = as.character(`Fasta headers`)
  ) %>%
  separate_rows(`Fasta headers`, sep = ";") %>%         # split multi-IDs only for lookup
  left_join(gene2homolog, by = c("Fasta headers" = "NVE_gene_model")) %>%
  group_by(row_id) %>%
  summarise(
    across(-c(`Fasta headers`, protein_homolog_match), first),
    protein_homolog_match = paste(unique(na.omit(protein_homolog_match)), collapse = ";"),
    .groups = "drop"
  ) %>%
  select(-row_id)

head(mass_spec_with_homolog)

# Rank by p-value
df <- mass_spec_with_homolog %>%
  arrange(desc(`log Difference`)) %>%
  mutate(rank = row_number())

head(df)

# Note: #1 RLRb. Ambiguous: #3 ALDH1A2;ALDH2, #5  ADGRL3;FBN2
# For simplicity:

df$protein_homolog_match [1] <- "RLRb"
df$protein_homolog_match [2] <- "TRIP11"
df$protein_homolog_match [3] <- "ALDH1A2"
df$protein_homolog_match [4] <- "LRP1B"
df$protein_homolog_match [5] <- "ADGRL3"

top_hits <- df %>% filter(rank <= 5)

print(top_hits)

library(ggplot2)
library(ggrepel)

ggplot(df, aes(x = rank, y = `log Difference`)) +
  geom_point(size = 1, alpha = 0.7, color = "black") +
  geom_text_repel(
    data = top_hits,
    aes(label = `protein_homolog_match`),
    size = 3.5,
    fontface = "italic",
    color = "red",
    max.overlaps = 100,
    segment.size = 0.2
  ) +
  labs(
    title = "Top Enriched Proteins in α-RLRb IP vs IgG Control",
    x = "Ranked Proteins (by log2 LFQ difference)",
    y = expression(Log[2]*"(RLRb / IgG) LFQ Intensity")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  )


# How many orders of magnitude? 

ratio <- top_hits$`ratio of LFQ intensity: RLRb vs IgG`[1]
log10(ratio)

# p value for RLRb

print(top_hits$`p-value`[1])

sessionInfo()

