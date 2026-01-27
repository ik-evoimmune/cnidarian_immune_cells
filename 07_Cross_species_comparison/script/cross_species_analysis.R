
# Reading the files and data wrangling ------------------------------------


library(tidyverse)

# Read orthogroups file
orthogroups<- read_tsv("~/immune_cells/cross_species/orthofinder/Orthogroups.tsv")
# Read cluster 1 markers 
nvec_de<- read.csv("~/immune_cells/cross_species/mast_results.csv")
styl_de<- read.csv("~/immune_cells/cross_species/DE_Spi/Spi_DE_2.3.csv")


# 1. Expand the orthogroups table: one row per gene per species
orthogroups_long <- orthogroups %>%
  separate_rows(Nvec, sep = ",\\s*") %>%
  separate_rows(Spi, sep = ",\\s*") %>%
  filter(Nvec != "", Spi != "") %>%
  mutate(
    Nvec = trimws(Nvec),
    Spi = trimws(Spi)
  )

# 2. Modify gene IDs format in the Nematostella differential expression table
nvec_de <- nvec_de %>%
  mutate(gene_id = gsub("-", "_", X))

# 3. Modify gene IDs format in the Stylophora differential expression table
styl_de <- styl_de %>%
  mutate(gene_id = gsub("^Spis_", "", gene),         # remove "Spis_"
         gene_id = sub("_(\\d+)$", ".\\1", gene_id)) # convert "_1" → ".1"

# 4. Join Nematostella differential expression results to orthogroups
orthogroups_de <- orthogroups_long %>%
  left_join(nvec_de, by = c("Nvec" = "gene_id"))

# 5. Join Stylophora differential expression results to orthogroups
orthogroups_de <- orthogroups_long %>%
  left_join(nvec_de, by = c("Nvec" = "gene_id")) %>%
  left_join(styl_de, by = c("Spi" = "gene_id"), suffix = c("_nvec", "_styl"))

# 6. Check how many matches succeeded
cat("Matched Nvec DE genes:\n")
print(table(is.na(orthogroups_de$avg_log2FC)))

cat("Matched Stylophora DE genes:\n")
print(table(is.na(orthogroups_de$log2FoldChange)))

# 7. Filter for for non NAs genes in both species
orthogroups_common_de <- orthogroups_de %>%
  filter(!is.na(avg_log2FC) & !is.na(log2FoldChange))

# 8. Scatter plot for comparison: all orthogroups from each group
library(ggplot2)

ggplot(orthogroups_common_de, aes(x = avg_log2FC, y = log2FoldChange)) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  labs(
    x = "log2FC (Nematostella)",
    y = "log2FC (Stylophora)",
    title = "Conserved DE genes across orthologs"
  )

# 9. Run a correlation test when p value in each dataset is < 0.05 (differentially expressed)

orthogroups_sig <- orthogroups_common_de %>%
  filter(p_val_adj < 0.05, padj < 0.05)

cor_result <- cor.test(
  orthogroups_sig$avg_log2FC,
  orthogroups_sig$log2FoldChange,
  method = "pearson"  # or "spearman"
)
print(cor_result)


# 10. Visualize the result


ggplot(orthogroups_sig, aes(x = avg_log2FC, y = log2FoldChange)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  theme_minimal() +
  labs(
    title = "Correlation of DE log2FC across orthologs",
    subtitle = paste0("Pearson r = ", round(cor_result$estimate, 2),
                      ", p = ", signif(cor_result$p.value, 2)),
    x = "log2FC (Nematostella)",
    y = "log2FC (Stylophora)"
  )

cor.test(orthogroups_sig$avg_log2FC, orthogroups_sig$log2FoldChange, method = "spearman")


# Correlation between Stylophora and Nematostella expressed orthol --------

# Visualize the disribution of orthogroups per gene
# How many orthologs (Stylophora proteins) map to each Nvec gene?

orthologs_per_gene <- orthogroups_de %>%
  filter(!is.na(Spi)) %>%
  count(Nvec, name = "n_orthologs")

ggplot(orthologs_per_gene, aes(x = pmin(n_orthologs, 20))) +
  geom_histogram(binwidth = 1, fill = "black", color = "white") +
  scale_x_continuous(
    breaks = 1:20,
    labels = c(1:19, "20+")
  ) +
  labs(
    x = "Number of Stylophora orthologs per Nvec gene",
    y = "Number of Nvec genes",
    title = "Distribution of ortholog counts per gene"
  ) + 
  geom_vline(
    xintercept = median(orthologs_per_gene$n_orthologs),
    linetype = 2
  ) +
  theme_classic() 

# Number of orthologs excedding 20 in Stylophora

sum(orthologs_per_gene$n_orthologs > 20)
max(orthologs_per_gene$n_orthologs)

# This can be interesting for future analysis.There is some massive expansion.  
# The vast majority of genes (~11-12k) have one ortholog

# Correcting for redundancy in orthofinder - taking the median
orthogroup_summary <- orthogroups_de %>%
  group_by(Orthogroup) %>%
  summarise(
    median_log2FC_nvec = median(avg_log2FC, na.rm = TRUE),
    median_log2FC_styl = median(log2FoldChange, na.rm = TRUE),
    n_nvec_genes = n_distinct(Nvec),
    n_styl_genes = n_distinct(Spi)
  ) %>%
  filter(!is.na(median_log2FC_nvec) & !is.na(median_log2FC_styl))

ggplot(orthogroup_summary, aes(x = median_log2FC_nvec, y = median_log2FC_styl)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  theme_minimal() +
  labs(
    title = "Orthogroup-level log2FC correlation",
    x = "Mean log2FC (Nematostella)",
    y = "Mean log2FC (Stylophora)"
  )


# 1. Aggregate by orthogroup (mean expression)
orthogroup_summary <- orthogroups_de %>%
  group_by(Orthogroup) %>%
  summarise(
    log2FC_Nvec = median(avg_log2FC, na.rm = TRUE),
    log2FC_Spi = median(log2FoldChange, na.rm = TRUE)
  ) %>%
  filter(!is.na(log2FC_Nvec) & !is.na(log2FC_Spi)) %>%
  mutate(
    category = case_when(
      log2FC_Nvec * log2FC_Spi > 1 ~ "conserved",
      log2FC_Nvec * log2FC_Spi < -1 ~ "divergent",
      TRUE ~ "neutral"
    )
  )

# Compute Pearson correlation
cor_res <- cor.test(
  orthogroup_summary$log2FC_Nvec,
  orthogroup_summary$log2FC_Spi,
  method = "pearson"
)

print(cor_res)
# Extract R and p-value
r_val <- round(cor_res$estimate, 2)
p_val <- signif(cor_res$p.value, 2)

# Generate the plot
ggplot(orthogroup_summary, aes(x = log2FC_Nvec, y = log2FC_Spi)) +
  stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", alpha = 0.3, contour_var = "density", bins = 12) +
  geom_point(aes(color = category), alpha = 0.6, size = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  scale_color_manual(values = c(conserved = "forestgreen", divergent = "firebrick", neutral = "gray40")) +
  scale_fill_viridis_c(option = "D", begin = 0.2, end = 0.9) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Conserved vs Divergent DE Genes Across Orthogroups",
    subtitle = paste0("Pearson r = ", r_val, ", p = ", p_val),
    x = "log2FC (Nematostella)",
    y = "log2FC (Stylophora)",
    color = "Expression Pattern",
    fill = "Density"
  )

# Improved version

# Density
library(ggplot2)
library(ggpubr)
library(ggpointdensity)
library(dplyr)

# Compute correlation (for annotation)
cor_res <- cor.test(orthogroup_summary$log2FC_Nvec, orthogroup_summary$log2FC_Spi)
r_val <- round(cor_res$estimate, 2)
p_val <- signif(cor_res$p.value, 2)

# Plot with ggpointdensity
ggplot(orthogroup_summary, aes(x = log2FC_Nvec, y = log2FC_Spi)) +
  geom_pointdensity(adjust = 1/3, size = 1.2) +  # density by point
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  stat_cor(
    method = "pearson",
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = min(orthogroup_summary$log2FC_Nvec, na.rm = TRUE),
    label.y = max(orthogroup_summary$log2FC_Spi, na.rm = TRUE),
    size = 5
  ) +
  scale_color_viridis_c(option = "D") +
  theme_pubr(base_size = 14) +
  labs(
    title = "Correlation of log2FC Across Orthogroups",
    subtitle = paste0("Pearson r = ", r_val, ", p = ", p_val),
    x = "log2FC (Nematostella)",
    y = "log2FC (Stylophora)",
    color = "Local density"
  )


# 1:1 orthogroups 
orthogroups_1to1 <- orthogroups_de %>%
  group_by(Orthogroup) %>%
  filter(n_distinct(Nvec) == 1 & n_distinct(Spi) == 1) %>%
  ungroup()

orthogroup_1to1_fc <- orthogroups_1to1 %>%
  select(Orthogroup, Nvec, Spi, avg_log2FC, log2FoldChange, p_val_adj, padj) %>%
  filter(!is.na(avg_log2FC) & !is.na(log2FoldChange)) %>%
  mutate(
    category = case_when(
      avg_log2FC * log2FoldChange > 1 ~ "conserved",
      avg_log2FC * log2FoldChange < -1 ~ "divergent",
      TRUE ~ "neutral"
    )
  )

library(ggpointdensity)
library(ggpubr)
library(ggplot2)

# Correlation
cor_res <- cor.test(orthogroup_1to1_fc$avg_log2FC, orthogroup_1to1_fc$log2FoldChange)
r_val <- round(cor_res$estimate, 2)
p_val <- signif(cor_res$p.value, 2)

# Plot
ggplot(orthogroup_1to1_fc, aes(x = avg_log2FC, y = log2FoldChange)) +
  geom_pointdensity(size = 1.2, adjust = 1/3) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  stat_cor(
    method = "pearson",
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = min(orthogroup_1to1_fc$avg_log2FC),
    label.y = max(orthogroup_1to1_fc$log2FoldChange),
    size = 5
  ) +
  scale_color_viridis_c(option = "D") +
  theme_pubr(base_size = 14) +
  labs(
    title = "Correlation of log2FC in 1:1 Orthologs",
    subtitle = paste0("Pearson r = ", r_val, ", p = ", p_val),
    x = "log2FC (Nematostella)",
    y = "log2FC (Stylophora)",
    color = "Local density"
  )

orthogroup_1to1_fc <- orthogroup_1to1_fc %>%
  filter(p_val_adj < 0.05 & padj < 0.05)


# Take the max instead of mean to visualizie inducible genes - immune response genes
orthogroup_summary2 <- orthogroups_de %>%
  group_by(Orthogroup) %>%
  summarise(
    max_log2FC_nvec = max(avg_log2FC, na.rm = TRUE),
    max_log2FC_styl = max(log2FoldChange, na.rm = TRUE),
    n_nvec_genes = n_distinct(Nvec),
    n_styl_genes = n_distinct(Spi)
  ) %>%
  filter(!is.na(max_log2FC_nvec) & !is.na(max_log2FC_styl))

orthogroup_summary2<- orthogroup_summary2 %>% filter()

ggplot(orthogroup_summary2, aes(x = max_log2FC_nvec, y = max_log2FC_styl)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  theme_minimal() +
  labs(
    title = "Orthogroup-level log2FC correlation",
    x = "Mean log2FC (Nematostella)",
    y = "Mean log2FC (Stylophora)"
  )

# Same as previos analysis (in the paper)
orthogroups_sig <- orthogroups_de %>%
  filter(p_val_adj < 0.05, padj < 0.05) %>%
  filter(!is.na(avg_log2FC), !is.na(log2FoldChange))

orthogroups_sig_max <- orthogroups_sig %>%
  mutate(abs_sum = abs(avg_log2FC) + abs(log2FoldChange)) %>%
  group_by(Nvec) %>%
  slice_max(abs_sum, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    category = case_when(
      avg_log2FC * log2FoldChange > 1 ~ "conserved",
      avg_log2FC * log2FoldChange < -1 ~ "divergent",
      TRUE ~ "neutral"
    )
  )

library(ggpubr)
library(ggpointdensity)
library(ggplot2)

# Compute correlation
cor_res <- cor.test(
  orthogroups_sig_max$avg_log2FC,
  orthogroups_sig_max$log2FoldChange,
  method = "pearson"
)

# Plot
ggplot(orthogroups_sig_max, aes(x = avg_log2FC, y = log2FoldChange)) +
  geom_pointdensity(size = 1.2, adjust = 1/3) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  stat_cor(
    method = "pearson",
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = min(orthogroups_sig_max$avg_log2FC),
    label.y = max(orthogroups_sig_max$log2FoldChange),
    size = 5
  ) +
  scale_color_viridis_c(option = "D") +
  theme_pubr(base_size = 14) +
  labs(
    title = "Correlation of Significant, Maximally Regulated Genes",
    x = "log2FC (Nematostella)",
    y = "log2FC (Stylophora)",
    color = "Local density"
  )

# Polished
ggplot(orthogroups_sig_max, aes(x = avg_log2FC, y = log2FoldChange)) +
  geom_pointdensity(size = 1.2, adjust = 1/3) +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "gray60", alpha = 0.2) +
  scale_color_viridis_c(option = "plasma", name = "Point density") +
  coord_fixed() +
  scale_x_continuous(limits = c(-6, 8), breaks = seq(-6, 8, 2)) +
  scale_y_continuous(limits = c(-5, 6), breaks = seq(-5, 5, 2)) +
  labs(
    title = "Correlation of Significant, Maximally Regulated Genes",
    x = "log2FC (Nematostella)",
    y = "log2FC (Stylophora)"
  ) +
  annotate(
    "text", x = -5.5, y = 5.5, size = 5,
    label = bquote(R^2 == .(round(cor_res$estimate^2, 2)) ~ "," ~ p < .(format.pval(cor_res$p.value, digits = 2, eps = 1e-16)))
  ) +
  theme_pubr(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

# Unfiltered 
# Filter out NA and keep strongest representative (no p-value filtering)
orthogroups_full_max <- orthogroups_de %>%
  filter(!is.na(avg_log2FC), !is.na(log2FoldChange)) %>%
  mutate(abs_sum = abs(avg_log2FC) + abs(log2FoldChange)) %>%
  group_by(Nvec) %>%
  slice_max(abs_sum, n = 1) %>%
  ungroup()

library(ggrepel)

top_labels <- orthogroups_full_max %>%
  filter(protein != "") %>%
  arrange(desc(abs(avg_log2FC * log2FoldChange))) %>%
  slice_head(n = 15)  # Adjust number of labels as needed

library(ggplot2)
library(ggpointdensity)
library(ggpubr)

# Correlation
cor_res <- cor.test(
  orthogroups_full_max$avg_log2FC,
  orthogroups_full_max$log2FoldChange,
  method = "pearson"
)

# Final polished plot
ggplot(orthogroups_full_max, aes(x = avg_log2FC, y = log2FoldChange)) +
  geom_pointdensity(size = 1.2, adjust = 1/3) +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "gray60", alpha = 0.2) +
  scale_color_viridis_c(option = "plasma", name = "Point density") +
  geom_text_repel(
    data = top_labels,
    aes(label = protein),
    size = 3.5,
    max.overlaps = 20,
    box.padding = 0.4,
    seed = 123
  ) +
  coord_fixed() +
  scale_x_continuous(limits = c(-6, 8), breaks = seq(-6, 8, 2)) +
  scale_y_continuous(limits = c(-5, 6), breaks = seq(-5, 5, 2)) +
  labs(
    title = "Correlation of log2FC Across All Orthologs",
    subtitle = bquote(R^2 == .(round(cor_res$estimate^2, 2)) ~ "," ~ p < .(format.pval(cor_res$p.value, digits = 2, eps = 1e-16))),
    x = "log2FC (Nematostella)",
    y = "log2FC (Stylophora)"
  ) +
  theme_pubr(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

ggsave("~/immune_cells/cnidarian_immune_cells/07_Cross_species_comparison/results/ortholog_correlation_full_labeled.png", width = 8, height = 6)



# Venn diagrams -----------------------------------------------------------

# Clean: account for NAs, duplicates etc...
# Load required packages
library(dplyr)
library(ggvenn)
library(gridExtra)
library(ggplot2)
library(svglite)

#1: Filter orthogroup table – remove NAs, keep unique gene–orthogroup mappings
orthogroups_filtered <- orthogroups_de %>%
  filter(!is.na(Nvec), !is.na(Spi)) %>%
  distinct(Nvec, Spi, Orthogroup)

#2: Keep only DE genes that are found in the orthogroup mapping
nvec_has_orthogroup <- unique(orthogroups_filtered$Nvec)
styl_has_orthogroup <- unique(orthogroups_filtered$Spi)

nvec_up_final <- intersect(nvec_up, nvec_has_orthogroup)
nvec_down_final <- intersect(nvec_down, nvec_has_orthogroup)
styl_up_final <- intersect(styl_up, styl_has_orthogroup)
styl_down_final <- intersect(styl_down, styl_has_orthogroup)

#3: Map DE genes to orthogroups
nvec_up_orthos <- orthogroups_filtered %>%
  filter(Nvec %in% nvec_up_final) %>%
  pull(Orthogroup) %>%
  unique()

styl_up_orthos <- orthogroups_filtered %>%
  filter(Spi %in% styl_up_final) %>%
  pull(Orthogroup) %>%
  unique()

nvec_down_orthos <- orthogroups_filtered %>%
  filter(Nvec %in% nvec_down_final) %>%
  pull(Orthogroup) %>%
  unique()

styl_down_orthos <- orthogroups_filtered %>%
  filter(Spi %in% styl_down_final) %>%
  pull(Orthogroup) %>%
  unique()

# Step 4: Create Venn diagrams
venn_up <- list(
  Nematostella_Up = nvec_up_orthos,
  Stylophora_Up = styl_up_orthos
)

venn_down <- list(
  Nematostella_Down = nvec_down_orthos,
  Stylophora_Down = styl_down_orthos
)

p1 <- ggvenn(venn_up,
             fill_color = c("steelblue", "darkorange"),
             stroke_size = 0.5,
             set_name_size = 5,
             text_size = 5) +
  ggtitle("Upregulated Orthogroups")

p2 <- ggvenn(venn_down,
             fill_color = c("skyblue", "tomato"),
             stroke_size = 0.5,
             set_name_size = 5,
             text_size = 5) +
  ggtitle("Downregulated Orthogroups")

# Step 5: Combine vertically and export
combined_plot <- grid.arrange(p1, p2, ncol = 1)

# Save to png
ggsave("~/immune_cells/cnidarian_immune_cells/07_Cross_species_comparison/results/orthogroup_venn_diagrams.png", combined_plot, width = 6, height = 10)

# Save to SVG
ggsave("~/immune_cells/figures_revised/orthogroup_venn_diagrams.svg", combined_plot, width = 6, height = 10, device = "svg")


ggsave("~/immune_cells/figures_revised/venn_up_down_shared.svg", combined_plot, width = 6, height = 10)


# Extract results

# Filter orthogroup table for only non-NA and unique mappings
orthos <- orthogroups_de %>%
  filter(!is.na(Nvec), !is.na(Spi)) %>%
  distinct(Nvec, Spi, Orthogroup)

# Upregulated
nvec_up_orthos <- orthos %>% filter(Nvec %in% nvec_up) %>% pull(Orthogroup) %>% unique()
styl_up_orthos <- orthos %>% filter(Spi %in% styl_up) %>% pull(Orthogroup) %>% unique()

# Downregulated
nvec_down_orthos <- orthos %>% filter(Nvec %in% nvec_down) %>% pull(Orthogroup) %>% unique()
styl_down_orthos <- orthos %>% filter(Spi %in% styl_down) %>% pull(Orthogroup) %>% unique()

# Print sizes
cat("UP:", "\n")
cat("  Nvec up orthogroups:", length(nvec_up_orthos), "\n")
cat("  Styl up orthogroups:", length(styl_up_orthos), "\n")
cat("  Shared up:", length(intersect(nvec_up_orthos, styl_up_orthos)), "\n\n")

cat("DOWN:", "\n")
cat("  Nvec down orthogroups:", length(nvec_down_orthos), "\n")
cat("  Styl down orthogroups:", length(styl_down_orthos), "\n")
cat("  Shared down:", length(intersect(nvec_down_orthos, styl_down_orthos)), "\n")

# Extract genes of interest
shared_up_orthos <- intersect(nvec_up_orthos, styl_up_orthos)
shared_down_orthos <- intersect(nvec_down_orthos, styl_down_orthos)

writeLines(shared_up_orthos, "~/immune_cells/cnidarian_immune_cells/07_Cross_species_comparison/results/shared_upregulated_orthogroups_new.txt")
writeLines(shared_down_orthos, "~/immune_cells/cnidarian_immune_cells/07_Cross_species_comparison/results/shared_downregulated_orthogroups_new.txt")

# Nematostella genes in shared upregulated orthogroups
shared_up_nvec_genes <- orthos %>%
  filter(Orthogroup %in% shared_up_orthos) %>%
  select(Orthogroup, Nvec) %>%
  distinct()

# Stylophora genes in shared upregulated orthogroups
shared_up_styl_genes <- orthos %>%
  filter(Orthogroup %in% shared_up_orthos) %>%
  select(Orthogroup, Spi) %>%
  distinct()


# Over-representation analysis  -------------------------------------------

# Over representation analysis to compare enriched processes in Stylophora vs. Nematostella

# Vector for Stylophora upregulated genes
Stylophora_df<- orthogroups_de %>% filter(log2FoldChange>1 & padj<0.5)
stylophora_vect<- Stylophora_df  %>% filter(log2FoldChange > 1 & padj< 0.05) %>% pull(Nvec)

# For Nematostella cluster 1
clust1<- read.csv("~/immune_cells/cross_species/mast_results.csv")
clust1_up<- clust1 %>% filter(avg_log2FC > 1 & p_val_adj < 0.05) %>% pull(X)
clust1_vect<- gsub("-","_",clust1_up)


gene_lists <- list(
  Nematostella = clust1_vect,
  Stylophora = stylophora_vect
)

#Load GO annotations 
TermName= read.csv("~/immune_cells/cnidarian_immune_cells/08_supplementary_files/Goterms_NewAnnotation.csv")
TermGene= read.csv("~/immune_cells/cnidarian_immune_cells/08_supplementary_files/GoName_NewAnnotation.csv")

library(clusterProfiler)
# Run GO enrichment for both
cc_result <- compareCluster(geneCluster = gene_lists, fun = "enricher", pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, TERM2GENE = TermGene, TERM2NAME = TermName)

library(dplyr)
dotplot(cc_result, showCategory = 10, font.size = 12) + ggtitle("GO Term Comparison: Nematostella vs Stylophora")


# Get shared and unique genes 
# Upregulated gene sets (padj < 0.05)
nvec_up <- orthogroups_de %>%
  filter(avg_log2FC > 1, p_val_adj < 0.05) %>%
  pull(Nvec) %>%
  unique()

styl_up <- orthogroups_de %>%
  filter(log2FoldChange > 1, padj < 0.05) %>%
  pull(Spi) %>%
  unique()

# Downregulated gene sets
nvec_down <- orthogroups_de %>%
  filter(avg_log2FC < -1, p_val_adj < 0.05) %>%
  pull(Nvec) %>%
  unique()

styl_down <- orthogroups_de %>%
  filter(log2FoldChange < -1, padj < 0.05) %>%
  pull(Spi) %>%
  unique()

shared_up <- intersect(nvec_up, orthogroups_de$Nvec[orthogroups_de$Spi %in% styl_up])
nvec_specific_up <- setdiff(nvec_up, shared_up)

length(shared_up)
length(nvec_specific_up)

# Merge with gene names 
gene_names_df<- readRDS("~/immune_cells/scRNAseq_analysis/annotaion/peptides_annotation.rds")
gene_names_df$gene_name<- gsub("-","_", gene_names_df$gene_name)

gene_names_shared<- gene_names_df %>% filter(gene_name %in% shared_up)
gene_names_shared$conserved<- "yes"

gene_names_specific<- gene_names_df %>% filter(gene_name %in% nvec_specific_up)
gene_names_specific$conserved<- "no"

comparison_styl_nvec<- rbind(gene_names_shared,gene_names_specific)

write.csv(comparison_styl_nvec ,"~/immune_cells/cnidarian_immune_cells/07_Cross_species_comparison/results/shared_vs_specific.csv")

# Convergence and divergence analysis 
# Focus only on DE genes in Nematostella
nvec_de_filtered <- orthogroups_de %>%
  filter(p_val_adj < 0.05, !is.na(avg_log2FC)) %>%
  mutate(
    classification = case_when(
      padj < 0.05 & avg_log2FC * log2FoldChange > 0 ~ "Conserved",
      padj < 0.05 & avg_log2FC * log2FoldChange < 0 ~ "Divergent",
      is.na(log2FoldChange) | padj >= 0.05 ~ "Nematostella-specific"
    )
  )

# Convert row names of target_df to a column
gene_names_df<- readRDS("~/immune_cells/scRNAseq_analysis/annotaion/peptides_annotation.rds")

# Merge the dataframes
nvec_de_filtered <- merge(
  nvec_de_filtered,           
  gene_names_df,           
  by.x = "Nvec",      
  by.y = "gene_name",      
  all.x = TRUE,            
  sort = FALSE             
)

# If any is conserved then conserved, if any is diverged then diverged 

orthogroup_classified <- orthogroups_de %>%
  filter(!is.na(avg_log2FC), !is.na(log2FoldChange)) %>%
  mutate(
    classification = case_when(
      p_val_adj < 0.05 & padj < 0.05 & avg_log2FC * log2FoldChange > 0 ~ "Conserved",
      p_val_adj < 0.05 & padj < 0.05 & avg_log2FC * log2FoldChange < 0 ~ "Divergent",
      p_val_adj < 0.05 & (is.na(padj) | padj >= 0.05) ~ "Nematostella-only",
      padj < 0.05 & (is.na(p_val_adj) | p_val_adj >= 0.05) ~ "Stylophora-only",
      TRUE ~ "Non-DE or NA"
    )
  )

nvec_gene_summary <- orthogroup_classified %>%
  filter(p_val_adj < 0.05) %>%
  group_by(Nvec) %>%
  summarise(
    classification = case_when(
      any(classification == "Conserved") ~ "Conserved",
      any(classification == "Divergent") ~ "Divergent",
      TRUE ~ "Nematostella-specific"
    )
  )

nvec_gene_summary <- nvec_gene_summary %>%
  left_join(nvec_de %>% mutate(gene_id = gsub("-", "_", X)), by = c("Nvec" = "gene_id"))

write.csv(nvec_gene_summary ,"~/immune_cells/cnidarian_immune_cells/07_Cross_species_comparison/results/orthogroup_classified.csv")


