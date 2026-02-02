library(Seurat)
library(Matrix)
library(ggplot2)
library(ggpubr)
library(MAST)
library(scCustomize)
library(tidyverse)
# Load counts
counts_matrix <- readRDS("~/immune_cells/scRNAseq_analysis/Expression_matrix/Nvec_sc_raw_counts.rds")
counts_matrix <- as.matrix(counts_matrix)

# Create Seurat object
seurat_object <- CreateSeuratObject(counts = counts_matrix)

# Read metadata tables
cellid_to_mc <- read.delim(
  "~/immune_cells/scRNAseq_analysis/Expression_matrix/CellID_to_MC_table.txt",
  header = FALSE
)
exp_class <- read.delim(
  "~/immune_cells/scRNAseq_analysis/Expression_matrix/Nvec_sc_experiment_classification.txt",
  header = FALSE
)

colnames(cellid_to_mc) <- c("CellName", "CellID")
colnames(exp_class)    <- c("CellName", "Condition")

# Merge metadata (by cell barcode)
meta <- merge(cellid_to_mc, exp_class, by = "CellName", all = FALSE)

# Keep only cells that exist in Seurat object
meta <- meta[meta$CellName %in% colnames(seurat_object), ]

# Put CellName as rownames so Seurat can align safely
rownames(meta) <- meta$CellName
meta$CellName <- NULL

# Make sure I am not missing anything unexpectedly
stopifnot(all(rownames(meta) %in% colnames(seurat_object)))

# Add once (Seurat will match by rownames)
seurat_object <- AddMetaData(seurat_object, metadata = meta)

stopifnot(!anyDuplicated(rownames(meta)))

stopifnot(!anyNA(meta$CellID), !anyNA(meta$Condition))

table(table(rownames(meta)))

unique(seurat_object@meta.data$Condition)

#Replace with the correct conditions

seurat_object$Condition <- dplyr::recode(seurat_object$Condition,
                                         "iHCl" = "NaCl",
                                         "tPIC" = "pIC")

seurat <- seurat_object
# Pre-processing and QC ---------------------------------------------------

setwd("~/immune_cells/scRNAseq_analysis/seurat_pipeline/")
#saveRDS(seurat, "seurat_unprocessed.rds")
seurat <- readRDS("seurat_unprocessed.rds")
# Visualize QC metrics as a violin plot
VlnPlot(seurat,
        features = c("nFeature_RNA", "nCount_RNA"),
        ncol = 2)

# nCount_RNA vs nFeature_RNA: depth–complexity relationship per cell - looks good.
# nFeature_RNA violin: gene detection comparability across identities - looks good.

# Assess relationship between sequencing depth (UMI counts) and transcriptome complexity (detected genes) per cell.
# Check if there is any difference between the biological replicates

FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


# Overall, looks clean. Minimal filtering is needed.


#-----------------------------#
# QC filtering thresholds
#-----------------------------#
MIN_FEATURES <- 500
MAX_FEATURES <- 3000
MAX_COUNTS   <- 10000

MIN_CELLS_PER_GENE <- 10
REMOVE_PATTERN     <- "orphan"   # Remove "orphan" genes which are irrelevant for this project.

#-----------------------------#
# Data before filtering
#-----------------------------#
qc_before <- seurat@meta.data[, c("nFeature_RNA", "nCount_RNA")]
qc_before$stage <- "Before"

#-----------------------------#
# Cell filtering
#-----------------------------#
seurat_filt <- subset(seurat,
                      subset = nFeature_RNA > MIN_FEATURES &
                        nFeature_RNA < MAX_FEATURES &
                        nCount_RNA   < MAX_COUNTS)

#-----------------------------#
# Gene filtering
#-----------------------------#
counts_mat <- GetAssayData(seurat_filt, assay = "RNA", layer = "counts")

keep_by_detection <- Matrix::rowSums(counts_mat > 0) >= MIN_CELLS_PER_GENE
keep_by_name      <- !grepl(REMOVE_PATTERN, rownames(seurat_filt), ignore.case = TRUE)

seurat_filt <- seurat_filt[keep_by_detection & keep_by_name, ]

stopifnot(!any(grepl(
  REMOVE_PATTERN, rownames(seurat_filt), ignore.case = TRUE
)))

#-----------------------------#
# Data after filtering
#-----------------------------#
qc_after <- seurat_filt@meta.data[, c("nFeature_RNA", "nCount_RNA")]
qc_after$stage <- "After"

qc_all <- rbind(qc_before, qc_after)

# Quick summaries
cat("Cells before:", nrow(qc_before), "\n")
cat("Cells after :", nrow(qc_after), "\n")
cat("Genes before:", nrow(seurat), "\n")
cat("Genes after :", nrow(seurat_filt), "\n")

#-----------------------------#
# Plots with thresholds
#-----------------------------#

# Scatter: nCount vs nFeature (Before/After) + threshold lines
p_scatter <- ggplot(qc_all, aes(x = nCount_RNA, y = nFeature_RNA)) +
  geom_point(size = 0.4, alpha = 0.35) +
  facet_wrap( ~ stage) +
  geom_vline(xintercept = MAX_COUNTS, linetype = "dashed") +
  geom_hline(yintercept = MIN_FEATURES, linetype = "dashed") +
  geom_hline(yintercept = MAX_FEATURES, linetype = "dashed") +
  theme_classic() +
  labs(title = "QC space before/after filtering", x = "nCount_RNA (UMIs per cell)", y = "nFeature_RNA (genes detected per cell)")

# Violin: nFeature_RNA + thresholds
p_vln_feat <- ggplot(qc_all, aes(x = stage, y = nFeature_RNA)) +
  geom_violin() +
  geom_jitter(width = 0.2,
              size = 0.25,
              alpha = 0.25) +
  geom_hline(yintercept = MIN_FEATURES, linetype = "dashed") +
  geom_hline(yintercept = MAX_FEATURES, linetype = "dashed") +
  theme_classic() +
  labs(title = "nFeature_RNA before/after", x = NULL, y = "nFeature_RNA")

# Violin: nCount_RNA + threshold
p_vln_count <- ggplot(qc_all, aes(x = stage, y = nCount_RNA)) +
  geom_violin() +
  geom_jitter(width = 0.2,
              size = 0.25,
              alpha = 0.25) +
  geom_hline(yintercept = MAX_COUNTS, linetype = "dashed") +
  theme_classic() +
  labs(title = "nCount_RNA before/after", x = NULL, y = "nCount_RNA")

# Print plots
p_scatter
p_vln_feat
p_vln_count

# Save filtered object
#saveRDS(seurat_filt, "seurat_filt.RDS")


# Normalization and transformation ----------------------------------------

# I used the sctransform R package (https://github.com/satijalab/sctransform) to preform normalization and transformation

# I ran the following command in R on the cluster:

# library(sctransform)

# seurat <- SCTransform(
#  seurat,
#  verbose = FALSE
# )

# Then I saved the SCT object:
# saveRDS(seurat, "seurat_SCT.rds")


# Dimension reduction and clustering  -------------------------------------

seurat_SCT <- readRDS("~/immune_cells/scRNAseq_analysis/seurat_pipeline/seurat_SCT.rds")

# These are now standard steps in the Seurat workflow for visualization and clustering
seurat <- RunPCA(seurat_SCT, verbose = FALSE)
# Inspect the elbow plot
ElbowPlot(seurat, ndims = 50)
# PC 1 to 30 explain most of the variance

seurat <- RunUMAP(seurat, dims = 1:30, verbose = FALSE)

seurat <- FindNeighbors(seurat, dims = 1:30, verbose = FALSE)

# Louvain clustering
# I intentionally used low resolution to avoid over-clustering and to focus on main cell types in early embryos
# This allowed the detection of neurons, gland cells, and even a tiny cnidocyte population

seurat <- FindClusters(seurat, verbose = FALSE, resolution = 0.2)

# Visualize
DimPlot(seurat, label = TRUE)
scCustomize::DimPlot_scCustom(
  seurat,
  label = T,
  figure_plot = T,
  colors_use = "Dark2",
  shuffle = T
)
# Save the clustered object
#saveRDS(seurat,
#        "~/immune_cells/scRNAseq_analysis/seurat_pipeline/seurat_clust.RDS")


#QC
DefaultAssay(seurat) <- "RNA"
FeaturePlot(seurat, features = c("nCount_RNA", "nFeature_RNA"))

# It looks relatively uniform:
# The UMAP structure was not shaped by sequencing depth or detected gene counts.
# I will proceed without further filtering.


# Analysis and visualization  ---------------------------------------------

# Read the clustered object

seurat <- readRDS("~/immune_cells/scRNAseq_analysis/seurat_pipeline/seurat_clust.RDS")

# Visualize replicates
# Extract replicate information from the cell names
seurat$Replicate <- ifelse(grepl("^Nvec01_", colnames(seurat)), "Replicate 1", "Replicate 2")

# Visualize UMAP, coloring cells by the 'Replicate' column
DimPlot(seurat, group.by = "Replicate")

# No batch effect - the two replicates are nearly identical.

# Visualize the distribution of cells by conditions
DimPlot(seurat, group.by = "Condition")

# poly(I:C) has a strong enrichment of cells in cluster 1

# Visualize metacells

DimPlot(seurat, group.by = "CellID")

# Difficult to see. Another method:

meta <- seurat@meta.data %>%
  select(CellID, seurat_clusters)

# contingency table
mc_cluster_tab <- table(meta$CellID, meta$seurat_clusters)

mc_frac <- as.data.frame(mc_cluster_tab)
colnames(mc_frac) <- c("CellID", "Cluster", "n_cells")

mc_frac <- mc_frac %>%
  group_by(CellID) %>%
  mutate(total_cells = sum(n_cells),
         frac = n_cells / total_cells) %>%
  ungroup()

mc_cluster1 <- mc_frac %>%
  filter(Cluster == "1") %>%
  arrange(desc(frac))

head(mc_cluster1, 20)

# Mostly the metacells from 141 and above which fits Arnau's analysis.

# Visualize expression of genes of interest

FeaturePlot(seurat, features = "mCherry-plus-strand")

VlnPlot(seurat,
        features = "mCherry-plus-strand",
        pt.size = 0.2,
        split.by = "Condition")

library(scCustomize)


# RlRb_1 "Nvec-vc1.1-XM-048731783.1"
FeaturePlot_scCustom(seurat, features = "Nvec-vc1.1-XM-048731783.1")

# RLRb_2 Nvec-vc1.1-XM-048731783.1
FeaturePlot_scCustom(seurat, features = "Nvec-vc1.1-XM-048731786.1")

# Visualize mCherry-plus-strand expression on UMAP grouped by 'Condition'
FeaturePlot_scCustom(seurat, features = "mCherry-plus-strand", split.by = "Condition")

# RLRb_1
FeaturePlot_scCustom(seurat, features = "Nvec-vc1.1-XM-048731783.1", split.by = "Condition")

# Check the correlation in the metacell expression matrix (uploaded to Zenodo).
f <- "~/working/scRNAseq Arnau/to_yehu/expression_matrices/Nvec_mc_umifrac.txt"
tab <- read.delim(header = T, f)

rows <- tab[c("mCherry_plus_strand", "Nvec_vc1.1_XM_048731783.1"), ]

df <- as.data.frame(t(rows))

# Create scatter plot with linear regression line and R-squared value
plot <- ggscatter(
  df,
  x = "Nvec_vc1.1_XM_048731783.1",
  y = "mCherry_plus_strand",
  add = "reg.line",
  conf.int = TRUE,
  cor.coef = TRUE,
  cor.method = "pearson",
  ggtheme = theme_pubr()
) +
  stat_cor(aes(label = paste(..rr.label.., sep = "~`, `~")), label.x = 0.1, label.y = 2) +
  labs(title = "", x = "RLRb UMI fraction", y = "mCherry UMI fraction") +
  theme(text = element_text(size = 12))

# Save as PNG
# ggsave("scatter_plot_mc.png", plot, width = 7, height = 5, dpi = 300)

print(plot)

# RLRb_1 against RLRb_2

rows <- tab[c("Nvec_vc1.1_XM_048731783.1", "Nvec_vc1.1_XM_048731786.1"), ]

df <- as.data.frame(t(rows))

# Create scatter plot with linear regression line and R-squared value
plot <- ggscatter(
  df,
  x = "Nvec_vc1.1_XM_048731783.1",
  y = "Nvec_vc1.1_XM_048731786.1",
  add = "reg.line",
  conf.int = TRUE,
  cor.coef = TRUE,
  cor.method = "pearson",
  ggtheme = theme_pubr()
) +
  stat_cor(aes(label = paste(..rr.label.., sep = "~`, `~")), label.x = 0.1, label.y = 2) +
  labs(title = "", x = "RLRb_1 UMI fraction", y = "RLRb_2 UMI fraction") +
  theme(text = element_text(size = 12))

print(plot)

# RLRb against GAPDH - as a sanity check

rows <- tab[c("Nvec_vc1.1_XM_048731783.1", "Nvec_vc1.1_XM_032382055.2"), ]

df <- as.data.frame(t(rows))

# Create scatter plot with linear regression line and R-squared value
plot <- ggscatter(
  df,
  x = "Nvec_vc1.1_XM_048731783.1",
  y = "Nvec_vc1.1_XM_032382055.2",
  add = "reg.line",
  conf.int = TRUE,
  cor.coef = TRUE,
  cor.method = "pearson",
  ggtheme = theme_pubr()
) +
  stat_cor(aes(label = paste(..rr.label.., sep = "~`, `~")), label.x = 0.1, label.y = 2) +
  labs(title = "", x = "RLRb_1 UMI fraction", y = "GAPDH UMI fraction") +
  theme(text = element_text(size = 12))

print(plot)

# No correlation.

# Visualize the expression of mCherry per condition

gene_expression <- FetchData(seurat, vars = c("mCherry-plus-strand", "Condition"))

# Create a boxplot
ggplot(gene_expression,
       aes(x = Condition, y = `mCherry-plus-strand`, fill = Condition)) +
  geom_boxplot(outlier.shape = NA) +  # Boxplot without outliers
  geom_jitter(
    width = 0.2,
    alpha = 0.6,
    color = "black",
    size = 0.5
  ) +  # Add jitter for individual cells
  theme_minimal() +
  labs(title = "Expression of GeneX across Conditions", x = "Condition", y = "Expression Level") +
  scale_fill_brewer(palette = "Set2")
unique(gene_expression$Condition)

# Cleaner version
ggplot(gene_expression,
       aes(x = Condition, y = `mCherry-plus-strand`, fill = Condition)) +
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(method = "wilcox.test",
                     label = "p.signif",
                     comparisons = list(c("iHCl", "Ctrl"), c("tPIC", "iHCl"), c("tPIC", "Ctrl"))) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "none",
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  labs(title = "GeneX Expression Across Conditions", x = "Condition", y = "Expression Level") +
  scale_fill_manual(values = c("#1f78b4", "#33a02c", "#e31a1c"))

# Visualize co-expression of RLRb and mCherry
umap_data <- FetchData(seurat, vars = c("umap_1", "umap_2", "mCherry-plus-strand", "Nvec-vc1.1-XM-048731783.1"))

umap_data$`mCherry-plus-strand` <- scales::rescale(umap_data$`mCherry-plus-strand`, to = c(0, 1))
umap_data$`Nvec-vc1.1-XM-048731783.1` <- scales::rescale(umap_data$`Nvec-vc1.1-XM-048731783.1`, to = c(0, 1))

ggplot(umap_data, aes(x = umap_1, y = umap_2)) +
  geom_point(aes(color = rgb(`mCherry-plus-strand`,`Nvec-vc1.1-XM-048731783.1`,0)), size = 1) +
  scale_color_identity() +
  theme_minimal() +
  labs(title = "Co-expression of mCherry and RLRb", x = "UMAP 1", y = "UMAP 2")


# Determine cell type proportions

# Create a contingency table of clusters and conditions
cluster_condition_table <- table(seurat$seurat_clusters, seurat$Condition)

# Convert to proportions
cluster_condition_proportions <- prop.table(cluster_condition_table, margin = 2)

# Plot using ggplot2
library(ggplot2)
data <- as.data.frame(as.table(cluster_condition_proportions))
colnames(data) <- c("Cluster", "Condition", "Proportion")

ggplot(data, aes(x = Condition, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_minimal() +
  labs(y = "Proportion", x = "Condition", fill = "Cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Activated versus control molecular profile ------------------------------
seurat <- readRDS("~/immune_cells/scRNAseq_analysis/seurat_pipeline/seurat_clust.RDS")

# Subset cells

Idents(seurat) <- "Condition"
subset_seurat <- subset(seurat, idents = c("iHCl", "tPIC"))

# Use MAST for differential expression analysis between poly(I:C) ("tPIC") cells vsrsus
# NaCl cells ("iHCl).
#"MAST" : Identifies differentially expressed genes between two groups of cells using a hurdle model tailored to scRNA-seq data. Utilizes the MAST package to run the DE testing.
# This can take a few minutes to run

#de_results <- FindMarkers(
#  seurat,
#  ident.1 = "tPIC",
#  ident.2 = "iHCl",
#  test.use = "MAST"
#)

# Save the MAST resuls
#saveRDS(de_results, "~/immune_cells/scRNAseq_analysis/seurat_pipeline/MAST_pIC_NaCl.rds")

de_results<- readRDS("~/immune_cells/scRNAseq_analysis/seurat_pipeline/MAST_pIC_NaCl.rds")
# Obtain annotations and load results

gene_names_df <- readRDS("~/immune_cells/scRNAseq_analysis/annotaion/peptides_annotation.rds")
de_results$gene_name <- rownames(de_results)
de_results <- readRDS("~/immune_cells/scRNAseq_analysis/seurat_pipeline/MAST_pIC_NaCl.rds")
de_results$gene_name <- rownames(de_results)

# Merge the dataframes
merged_df <- merge(
  de_results,
  gene_names_df,
  by.x = "gene_name",
  by.y = "gene_name",
  all.x = TRUE,
  sort = FALSE
)

# Restore row names
rownames(merged_df) <- merged_df$gene_name
merged_df$gene_name <- NULL

# View the result
head(merged_df)

# Visualize top genes using a heatmap
top_genes <- rownames(de_results[de_results$p_val_adj < 0.05 &
                                   abs(de_results$avg_log2FC) > 2, ])
DoHeatmap(subset_seurat, features = top_genes) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# Volcano plot
df <- merged_df
library(tidyverse)

filtered_df <- df %>%
  group_by(protein) %>%
  filter(abs(avg_log2FC) == max(abs(avg_log2FC))) %>%
  ungroup() %>% na.omit()


colors <- rep("#BFBFBF", nrow(filtered_df))
names(colors) <- rep("Negative", nrow(filtered_df))

colors[which(filtered_df$avg_log2FC >= 1 &
               filtered_df$p_val_adj < 0.05)] <- "#E69F00"
names(colors) [which(filtered_df$avg_log2FC >= 1 &
                       filtered_df$p_val_adj < 0.05)] <- "UP"

colors[which(filtered_df$avg_log2FC <= -1 &
               filtered_df$p_val_adj < 0.05)] <- "#990099"
names(colors) [which(filtered_df$avg_log2FC <= -1 &
                       filtered_df$p_val_adj < 0.05)] <- "DOWN"

library(EnhancedVolcano)

volcano_plot <- EnhancedVolcano(
  filtered_df,
  lab = filtered_df$protein,
  x = "avg_log2FC",
  y = "p_val_adj",
  pCutoff = 0.05,
  FCcutoff = 1,
  colCustom = colors,
  colAlpha = 0.3,
  title = "",
  subtitle = "",
  caption = "",
  drawConnectors = TRUE,
  widthConnectors = 0.75,
  ylim = c(NA, 100),
  max.overlaps = 15,
  labSize = 3.0
)

# Save the plot as a PDF
ggsave(
  filename = "~/immune_cells/cnidarian_immune_cells/03_scRNAseq/02_seurat/scRNA-seq_analysis_files/figure-markdown_strict/MAST_volcano_corrected.png",
  # File name
  plot = volcano_plot,
  # The plot object
  device = "png",
  # File format
  width = 8,
  height = 6,
  # Dimensions in inches
  dpi = 96                      # Resolution (not relevant for PDFs but ensures sharpness)
)


# Immune cluster characterization  ----------------------------------------
# =============================================================================
# Cluster 1 markers, visualization, and enrichment analysis
# =============================================================================

# ---- Libraries ----
library(Seurat)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)

# ---- Paths ----
seurat_rds <- "~/immune_cells/scRNAseq_analysis/seurat_pipeline/seurat_clust.RDS"
go_dir     <- "~/immune_cells/scRNAseq_analysis/annotaion/GOseq"
out_dir    <- "~/immune_cells/mCherry_RLRb_FACS/Differential_expression/figs_publication"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load Seurat object ----
seurat_obj <- readRDS(seurat_rds)

# =============================================================================
# 1) Differential expression: Cluster 1 markers
# =============================================================================
cluster_id <- 1

cluster1_markers <- FindMarkers(
  object = seurat_obj,
  ident.1 = cluster_id,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Quick peek at top markers
head(cluster1_markers)

# =============================================================================
# 2) Heatmap: top markers for Cluster 1
# =============================================================================
n_top_markers <- 50
top_markers <- rownames(cluster1_markers)[seq_len(min(n_top_markers, nrow(cluster1_markers)))]

heatmap_plot <- DoHeatmap(
  object   = seurat_obj,
  features = top_markers,
  group.by = "seurat_clusters",
  size     = 4
) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  ggtitle(sprintf("Top %d Markers for Cluster %s", length(top_markers), cluster_id)) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title     = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "right",
    legend.title    = element_blank()
  )

print(heatmap_plot)

#ggsave(
#  filename = file.path(out_dir, "cluster1_markers_heatmap.png"),
#  plot     = heatmap_plot,
#  width    = 8,
#  height   = 6,
#  dpi      = 300
#)

# =============================================================================
# 3) Expression visualization: selected immune genes across clusters
# =============================================================================
genes_of_interest <- c(
  "mCherry-plus-strand",
  "Nvec-vc1.1-XM-048731786.1",
  "Nvec-vc1.1-XM-032361976.2",
  "Nvec-vc1.1-XM-048730161.1"
)

vln_plot <- VlnPlot(
  object   = seurat_obj,
  features = genes_of_interest,
  pt.size  = 0,
  group.by = "seurat_clusters",
  ncol     = 1,
  combine  = TRUE
) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(vln_plot)

# =============================================================================
# 5) Prepare ranked gene list (for GSEA-style inputs)
# =============================================================================
gene_list <- cluster1_markers$avg_log2FC
names(gene_list) <- gsub("-", "_", rownames(cluster1_markers))

gene_list <- na.omit(gene_list)
gene_list <- sort(gene_list, decreasing = TRUE)

# =============================================================================
# 7) Load GO annotation tables (TERM2GENE / TERM2NAME)
# =============================================================================
term2gene_path <- file.path(go_dir, "GoName_NewAnnotation.csv")
term2name_path <- file.path(go_dir, "Goterms_NewAnnotation.csv")

TermGene <- read.csv(term2gene_path, header = TRUE, check.names = FALSE)
TermName <- read.csv(term2name_path, header = TRUE, check.names = FALSE)

stopifnot(is.data.frame(TermGene), is.data.frame(TermName))

# ---- Over-representation analysis ----
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# Filter significant upregulated genes
keep <- cluster1_markers$p_val_adj < 0.05 & cluster1_markers$avg_log2FC > 0.5
genes <- rownames(cluster1_markers[keep, , drop = FALSE])
genes <- gsub("-", "_", genes)
genes <- na.omit(genes)

# Load annotation tables
go_dir <- "~/immune_cells/scRNAseq_analysis/annotaion/GOseq"
TermGene <- read.csv(file.path(go_dir, "GoName_NewAnnotation.csv"), header = TRUE, check.names = FALSE)
TermName <- read.csv(file.path(go_dir, "Goterms_NewAnnotation.csv"), header = TRUE, check.names = FALSE)

# Run ORA
Results <- enricher(
  gene      = genes,
  TERM2GENE = TermGene,
  TERM2NAME = TermName,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.2,
  minGSSize = 10,
  maxGSSize = 500
)

# Visualize 
p = barplot(Results, showCategory = 15)

# Rotate the labels by 45 degrees
p + theme(axis.text.y = element_text(angle = 45, hjust = 1))

# Save table 
write.table(
  as.data.frame(Results),
  file = "~/immune_cells/mCherry_RLRb_FACS/Differential_expression/figs_publication/mCherry_upregulated_ORA.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# Finding activation and basal condition specific markers within the immune cluster

# ==========================================================
# Cluster 1 "immune" DE (tPIC vs Ctrl) + GO-based labeling
# Assign each DE gene to TF / Receptor / Effector / Other
# - Produce: (1) DotPlot of top genes per category, (2) compareCluster ORA
# ==========================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(clusterProfiler)
})

# -------------------------
# Inputs 
# -------------------------
seurat_obj <- readRDS("~/immune_cells/scRNAseq_analysis/seurat_pipeline/seurat_clust.RDS")
gene_names_df <- readRDS("~/immune_cells/scRNAseq_analysis/annotaion/peptides_annotation.rds")

TermGene <- read.csv(
  file = "~/immune_cells/scRNAseq_analysis/annotaion/GOseq/GoName_NewAnnotation.csv",
  header = TRUE,
  check.names = FALSE
)

TermName <- read.csv(
  file = "~/immune_cells/scRNAseq_analysis/annotaion/GOseq/Goterms_NewAnnotation.csv",
  header = TRUE,
  check.names = FALSE
)

# -------------------------
# Parameters I am interested in
# -------------------------
cluster_id   <- "1"
cond1        <- "tPIC"
cond2        <- "Ctrl"

# DE settings (FindMarkers)
logfc_thresh <- 0.25
min_pct      <- 0.10

# What you call "strong" DE for downstream summaries
de_fc_cutoff <- 0.5
padj_cutoff  <- 0.05

# How many representatives to show per category
top_n        <- 3

# Outputs (same as your original intent)
out_tpic_csv <- "~/immune_cells/figures_revised/tpic_up_df.csv"
out_ctrl_csv <- "~/immune_cells/figures_revised/ctrl_up_df.csv"
out_portrait <- "~/immune_cells/figures_revised/molecular_portrait.csv"

# -------------------------
# Helper: pick top N genes per category (by absolute effect size)
# -------------------------
get_top_genes <- function(df, direction, n = 3) {
  df %>%
    group_by(category) %>%
    slice_max(order_by = abs(avg_log2FC), n = n, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(condition = direction)
}

# ==========================================================
# 1) Subset the immune cluster and run DE between conditions
# ==========================================================
immune_cells <- subset(seurat_obj, idents = cluster_id)

# Use Condition as the identity class (so FindMarkers compares conditions)
Idents(immune_cells) <- immune_cells$Condition

de_genes <- FindMarkers(
  immune_cells,
  ident.1 = cond1,
  ident.2 = cond2,
  logfc.threshold = logfc_thresh,
  min.pct = min_pct
)

# ==========================================================
# 2) Build a compact "gene -> GO term names" lookup table
#    (so each gene gets a single text field we can keyword-scan)
# ==========================================================
gene_go <- TermGene %>%
  inner_join(TermName, by = "Go_Term") %>%
  transmute(
    gene = Ids,   # TERM2GENE IDs (underscore format)
    GO_name = Name
  ) %>%
  group_by(gene) %>%
  summarise(
    GO_combined = paste(unique(GO_name), collapse = "; "),
    .groups = "drop"
  )

# ==========================================================
# 3) Attach GO text + peptide annotations to the DE table
#    Then label each gene into coarse functional buckets.
#
# ==========================================================
de_annot <- de_genes %>%
  tibble::rownames_to_column("gene_raw") %>%
  mutate(
    # Normalize for GO join: GO uses underscores
    gene = gsub("-", "_", gene_raw)) %>%
  left_join(gene_go, by = "gene") %>%
  mutate(
    category = dplyr::case_when(
      grepl("transcription factor|DNA-binding|regulation of transcription",
            GO_combined, ignore.case = TRUE) ~ "TF",
      
      grepl("receptor activity|signaling receptor|G protein-coupled receptor|membrane",
            GO_combined, ignore.case = TRUE) ~ "Receptor",
      
      grepl("immune|defense|cytokine|interferon|effector",
            GO_combined, ignore.case = TRUE) ~ "Effector",
      
      TRUE ~ "Other"
    )
  ) %>%
  left_join(gene_names_df, by = c("gene_raw" = "gene_name"))

message("Category counts:")
print(table(de_annot$category, useNA = "ifany"))

message("Peptide annotation coverage:")
message("  annotated rows: ", sum(!is.na(de_annot$protein) | !is.na(de_annot$domain), na.rm = TRUE))
message("  total rows:     ", nrow(de_annot))
# ==========================================================
# 4) Define the "up" gene sets
#    (strong log2foldchange + significant adjusted p-value)
# ==========================================================
tpic_up <- de_annot %>%
  filter(avg_log2FC >  de_fc_cutoff, p_val_adj < padj_cutoff)

ctrl_up <- de_annot %>%
  filter(avg_log2FC < -de_fc_cutoff, p_val_adj < padj_cutoff)

# Save tables
#write.csv(tpic_up, out_tpic_csv, row.names = FALSE)
#write.csv(ctrl_up, out_ctrl_csv, row.names = FALSE)

# ==========================================================
# 5) Quick visual sanity check: DotPlot of representative genes
#    We show the strongest hits per category for each direction.
# ==========================================================
top_genes <- bind_rows(
  get_top_genes(tpic_up, cond1, n = top_n),
  get_top_genes(ctrl_up, cond2, n = top_n)
) %>%
  arrange(condition, category)

# DotPlot needs Seurat's feature names -> use gene_raw
genes_to_plot <- unique(top_genes$gene_raw)

p_dot <- DotPlot(immune_cells, features = genes_to_plot, group.by = "Condition") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = paste0("Cluster ", cluster_id, ": top DE genes by GO-category (", cond1, " vs ", cond2, ")"))

print(p_dot)

# ==========================================================
# 6) ORA: compare GO enrichment between the two "up" sets
#    Here I use underscore-normalized IDs to match TERM2GENE.
# ==========================================================
gene_lists <- list(
  tPIC_up = unique(tpic_up$gene),
  Ctrl_up = unique(ctrl_up$gene)
)

cc_result <- compareCluster(
  geneCluster = gene_lists,
  fun = "enricher",
  TERM2GENE = TermGene,
  TERM2NAME = TermName,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2
)

p_cc <- dotplot(cc_result, showCategory = 10, font.size = 12) +
  ggtitle("GO term enrichment comparison: tPIC-up vs Ctrl-up")

print(p_cc)

# ==========================================================
# 7) Save the full annotated DE table ("molecular portrait")
# ==========================================================
write.csv(de_annot, out_portrait, row.names = FALSE)

message(
  "Finished.\n",
  "Wrote:\n",
  "  - ", out_tpic_csv, "\n",
  "  - ", out_ctrl_csv, "\n",
  "  - ", out_portrait
)


# Comparison with the bulk-RNA-seq results --------------------------------
# Compute module score for the mCherry positive cells and visualize the bulk-RNA-seq.
# Where the mCherry+ genes are expressed?
seurat<- readRDS("~/immune_cells/scRNAseq_analysis/seurat_pipeline/seurat_clust.RDS")
# Load the mCherry+ upregulated genes
gene_set <- read.delim(
  "~/immune_cells/mCherry_RLRb_FACS/Differential_expression/mCherry_upregulated.txt"
)  # Replace the the names from underscore to dash
gene_set <- rownames(gene_set)
gene_set <- gsub("_", "-", gene_set)
# Check if genes exist in the dataset
gene_set <- intersect(gene_set, rownames(seurat_obj))

# Add module score for the gene set
seurat_obj <- AddModuleScore(seurat, features = list(gene_set), name = "ProcessScore")

# Visualize the module score on UMAP
aggregated_score <- FeaturePlot_scCustom(
  seurat_obj,
  features = "ProcessScore1",
  reduction = "umap",
  pt.size = 0.5,
  colors_use = viridis_inferno_dark_high
) +  ggtitle("") + theme(axis.title.x = element_blank(), axis.title.y = element_blank())


# Display the plot
print(aggregated_score)

# Save the plot as high-resolution output
ggsave(
  "aggregated_mCherry_score_umap.png",
  plot = aggregated_score,
  width = 8,
  height = 6,
  dpi = 300
)

# Identification of marker genes and annotation ---------------------------
# Cell type markers S6
# All marker genes
cellType_markers <- read.csv(
  "~/immune_cells/scRNAseq_analysis/seurat_pipeline/figs_publication/supplementary/cellType_markers3.csv"
)

# Read the clustered Seurat object
#seurat<- readRDS("~/immune_cells/scRNAseq_analysis/seurat_pipeline/seurat_clust.RDS")
# Add immune genes


# Replace underscores with dashes in the entire data frame
df_markers <- data.frame(lapply(cellType_markers, function(x) {
  if (is.character(x) || is.factor(x)) {
    gsub("_", "-", as.character(x)) # Replace in character or factor columns
  } else {
    x # Leave other columns unchanged
  }
}), stringsAsFactors = FALSE)


df_markers <- df_markers[1:17, ]

# Plotting

library(patchwork)

visualize_markers <- function(seurat, df_markers) {
  plots <- lapply(seq_len(nrow(df_markers)), function(i) {
    gene <- df_markers$DToL[i]
    model <- df_markers$ID[i]
    name <- df_markers$Marker[i]
    
    # Check if the gene exists in the Seurat object
    if (gene %in% rownames(seurat)) {
      FeaturePlot_scCustom(seurat, features = gene) +
        ggtitle(paste(name, ":", model)) +
        theme(plot.title = element_text(size = 8)) +
        theme(axis.title = element_blank())
    } else {
      message(paste("Skipping", gene, "as it is not found in the Seurat object"))
      NULL
    }
  })
  
  # Remove NULL entries from the plot list
  plots <- Filter(Negate(is.null), plots)
  
  # Combine plots if any are available
  if (length(plots) > 0) {
    wrap_plots(plots, ncol = 4)
  } else {
    message("No valid markers found for visualization.")
  }
}

visualize_markers(seurat = seurat, df_markers = df_markers)

S6 <- visualize_markers(seurat = seurat, df_markers = df_markers)


ggsave(
  "S6_marker_genes.png",
  plot = S6,
  width = 12,
  height = 9,
  dpi = 300
)



# Annotate cell types and re-print figure

Idents(seurat) <- "seurat_clusters"

table(Idents(seurat))

# Convert to characters

Idents(seurat) <- factor(as.character(seurat$seurat_clusters))
levels(Idents(seurat))

cell_types <- c(
  "0" = "Ectoderm",
  "1" = "Immune",
  "2" = "Ect.oral",
  "3" = "Nemtocye",
  "4" = "Ect.aboral",
  "5" = "Gland",
  "6" = "Neuronal"
)

cell_types_num <- setNames(paste0(names(cell_types), ": ", cell_types), names(cell_types))

seurat <- RenameIdents(seurat, cell_types_num)


table(Idents(seurat))
levels(Idents(seurat))

DimPlot(seurat)


# Save as RDS file (available on Zenodo)
#saveRDS(
#  seurat,
#  "~/immune_cells/scRNAseq_analysis/seurat_pipeline/seurat_obj2_annotated.rds"
#)


# Visualizing the WGCNA metacell results  ---------------------------------

library(Seurat)
library(dplyr)

seurat<- readRDS("~/immune_cells/scRNAseq_analysis/seurat_pipeline/seurat_clust.RDS")
modules<- read.delim("/sci/labs/yehum79/itamar273/scRNAseq Arnau/to_yehu/gene_modules/WGCNA_gmod_annotation.txt")

# seurat: Clustered Seurat object
# modules:  WGCNA data.frame with columns gene_id, gene_module, membership_score, ...

# Make sure clustering identities are the active identities
Idents(seurat) <- "seurat_clusters"


# Apply the same order as the WGCNA compact heatmap (Supplementary Fig. 7a)

gs_order <- c(
  "turquoise", "pink", "grey60", "salmon", "cyan",
  "grey", "green", "yellow", "tan", "greenyellow",
  "blue", "magenta", "black", "red", "purple",
  "lightcyan", "brown", "lightgreen", "midnightblue"
)

GS_names <- paste0("GS", seq_along(gs_order))

modules$gene_id<- gsub("_","-", modules$gene_id)

library(dplyr)

gene_sets_df <- modules %>%
  filter(gene_module %in% gs_order) %>%
  filter(gene_id %in% rownames(seurat)) %>%
  distinct(gene_module, gene_id) %>%
  group_by(gene_module) %>%
  summarise(genes = list(gene_id), .groups = "drop")

gene_sets <- setNames(gene_sets_df$genes, gene_sets_df$gene_module)

gene_sets <- gene_sets[gs_order]

stopifnot(
  length(gene_sets) == 19,
  all(names(gene_sets) == gs_order)
)


seurat <- AddModuleScore(
  object   = seurat,
  features = gene_sets,
  name     = "GS_",
  assay    = DefaultAssay(seurat),
  slot     = "data"
)

meta <- seurat[[]]

old <- paste0("GS_", seq_along(gs_order))
new <- paste0("GS",  seq_along(gs_order))

colnames(meta)[match(old, colnames(meta))] <- new
seurat[[]] <- meta

all(new %in% colnames(seurat[[]]))  # TRUE



vln_plot <- VlnPlot(
  seurat,
  features = new,
  group.by = "seurat_clusters",
  pt.size = 0.1,
  ncol = 7,
  alpha = 0.3,
  cols = rep("grey",7),
) 

print(vln_plot)

# Conduct over-representation analysis on individual gene sets

library(dplyr)
library(clusterProfiler)

# ---- Load annotation once ----
setwd("~/immune_cells/scRNAseq_analysis/annotaion/GOseq/")

TermGene <- read.csv("GoName_NewAnnotation.csv", header = TRUE, check.names = FALSE)
TermName <- read.csv("Goterms_NewAnnotation.csv", header = TRUE, check.names = FALSE)

out_dir  <- "~/immune_cells/Supp_tables/ORA_GS14_GS17/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


stopifnot(is.data.frame(TermGene), is.data.frame(TermName))

# ---- Define which modules correspond to GS14–GS17 ----
gs_to_module <- c(
  GS14 = "red",
  GS15 = "purple",
  GS16 = "lightcyan",
  GS17 = "brown"
)

modules$gene_id<- gsub("-","_", modules$gene_id)

# ---- Run ORA + barplot for each ----
ora_results <- list()
ora_plots   <- list()

for (gs in names(gs_to_module)) {
  
  mod <- gs_to_module[[gs]]
  
  genes <- modules %>%
    filter(gene_module == mod) %>%
    pull(gene_id) %>%
    unique()
  
  res <- enricher(
    gene          = genes,
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.2,
    minGSSize     = 10,
    maxGSSize     = 500,
    TERM2GENE     = TermGene,
    TERM2NAME     = TermName
  )
  
  ora_results[[gs]] <- res
  
  # only make a plot if enrichment exists
  if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
    ora_plots[[gs]] <- barplot(
      res,
      showCategory = 15,
      title = paste0(gs, " (", mod, ")")
    )
  } else {
    ora_plots[[gs]] <- NULL
  }
}

print(ora_plots$GS14)
print(ora_plots$GS15)
print(ora_plots$GS16)
print(ora_plots$GS17)

sessionInfo()
