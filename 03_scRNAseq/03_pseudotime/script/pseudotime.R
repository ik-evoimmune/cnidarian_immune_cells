# Trajectory inference across conditions: differential expression and differential progression
suppressPackageStartupMessages({
  library(slingshot)
  library(SingleCellExperiment)
  library(RColorBrewer)
  library(scales)
  library(viridis)
  library(UpSetR)
  library(pheatmap)
  library(msigdbr)
  library(fgsea)
  library(knitr)
  library(ggplot2)
  library(gridExtra)
  library(tradeSeq)
})

# Load data
seurat<- readRDS("~/immune_cells/scRNAseq_analysis/seurat_pipeline/seurat_clust.RDS")
sce<- Seurat::as.SingleCellExperiment(seurat, assay = "RNA")

# Take 1:49 dimensions to avoid a bug in slingshot 
reducedDim(sce, 'PCA_corrected') <- reducedDim(sce, 'PCA')[, 1:49]


# Re-clustering of the immune cluster (cluster 1) -------------------------



# Reclustering of cluster 1 
library(Seurat)
cluster_of_interest <- "1"  
subset_obj <- subset(seurat, idents = cluster_of_interest)

subset_obj <- NormalizeData(subset_obj)
subset_obj <- FindVariableFeatures(subset_obj)
subset_obj <- ScaleData(subset_obj)

subset_obj <- RunPCA(subset_obj)
ElbowPlot(subset_obj)
subset_obj <- RunUMAP(subset_obj, dims = 1:15)  

subset_obj <- FindNeighbors(subset_obj, dims = 1:15)
subset_obj <- FindClusters(subset_obj, resolution = 0.2)  
DimPlot(subset_obj, reduction = "pca", label = TRUE)

# Renaming the clusters

# Extract current identities
current_clusters <- Idents(subset_obj)

#Define the mapping (e.g., 2 → 0, 0 → 1, 1 → 2, 3 → 3)
new_cluster_mapping <- c("2" = "0", "0" = "1", "1" = "2", "3" = "3")

#Rename the clusters based on the mapping
new_clusters <- as.character(new_cluster_mapping[as.character(current_clusters)])

#Update the identities in the Seurat object
Idents(subset_obj) <- factor(new_clusters, levels = unique(new_clusters))

# Number of cells per cluster 
table(Idents(subset_obj))

# PCA plot of re-clustering 
library(ggplot2)
library(RColorBrewer)


# Plotting 
# Reorder cluster identities
subset_obj <- subset_obj
Idents(subset_obj) <- factor(Idents(subset_obj), levels = c("0", "1", "2", "3"))

# Generate the plot
library(RColorBrewer)

p <- DimPlot(
  subset_obj, 
  reduction = "pca", 
  label = TRUE,
  label.size = 6,                
  pt.size = 1.5,                 
  cols = brewer.pal(4, "Set2"),  
  alpha = 0.5                   
)

# Customize theme
p <- p + theme_minimal() + 
  theme(
    axis.text = element_text(size = 14),       
    axis.title = element_text(size = 16),     
    legend.title = element_text(size = 14),   
    legend.text = element_text(size = 12)   
  )

# Display the plot
print(p)



# Finding markers for each sub-cluster 

# Find markers for each sub-cluster compared to all other sub-clusters
markers <- FindAllMarkers(subset_obj, 
                          group.by = "seurat_clusters",
                          min.pct = 0.25,             
                          logfc.threshold = 0.25)     

head(markers)

# Add the protein annotation 
gene_names_df<- readRDS("~/immune_cells/scRNAseq_analysis/annotaion/peptides_annotation.rds")
markers$gene_name <- rownames(markers)

# Merge the dataframes
merged_df <- merge(
  markers,               
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
write.csv(merged_df, "~/immune_cells/cnidarian_immune_cells/03_scRNAseq/03_pseudotime/results/S11_table.csv")
# Visualize top markers 
# Top 10 markers for visualization
library(tidyverse)
top10_markers <- merged_df %>% group_by(cluster) %>% top_n(10, wt = avg_log2FC)

# Plot the top markers
DotPlot(subset_obj, features = top10_markers$gene) + 
  RotatedAxis()  # Rotate axis labels for better readability

# Visualizing a heatmap for the top markers
DoHeatmap(subset_obj, features = top10_markers$gene) + 
  scale_fill_viridis()  # Optional, for a nicer color scale

# Reshape the data into a tidy format suitable for ggplot2
dotplot_data <- merged_df %>%
  select(protein, cluster, avg_log2FC, pct.1) %>%  
  mutate(cluster = factor(cluster))  

library(ggplot2)
library(dplyr)


#Filter proteins with non-empty annotations and pct.1 > 0.1
filtered_data <- merged_df %>%
  filter(!is.na(protein) & protein != "") %>%  # Keep only rows with non-empty protein annotations
  filter(pct.1 > 0.1)  # Only keep proteins expressed in >10% of cells

#Select top 10 markers per cluster based on avg_log2FC
top_markers <- filtered_data %>%
  group_by(cluster) %>%
  top_n(10, wt = avg_log2FC) %>%
  ungroup()  # Remove grouping after filtering

#Create the DotPlot using ggplot2
ggplot(top_markers, aes(x = protein, y = cluster, size = pct.1, color = avg_log2FC)) +
  geom_point() +                              
  scale_size(range = c(2, 6)) +                
  scale_color_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  theme_minimal() +                            
  theme(axis.text.x = element_text(angle = 90, hjust = 1),  
        axis.text = element_text(size = 12),    
        axis.title = element_text(size = 14),   
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 10))  

# Optimized 

library(dplyr)
library(ggplot2)
library(viridis)

#Filter proteins with non-empty annotations and pct.1 > 0.1
filtered_data <- merged_df %>%
  filter(!is.na(protein) & protein != "") %>%
  filter(pct.1 > 0.1)

#Select top 10 markers per cluster based on avg_log2FC
top_markers <- filtered_data %>%
  group_by(cluster) %>%
  top_n(10, wt = avg_log2FC) %>%
  ungroup()

#Order proteins within each cluster by avg_log2FC for better grouping
top_markers <- top_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  mutate(protein_ordered = factor(protein, levels = unique(protein))) %>%
  ungroup()

#Sort the clusters for better control of y-axis order
top_markers$cluster <- factor(top_markers$cluster, levels = sort(unique(top_markers$cluster)))

# Select top proteins
top_proteins <- top_markers$protein %>% unique()

# Filter full data for those proteins across all clusters
full_marker_data <- filtered_data %>%
  filter(protein %in% top_proteins)

# Reorder factor levels for plotting
full_marker_data <- full_marker_data %>%
  mutate(protein_ordered = factor(protein, levels = unique(top_markers$protein)))

# Plot
ggplot(full_marker_data, aes(x = protein_ordered, y = factor(cluster), size = pct.1, color = avg_log2FC)) +
  geom_point() +
  scale_size(range = c(2, 6)) +
  scale_color_viridis(option = "plasma", direction = -1) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title = element_blank(),
    panel.grid.major = element_line(color = "grey85")
  ) +
  guides(
    color = guide_colorbar(title = "avg_log2FC"),
    size = guide_legend(title = "pct.1")
  )


# Pseudotime analysis  ----------------------------------------------------

# Slingshot trajectory analysis using the subset 
sce<- as.SingleCellExperiment(subset_obj, assay = "RNA")

library(SingleCellExperiment)

# Extract the correct cluster identities from the Seurat object
correct_clusters <- as.character(Idents(subset_obj))

# Assign these identities to the SCE object
sce$ident <- factor(correct_clusters, levels = c("0", "1", "2", "3"))


# Verify the changes
table(sce$ident)


reducedDim(sce, 'PCA_corrected') <- reducedDim(sce, 'PCA')[, 1:49]


slingshot <- slingshot::slingshot(
  sce, 
  reducedDim = 'PCA_corrected', 
  clusterLabels = sce$seurat_clusters
)
slingshot::slingLineages(slingshot)

# Fetching pseudotime
sce$pseudotime <- slingshot::slingPseudotime(slingshot)[, 'Lineage1']

# Fetching principal curve in PCA space
pca_curve <- slingshot::slingCurves(slingshot, as.df = TRUE)
colnames(pca_curve) <- paste0('PC', 1:ncol(pca_curve))
head(pca_curve)

umap_curve <- slingshot::embedCurves(slingshot, 'UMAP', smoother = 'loess', span = 0.1) |> 
  slingshot::slingCurves(as.df = TRUE)
umap_curve <- umap_curve[order(umap_curve$Order), ]
head(umap_curve)

# Plotting
library(tibble)
library(ggplot2)
df <- tibble(
  PC1 = reducedDim(sce, 'PCA_corrected')[,1], 
  PC2 = reducedDim(sce, 'PCA_corrected')[,2], 
  UMAP1 = reducedDim(sce, 'UMAP')[,1], 
  UMAP2 = reducedDim(sce, 'UMAP')[,2], 
  annotation = sce$ident, 
  pseudotime = sce$pseudotime,
  condition = sce$Condition
)

# Plotting paramaters p <- DimPlot(subset_obj, reduction = "pca", label = TRUE, 
# label.size = 6,              # Increase label size
# pt.size = 1.5,               # Adjust point size
# cols = brewer.pal(4, "Set2"), # Use the Set2 palette for 4 clusters
# alpha = 0.5)     

cowplot::plot_grid(
  df |> 
    ggplot() + 
    geom_point(aes(PC1, PC2, col = annotation)) + 
    geom_path(data = pca_curve, aes(x = PC1, y = PC2)) + 
    theme_bw() + 
    coord_fixed(),
  df |> 
    ggplot() + 
    geom_point(aes(UMAP1, UMAP2, col = annotation)) + 
    geom_path(data = umap_curve, aes(x = umap_1, y = umap_2)) + 
    theme_bw() + 
    coord_fixed(),
  df |> 
    ggplot() + 
    geom_point(aes(PC1, PC2, col = pseudotime)) + 
    geom_path(data = pca_curve, aes(x = PC1, y = PC2)) + 
    theme_bw() + 
    coord_fixed(),
  df |> 
    ggplot() + 
    geom_point(aes(UMAP1, UMAP2, col = pseudotime)) + 
    geom_path(data = umap_curve, aes(x = umap_1, y = umap_2)) + 
    theme_bw() + 
    coord_fixed()
)

df |>
  ggplot() +
  geom_point(aes(PC1, PC2, col = annotation)) +
  geom_path(data = pca_curve, aes(x = PC1, y = PC2)) +
  theme_bw() +
  coord_fixed() +
  facet_wrap(~condition)

# fig 6D optimized

library(ggplot2)
library(viridis)

p<- df |> 
  ggplot(aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = pseudotime), size = 1, alpha = 0.8) +         # Tighter point size and transparency
  geom_path(data = pca_curve, aes(x = PC1, y = PC2), color = "black", linewidth = 0.8) + 
  scale_color_viridis(name = "Pseudotime", option = "viridis") +      # Better color scale
  coord_fixed() +
  theme_minimal(base_size = 14) +                                     # Cleaner theme
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

print(p)

ggsave("~/immune_cells/figures_revised/pseudotime_pca_plot.png", plot = p, width = 6, height = 8, units = "in", dpi = 600)


# Match the style to the previous plots
library(ggplot2)
library(RColorBrewer)

# Define the color palette
set2_palette <- brewer.pal(4, "Set2")

# Update the ggplot2-based plot to match DimPlot appearance
df |>
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(col = annotation), size = 1.5, alpha = 0.5) +
  geom_path(data = pca_curve, aes(x = PC1, y = PC2)) +        
  scale_color_manual(values = set2_palette) +                  
  theme_minimal() +                                            
  theme(
    axis.text = element_text(size = 14),       
    axis.title = element_text(size = 16),   
    legend.title = element_text(size = 14),   
    legend.text = element_text(size = 12)     
  ) +
  coord_fixed()  


# By conditions plot
library(ggplot2)
library(RColorBrewer)

# Define the color palette
set2_palette <- brewer.pal(4, "Set2")

# Update the faceted ggplot2 plot to match DimPlot appearance
p2<- df |>
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(col = annotation), size = 1.5, alpha = 0.5) +  
  geom_path(data = pca_curve, aes(x = PC1, y = PC2)) +    
  scale_color_manual(values = set2_palette) +                  
  theme_minimal() +                                            
  theme(
    axis.text = element_text(size = 14),       
    axis.title = element_text(size = 16),     
    legend.title = element_text(size = 14),   
    legend.text = element_text(size = 12),    
    strip.text = element_text(size = 14)      
  ) +
  coord_fixed() +  
  facet_wrap(~condition)  

print(p2)

ggsave("~/immune_cells/figures_revised/pseudotime_by_condition_pca_plot.png", plot = p2, width = 6, height = 8, units = "in", dpi = 600)

# Normalize the data
# Compute size factors for normalization
sce<- scuttle::logNormCounts(sce) 

# Ensure 'ident' is a factor or integer vector
sce$ident <- factor(sce$ident, levels = c(0, 1, 2, 3))

# Replace with new labels
levels(sce$ident) <- c("1", "2", "0", "1")
library(tidyverse)
# Genes of interest as a function of pseudotime
# This order: mCherry, MPEG1, RLRa, IFI44L
genes <- c('mCherry-plus-strand',"Nvec-vc1.1-XM-048733868.1","Nvec-vc1.1-XM-032385537.2","Nvec-vc1.1-XM-001625557.3")
fitExprs <- logcounts(sce[genes, ]) |> # ----------------------------------- Get norm. counts for genes of interest
  as.matrix() |> 
  t() |> 
  as_tibble() |> 
  mutate(  # ----------------------------------------------------------------- Add information for each cell
    cellID = colnames(sce), 
    annotation = factor(sce$ident), 
    pseudotime = sce$pseudotime
  ) |> 
  pivot_longer(contains(genes), names_to = 'gene', values_to = 'obs_expr') |> # - Pivot in "long" tidy format 
  mutate(gene = factor(gene, genes)) |> 
  group_by(gene) |> # ------------------------------------------------------- Group rows by genes
  nest(.key = 'data') |> # -------------------------------------------------- For each gene, extract the subtable into a column named data
  mutate(
    gamModel = map(data, ~mgcv::gam(obs_expr ~ s(pseudotime, bs = "cs"), data = .)), 
    gamFitted_expr = map(gamModel, predict) # ------------------------------ For each gene, fit the expression values ~ pseudotime with a GAM
  ) |> 
  dplyr::select(-ends_with('Model')) |> 
  unnest(c(data, ends_with('_expr'))) # -------------------------------------- Unnest all the modelled expressions
ggplot(fitExprs) + 
  geom_point(aes(x = pseudotime, y = obs_expr, col = annotation), alpha = 0.5) + 
  geom_line(aes(x = pseudotime, y = gamFitted_expr), col = 'white', size = 2, alpha = 0.5) + 
  geom_line(aes(x = pseudotime, y = gamFitted_expr), col = '#af2d0c', size = 1) +
  theme_bw() + 
  facet_grid(gene~., scales = 'free') + 
  labs(y = 'logcounts') + 
  ggtitle('Fitted models of pseudotime-dependent gene expression')

setwd("~/immune_cells/scRNAseq_analysis/pseudotime_analysis/figs_publication/")


# Define the color palette
set2_palette <- brewer.pal(4, "Set2")

# Update the plot to match the DimPlot appearance
ggplot(fitExprs) +
  geom_point(aes(x = pseudotime, y = obs_expr, col = annotation), alpha = 0.5, size = 1.5) +  
  geom_line(aes(x = pseudotime, y = gamFitted_expr), col = 'white', size = 2, alpha = 0.5) +  
  geom_line(aes(x = pseudotime, y = gamFitted_expr), col = '#af2d0c', size = 1) +             
  scale_color_manual(values = set2_palette) +                                                
  theme_minimal() +                                                                         
  theme(
    axis.text = element_text(size = 14),       
    axis.title = element_text(size = 16),     
    legend.title = element_text(size = 14),   
    legend.text = element_text(size = 12),    
    strip.text = element_text(size = 14)      
  ) +
  facet_wrap(~gene, scales = 'free', ncol = 2) +      
  labs(
    y = 'logcounts',                           
    x = 'Pseudotime',                          
    title = 'Fitted models of pseudotime-dependent gene expression'
  )


# Over-representation analysis for the re-clustered data  -----------------


# Loading Nematostella annotation files
setwd("~/immune_cells/scRNAseq_analysis/annotaion/GOseq/")
list.files()
TermGene  <- read.csv(file = 'GoName_NewAnnotation.csv',header=TRUE, check.names=FALSE)
head(TermGene)
is.data.frame(TermGene)

TermName  <- read.csv(file = 'Goterms_NewAnnotation.csv',header=TRUE, check.names=FALSE)
head(TermName)
is.data.frame(TermName)

#Running ORA 
library(clusterProfiler)

# Visualize cluster 0 and 1 (ORA analysis) 

df1 <- merged_df %>% filter (cluster == 0) 
rownames(df1) <- gsub("-", "_", rownames(df1))
keep<- df1$p_val_adj <0.05 & df1$avg_log2FC> 0.5
# 129 genes

df1<- df1[keep,] 
clust0<- rownames(df1)
nrow(df1)
clust0<- rownames(df1)


df2 <- merged_df %>% filter (cluster == 1) 
rownames(df2) <- gsub("-", "_", rownames(df2))
keep<- df2$p_val_adj <0.05 & df2$avg_log2FC> 0.5
df2<- df2[keep,] 
nrow(df2)
clust1<- rownames(df2)

gene_lists <- list(
  Condition1 = clust0,
  Condition2 = clust1
)

library(clusterProfiler)

cc_result <- compareCluster(geneCluster = gene_lists, fun = "enricher", pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, TERM2GENE = TermGene, TERM2NAME = TermName)

dotplot(cc_result, showCategory = 8, font.size = 12) + ggtitle("Enrichment Comparison")

res1 <- enricher(gene_lists[[1]], pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, TERM2GENE = TermGene, TERM2NAME = TermName)
res2 <- enricher(gene_lists[[2]], pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, TERM2GENE = TermGene, TERM2NAME = TermName)

res1@result$group <- "Condition1"
res2@result$group <- "Condition2"

combined <- rbind(res1@result, res2@result)
top_terms <- combined %>%
  group_by(group) %>%
  slice_min(p.adjust, n = 10)  


library(ggplot2)

ggplot(top_terms, aes(x = -log10(p.adjust), y = fct_reorder(Description, p.adjust), color = group)) +
  geom_point(size = 3) +
  facet_wrap(~ group, scales = "free_y") +
  labs(x = "-log10 adjusted p-value", y = "GO Term", title = "GO Enrichment") +
  theme_bw()

getwd()
sink("sessionInfo.txt")
sessionInfo()
sink()
