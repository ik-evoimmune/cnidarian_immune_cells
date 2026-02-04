setwd(
  "~/immune_cells/cnidarian_immune_cells/01_Immuno_FACS_bulkRNAseq/02_differential_expression_analysis/"
)
#Read count matrix
f <- "data/immuno_count_matrix.txt"
count_matrix <- read.delim(f,
                           header = T,
                           row.names = 1,
                           skip = 1)
colnames(count_matrix) <- c(
  "Chr",
  "Start",
  "End",
  "Strand",
  "Length",
  "RLRb_low1",
  "RLRb_low2",
  "RLRb_low3",
  "RLRb_low4",
  "RLRb_high1",
  "RLRb_high2",
  "RLRb_high3",
  "RLRb_high4"
)
library(dplyr)
colnames(count_matrix)
count_matrix <- count_matrix %>% dplyr::select(-Chr, -Start, -End, -Strand)
countData <- as.matrix(count_matrix %>% dplyr::select(-Length))

#create coldata
SampleName <- colnames(countData)
Condition <- c(rep("RLRb_low", 4), rep("RLRb_high", 4))
colData <- data.frame(SampleName, Condition)
rownames(colData) <- colnames(countData)

all(rownames(colData) %in% colnames(countData))
all(rownames(colData) == colnames(countData))

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ Condition)
#Filtering
keep <- rowSums(counts(dds)) >= 50
table(keep)
dds <- dds[keep, ]

#Running the analysis
dds <- DESeq(dds)
res <- results(dds)
# Define the contrasts correctly
res <- results(dds, contrast = c("Condition", "RLRb_high", "RLRb_low"))

res <- na.omit(res)
keep <- res$padj < 0.05 & res$log2FoldChange > 1
res_filt <- res[keep, ]
nrow(res_filt)
write.table(res_filt, "data/RLRb_high_upregulated_filt.txt", sep = "\t")
sum(res$padj < 0.05 & res$log2FoldChange > 1 , na.rm = T)
# 2048 upregulated
sum(res$padj < 0.05 & res$log2FoldChange < -1 , na.rm = T)
# 384 downregulated
write.csv(res, file = "S2_table.csv", row.names = TRUE)

# Exploratory andlysis
# PCA plot
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "Condition")

# Assuming 'vsd' is your variance-stabilized or rlog-transformed DESeq2 object
pcaData <- plotPCA(vsd, intgroup = "Condition", returnData = TRUE)

# Calculate percentage of variance explained by each principal component
percentVar <- round(100 * attr(pcaData, "percentVar"))
# Save files:
saveRDS(percentVar, "data/percentVar.rds")
saveRDS(pcaData, "data/pcaData.rds")

# Add homology table
peptides_anotation <- read.csv("../../09_supplementary_files/Dictionary_Kozlovski_et.al.csv")
res_df <- as.data.frame(res)
rownames(res_df) <- gsub("-", "_", rownames(res_df))
res_df$gene_name <- rownames(res_df)

# Heatmap
library(ComplexHeatmap)
ntd <- normTransform(dds, f = log2, pc = 1)
# Show only the top 100 genes
top.genes <- order(rowVars(assay(ntd)), decreasing = TRUE)[1:100]
top.ntd <- ntd[top.genes, ]
assay(top.ntd) <- assay(top.ntd) - rowMeans(assay(top.ntd))
col.anno <- HeatmapAnnotation(RLRb_level = ntd$Condition)

Heatmap(
  assay(top.ntd),
  show_row_names = FALSE,
  show_column_names = FALSE,
  name = "Expression",
  top_annotation = col.anno
)


# Rename the RLRs
#RLRa
target_gene <- "Nvec-vc1.1-XM-032364224.2"
new_protein_name <- "RLRa"
peptides_anotation$protein[peptides_anotation$gene_name == target_gene] <- new_protein_name

# RLRb
target_gene <- "Nvec-vc1.1-XM-048731783.1"
new_protein_name <- "RLRb"
peptides_anotation$protein[peptides_anotation$gene_name == target_gene] <- new_protein_name

# RLRb paralog
target_gene <- "Nvec-vc1.1-XM-048731786.1"
new_protein_name <- "RLRb_2"
peptides_anotation$protein[peptides_anotation$gene_name == target_gene] <- new_protein_name

# Merge the dataframes
merged_df1 <- merge(
  res_df,
  peptides_anotation,
  by.x = "gene_name",
  by.y = "DToL_gene_model",
  all.x = TRUE,
  sort = FALSE
)

# Restore row names
rownames(merged_df1) <- merged_df1$gene_name
merged_df1$gene_name <- NULL


# Save as table S2
write.csv(merged_df1, "data/S2_table.csv")

library(EnhancedVolcano)
EnhancedVolcano(merged_df1,
                lab = merged_df1$protein,
                x = 'log2FoldChange',
                y = 'padj')


# Take only one copy of a gene that is the most highly expressed
df <- merged_df1

filtered_df <- df %>%
  group_by(protein) %>%
  filter(abs(log2FoldChange) == max(abs(log2FoldChange))) %>%
  ungroup() %>% na.omit()

EnhancedVolcano(filtered_df,
                lab = filtered_df$protein,
                x = 'log2FoldChange',
                y = 'padj')

# Improved visualization
library(EnhancedVolcano)

colors <- rep("#C0C0C0", nrow(filtered_df))
names(colors) <- rep("Negative", nrow(filtered_df))

#Assign colors to down vs. up genes

colors[which(filtered_df$log2FoldChange >= 1 &
               filtered_df$padj < 0.05)] <- "#D55E00"
names(colors) [which(filtered_df$log2FoldChange >= 1 &
                       filtered_df$padj < 0.05)] <- "UP"

colors[which(filtered_df$log2FoldChange <= -1 &
               filtered_df$padj < 0.05)] <- "#0072B2"
names(colors) [which(filtered_df$log2FoldChange <= -1 &
                       filtered_df$padj < 0.05)] <- "DOWN"

volcano_plot <- EnhancedVolcano(
  filtered_df,
  lab = filtered_df$protein,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 1,
  selectLab = c(
    "GBP3",
    "IRF8",
    "IFI30",
    "IRF2",
    "RLRb",
    "OAS1",
    "IFI44",
    "GBP6",
    "FOS",
    "TRIM56",
    "TRIM45",
    "CD38",
    "TRAF4",
    "JUN",
    "MAFA",
    "CLEC4A",
    "PTPRD",
    "LBP",
    "TRAF4"
  ),
  colCustom = colors,
  colAlpha = 0.3,
  title = "",
  subtitle = "",
  caption = "",
  drawConnectors = TRUE,
  widthConnectors = 0.75,
  ylim = c(NA, 150),
  max.overlaps = 33,
  labSize = 3.0
)

# Save the plot as a PDF
ggsave(
  filename = "data/volcano_plot1_RLRb.pdf",
  plot = volcano_plot,
  device = "pdf",
  width = 8,
  height = 6,
  dpi = 300
)

library(EnhancedVolcano)

# Generate the plot
volcano_plot2 <- EnhancedVolcano(
  filtered_df,
  lab = filtered_df$protein,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 1,
  selectLab = c(
    "GBP3",
    "IRF8",
    "IFI30",
    "IRF2",
    "RLRb",
    "OAS1",
    "IFI44",
    "GBP6",
    "FOS",
    "TRIM56",
    "TRIM45",
    "CD38",
    "TRAF4",
    "JUN",
    "MAFA",
    "CLEC4A",
    "PTPRD",
    "LBP",
    "TRAF4"
  ),
  colCustom = colors,
  colAlpha = 0.2,
  title = NULL,
  subtitle = NULL,
  caption = NULL,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = "grey30",
  boxedLabels = TRUE,
  labSize = 3.5,
  max.overlaps = Inf,
  pointSize = 2.0,
  axisLabSize = 12,
  labCol = "black",
  labFace = "plain",
  ylim = c(NA, 150)
)

ggsave(
  filename = "data/volcano_plot2.pdf",
  plot = volcano_plot2,
  device = "pdf",
  width = 7.5,
  height = 6,
  useDingbats = FALSE
)




# Gene set enrichment analysis  -------------------------------------------

merged_df = merged_df1
dim(merged_df)
head(merged_df)

# we want the log2 fold change
original_gene_list <- merged_df1$log2FoldChange

# name the vector
names(original_gene_list) <- rownames(merged_df)
head(original_gene_list)
# omit any NA values
gene_list <- na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
# Loading Nematostella annotation files
setwd("~/immune_cells/scRNAseq_analysis/annotaion/GOseq/")
list.files()
TermGene  <- read.csv(file = 'GoName_NewAnnotation.csv',
                      header = TRUE,
                      check.names = FALSE)
head(TermGene)
is.data.frame(TermGene)

TermName  <- read.csv(file = 'Goterms_NewAnnotation.csv',
                      header = TRUE,
                      check.names = FALSE)
head(TermName)
is.data.frame(TermName)

#Running GSEA
library(clusterProfiler)
gse <- GSEA(
  gene_list,
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  TERM2GENE = TermGene,
  TERM2NAME = TermName,
  verbose = TRUE,
  seed = FALSE
)
#Save table
Results_GSEA <- as.data.frame(gse)
write.csv(Results_GSEA,
          file = "~/immune_cells/cnidarian_immune_cells/01_Immuno_FACS_bulkRNAseq/02_differential_expression_analysis/data/S3_table.csv",
          quote = FALSE,
          row.names = FALSE)

library(enrichplot)
# Representative processes
gseaplot2(
  gse,
  geneSetID = c(
    "GO:2001235",
    "GO:0006956",
    "GO:0030414",
    "GO:0006955",
    "GO:0035821"
  )
)


library(clusterProfiler)
library(enrichplot)

# Optimized GSEA plot
gseaplot2(
  gse,
  geneSetID = c(
    "GO:2001235",
    "GO:0006956",
    "GO:0030414",
    "GO:0006955",
    "GO:0035821"
  ),
  base_size = 16,
  title = "Selected GO Enrichments in GSEA",
  subplots = 1:3,
  rel_heights = c(1.5, 0.5, 1),
  color = c("#D55E00", "#0072B2", "#009E73", "#F0E442", "#CC79A7"),
  pvalue_table = TRUE
)

# Save the table of immune processes
library(gridExtra)
library(grid)
library(Cairo)

# Extract selected terms with additional metrics
selected_ids <- c("GO:2001235",
                  "GO:0006956",
                  "GO:0030414",
                  "GO:0006955",
                  "GO:0035821")

table_df <- gse@result[gse@result$ID %in% selected_ids, c("Description", "NES", "pvalue", "p.adjust", "setSize")]

# Clean and format
table_df$pvalue <- signif(table_df$pvalue, 3)
table_df$p.adjust <- signif(table_df$p.adjust, 3)
table_df$NES <- round(table_df$NES, 2)

# Rename for clarity
colnames(table_df) <- c("GO Term", "NES", "p-value", "Adjusted p-value", "Gene Set Size")


library(gridExtra)
library(grid)
library(Cairo)

# Build table grob with left alignment and larger font
table_grob <- tableGrob(
  table_df,
  rows = NULL,
  theme = ttheme_default(
    core = list(fg_params = list(
      hjust = 0,
      x = 0.05,
      fontsize = 11
    )),
    colhead = list(fg_params = list(
      fontface = "bold", fontsize = 12
    )),
    padding = unit(c(5, 5), "mm")
  )
)

# Widen left column if clipped
table_grob$widths[1] <- max(unit(4, "cm"), table_grob$widths[1])

# Save as high-quality PNG and PDF
CairoPNG(
  "GSEA_metrics_table.png",
  width = 2000,
  height = 700,
  res = 300
)
grid.newpage()
grid.draw(table_grob)
dev.off()

CairoPDF("GSEA_metrics_table.pdf",
         width = 16,
         height = 2.8)
grid.newpage()
grid.draw(table_grob)
dev.off()


# RT-qPCR validation ------------------------------------------------------

library(pcr)
#Loat CT values
res_tab <- read.csv(
  "~/immune_cells/cnidarian_immune_cells/01_Immuno_FACS_bulkRNAseq/02_differential_expression_analysis/data/CT_table.csv"
)

## add grouping variable
group_var <- rep(c('Positive', 'Negative'), each = 3)
# calculate all values and errors in one step
## mode == 'separate_tube' default

res1 <- pcr_analyze(res_tab,
                    group_var = group_var,
                    reference_gene = 'HKG',
                    reference_group = "Negative")
res1
library(ggpubr)

# Run t-test

tst1 <- pcr_test(res_tab,
                 group_var = group_var,
                 reference_gene = 'HKG',
                 reference_group = 'Negative',
                 test = 't.test')
knitr::kable(tst1, caption = 'pIC injection: t-test summary')

library(dplyr)
library(ggplot2)
library(forcats)

# --- Significance stars (trim whitespace just in case) ---
tst1_sig <- tst1 %>%
  mutate(sig = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE ~ "ns"
  )) %>%
  mutate(sig = trimws(sig)) %>%
  select(gene, estimate, sig)

# Order genes by |effect| 
gene_order <- tst1_sig %>% arrange(desc(abs(estimate))) %>% pull(gene)

sum_df <- res1 %>%
  transmute(group, gene, relative_expression,
            sem = error) %>%
  left_join(tst1_sig %>% select(gene, sig), by = "gene") %>%
  mutate(
    gene  = factor(gene,  levels = gene_order),
    group = factor(group, levels = c("Negative","Positive"))
  )

# y-position for stars (slightly above the tallest mean+SEM per gene)
star_df <- sum_df %>%
  group_by(gene) %>%
  summarise(
    y_star = max(relative_expression + sem, na.rm = TRUE) * 1.09,
    sig    = dplyr::first(sig),
    .groups = "drop"
  )

pd <- position_dodge(width = 0.55)

p <- ggplot(sum_df, aes(x = gene, y = relative_expression, color = group, fill = group)) +
  geom_point(position = pd, size = 3.5, alpha = 0.95) +
  geom_errorbar(aes(ymin = relative_expression - sem, ymax = relative_expression + sem),
                position = pd, width = 0.16, linewidth = 1) +
  geom_text(data = star_df, aes(x = gene, y = y_star, label = sig),
            inherit.aes = FALSE, size = 6, vjust = 0) +
  # Custom stronger colors 
  scale_color_manual(values = c("Negative" = "#4D4D4D", "Positive" = "#D73027")) +
  scale_fill_manual(values = c("Negative" = "#4D4D4D", "Positive" = "#D73027")) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(x = NULL, y = "Relative expression (ΔΔCt)", color = NULL, fill = NULL) +
  theme_classic(base_size = 16) +   # bigger font
  theme(
    legend.position = "top",
    legend.box      = "horizontal",
    axis.text.x     = element_text(angle = 40, hjust = 1),
    axis.line       = element_line(linewidth = 0.7),
    plot.margin     = margin(6, 14, 12, 6)
  ) +
  coord_cartesian(clip = "off")

p

sessionInfo()
