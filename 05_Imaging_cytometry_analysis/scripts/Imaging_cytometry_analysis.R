
# PCA analysis  -----------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(ggpubr)

setwd("C:/Users/Itamar/OneDrive - huji.ac.il/Documents/Yehu Moran Antiviral Evolution/Imaging cytometry analysis/RLRb imaging CSV files/")


RLRb1<- read.csv("RLRb_exp_1.csv", header = T, skip = 1)
RLRb2<- read.csv("RLRb_exp_2.csv", header = T, skip = 1)
RLRb3<- read.csv("RLRb_exp_3.csv", header = T, skip = 1)
IgG1<- read.csv("IgG_exp_1.csv", header = T, skip = 1)
IgG2<- read.csv("IgG_exp_2.csv", header = T, skip = 1)
IgG3<- read.csv("IgG_exp_3.csv", header = T, skip = 1)


#Add a column for replicate number and status (IgG vs. RLRb)
RLRb1$status<- "RLRb1"
RLRb1$replicate<- 1

RLRb2$status<- "RLRb2"
RLRb2$replicate<- 2

RLRb3$status<- "RLRb3"
RLRb3$replicate<- 3

IgG1$status<- "IgG1"
IgG1$replicate<- 1

IgG2$status<- "IgG2"
IgG2$replicate<- 2

IgG3$status<- "IgG3"
IgG3$replicate<- 3


# Compare the to 15% of RLRb expressing cells against the rest

combined.rlrb<- rbind(RLRb1,RLRb2,RLRb3)
dim(combined.rlrb)
colnames(combined.rlrb)

#1: Split into top 15% vs the rest
threshold <- quantile(combined.rlrb$Intensity_MC_Ch02, 0.85) # Top 15% threshold
df <- combined.rlrb %>%
  mutate(Group = ifelse(Intensity_MC_Ch02 >= threshold, "Top 15%", "Rest"))
dim(df)
#2: Select numeric data and remove zero-variance columns
# Ensure `Group` is retained for alignment
numeric_data <- df %>%
  select(-status, -Intensity_MC_Ch02, -replicate) %>% na.omit() %>% 
  select(Group ,where(~ is.numeric(.) && var(.) > 0.01))

#3: Perform PCA on scaled data
pca_result <- prcomp(select(numeric_data, -Group), scale. = TRUE)

#4: Add PCA results and Group column
pca_df <- as.data.frame(pca_result$x) %>%
  bind_cols(Group = numeric_data$Group)

# Step 5: Visualize PCA by Group
ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.7, shape = 16) +  
  labs(
    title = "PCA Plot (Top 15% vs Rest)",
    x = "Principal Component 1",
    y = "Principal Component 2"
  ) +
  theme_minimal() +
  theme(legend.position = "top") +
  facet_wrap(~ Group, ncol = 1) +
  geom_jitter(width = 0.1, height = 0.1, alpha = 0.5, size = 2)

# Show only the centroids representing each of the replicates 

#1: Split data into Top 15% vs Rest for each status
df_split <- df %>%
  group_by(status) %>%
  mutate(
    Group = ifelse(Intensity_MC_Ch02 >= quantile(Intensity_MC_Ch02, 0.85), "Top 15%", "Rest")
  )

#2: Select numeric data and remove zero-variance columns

numeric_data <- df_split %>%
  select( -Intensity_MC_Ch02, -replicate, -Object.Number) %>% na.omit() %>% 
  select(Group, status, where(~ is.numeric(.) && var(.) > 0.01))
colnames(numeric_data)

numeric_data_clean <- numeric_data %>%
  ungroup() %>%
  select(-Group, -status)

#3 Standardize the data before PCA
numeric_data_scaled <- scale(numeric_data_clean)
dim(numeric_data_clean)

#4: Perform PCA on the scaled numeric data
pca_result <- prcomp(numeric_data_scaled, scale. = TRUE)

# Step 5: Add PCA results and status + Group columns for visualization
pca_df <- as.data.frame(pca_result$x) %>%
  bind_cols(Group = numeric_data$Group, status = numeric_data$status)

# Calculate centroids (means) for each group
centroids <- pca_df %>%
  group_by(Group, status) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = 'drop')

# Now plot the PCA with only the centroids (dots for the centroids)
ggplot() +
  geom_point(data = centroids, aes(x = PC1, y = PC2, color = Group), size = 4, shape = 16) +
  labs(x = "PC1", y = "PC2", title = "PCA - Group Centroids Only") +
  theme_minimal()


# Variance explained
# Extract variance explained
variances <- pca_result$sdev^2  # squared singular values (standard deviations)
variance_explained <- variances / sum(variances) * 100  # percentage variance explained

# Create PCA plot with variance explained on the axes
ggplot() +
  geom_point(data = centroids, aes(x = PC1, y = PC2, color = Group), size = 4, shape = 16) +
  labs(
    x = paste("PC1 (", round(variance_explained[1]), "%)", sep = ""),
    y = paste("PC2 (", round(variance_explained[2]), "%)", sep = ""),
    title = "PCA - Group Centroids Only"
  ) +
  theme_minimal()

# Define custom colors
centroids$Group
custom_colors <- c("Rest" = "cyan", "Top 15%" = "magenta") 

# Create the plot
pca_plot1 <- ggplot(data = centroids) +
  geom_point(
    aes(x = PC1, y = PC2, color = Group),
    size = 4,
    shape = 16
  ) +
  scale_color_manual(values = custom_colors) +
  labs(
    x = paste("PC1 (", round(variance_explained[1]), "%)", sep = ""),
    y = paste("PC2 (", round(variance_explained[2]), "%)", sep = ""),
    title = "Morphological Data",
    color = "Group"  
  ) +
  theme_minimal(base_size = 14) + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"), 
    axis.title = element_text(face = "bold"),             
    legend.position = "right"                             
  )

# Display the plot
print(pca_plot1)

# PCA for transcriptomic data
pcaData<- readRDS("pcaData.rds")
percentVar<- readRDS("percentVar.rds")
# Define custom colors for Conditions
custom_colors <- c("RLRb_low" = "cyan", "RLRb_high" = "magenta") # Replace with actual condition names

# Create PCA plot
pca_plot2 <- ggplot(pcaData, aes(x = PC1, y = PC2, label = rownames(pcaData), color = Condition)) +
  geom_point(size = 4, shape = 16) +                            # Plot PCA points with larger size
  scale_color_manual(values = custom_colors) +                  # Apply custom colors
  labs(
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance"),
    title = "Transcriptomic Data",
    color = "Condition"                                          # Add legend title
  ) +
  theme_minimal(base_size = 14) +                               # Use larger base font
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),      # Center and bold title
    axis.title = element_text(face = "bold"),                  # Bold axis labels
    legend.position = "right",                                 # Place legend on the right
    legend.title = element_text(face = "bold")                 # Bold legend title
  )

# Display the plot
print(pca_plot2)

# Plot them together in 1 column
library(patchwork)
# Combine both PCA plots into one column
combined_plot <- pca_plot1 / pca_plot2

# Display the combined plot
print(combined_plot)

# Extract the loadings
# Extract PCA loadings (eigenvectors)
loadings <- pca_result$rotation

# Compute the absolute values of the loadings for the first PC and second PC
pc1_loadings <- abs(loadings[, 1])  # Absolute values for PC1
pc2_loadings <- abs(loadings[, 2])  # Absolute values for PC2

# Rank features by their contribution to the first principal component (PC1)
pc1_contrib <- sort(pc1_loadings, decreasing = TRUE)
top_features_pc1 <- names(pc1_contrib)[1:20]  # Top 10 features contributing to PC1

# Rank features by their contribution to the second principal component (PC2)
pc2_contrib <- sort(pc2_loadings, decreasing = TRUE)
top_features_pc2 <- names(pc2_contrib)[1:20]  # Top 10 features contributing to PC2

# Print the top features for both PC1 and PC2
cat("Top 10 features contributing to PC1:\n", top_features_pc1, "\n")
cat("Top 10 features contributing to PC2:\n", top_features_pc2, "\n")


# Individual plots of IDEAS data --------

# Combine all data 

combined.data<- rbind(IgG1,IgG2,IgG3,RLRb1,RLRb2,RLRb3)

#Get the median intensity of antibody staining  for each condition 
dim(combined.data)
colnames(combined.data)

median.intensity <- combined.data |>
  na.omit() |>
  dplyr::filter(Intensity_MC_Ch02 > 0) |>
  dplyr::group_by(status, replicate) |>
  dplyr::summarise(
    median_fluoresence_intensity = median(Intensity_MC_Ch02),
    .groups = "drop"
  )


print(median.intensity)

# For the actual comparison I used the statistics computed based on the gating in IDEAS software

#Area median from software MO1 (brightfield)


Area_dat<- data.frame(Ab = c(rep("IgG",3),rep("RLRb",3)), value = c(66.44,64,62.78,78.44,76.06,76.33))
t.test(data = Area_dat, value~Ab)

#p= 001134

library('Hmisc')

# Calculate mean and standard deviation for each group
summary_data <- Area_dat %>%
  group_by(Ab) %>%
  summarise(
    mean = mean(value),
    sd = sd(value)
  )

my_comparisons <- list( c("IgG", "RLRb"))


P<- ggplot(Area_dat, aes(x = Ab, y = value, color = Ab)) +
  geom_point(size = 4, position = position_jitter(width = 0.2), alpha = 0.8) +  # Show individual points
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "black", aes(ymin = ..y.., ymax = ..y..)) +  # Add mean line
  # Add standard deviation error bars
  geom_errorbar(data = summary_data, aes(y = mean, ymin = mean - sd, ymax = mean + sd), 
                width = 0.2, color = "black") +  theme_minimal() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("black", "red")) +
  stat_compare_means(
    comparisons = my_comparisons,
    label = "p.signif",
    method = "t.test", method.args = list(alternative = "two.sided"), label.y = 80, size = 8, # Increase size of the asterisks
    vjust = 0.5, # Adjust vertical alignment of labels
    tip.length = 0.05 # Adjust line tip length
  ) + theme(legend.position = "none")  + ylab("Area") + xlab("") +
  theme(
    plot.title = element_text(size = 16, face = "bold"),  # Title font size
    axis.title.x = element_text(size = 14),              # X-axis title font size
    axis.title.y = element_text(size = 14),              # Y-axis title font size
    axis.text.x = element_text(size = 12),               # X-axis text font size
    axis.text.y = element_text(size = 12)                # Y-axis text font size
  ) 

# Circularity median from software MO2

circularity_dat_M02<- data.frame(Ab = c(rep("IgG",3),rep("RLRb",3)), value = c(5.367,5.428,8.499,9.577,9.689,9.664))
t.test(data = circularity_dat_M02, value~Ab)

# p=0.08972

summary_data <- circularity_dat_M02 %>%
  group_by(Ab) %>%
  summarise(
    mean = mean(value),
    sd = sd(value)
  )


P1<- ggplot(circularity_dat_M02, aes(x = Ab, y = value, color = Ab)) +
  geom_point(size = 4, position = position_jitter(width = 0.2), alpha = 0.8) +  # Show individual points
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "black", aes(ymin = ..y.., ymax = ..y..)) +  # Add mean line
  # Add standard deviation error bars
  geom_errorbar(data = summary_data, aes(y = mean, ymin = mean - sd, ymax = mean + sd), 
                width = 0.2, color = "black") +  theme_minimal() +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("black", "red")) +
  stat_compare_means(
    comparisons = my_comparisons,
    label = "p.signif",
    method = "t.test", method.args = list(alternative = "two.sided") , label.y = 11, size = 8, # Increase size of the asterisks
    vjust = 0.5, # Adjust vertical alignment of labels
    tip.length = 0.05 # Adjust line tip length
  ) + theme(legend.position = "none")  + ylab("Circularity index") + xlab("") +
  theme(
    plot.title = element_text(size = 16, face = "bold"),  # Title font size
    axis.title.x = element_text(size = 14),              # X-axis title font size
    axis.title.y = element_text(size = 14),              # Y-axis title font size
    axis.text.x = element_text(size = 12),               # X-axis text font size
    axis.text.y = element_text(size = 12)                # Y-axis text font size
  ) 

# Granularity (intensity ch06)

Intensity_dat_ch06<- data.frame(Ab = c(rep("IgG",3),rep("RLRb",3)), value = c(68183.09,42904.67,39866.03,110784.7,112058.12,82919.7))
t.test(data = Intensity_dat_ch06, value~Ab)

summary_data <- Intensity_dat_ch06 %>%
  group_by(Ab) %>%
  summarise(
    mean = mean(value),
    sd = sd(value)
  )

P2<- ggplot(Intensity_dat_ch06, aes(x = Ab, y = value, color = Ab)) +
  geom_point(size = 4, position = position_jitter(width = 0.2), alpha = 0.8) +  # Show individual points
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "black", aes(ymin = ..y.., ymax = ..y..)) +  # Add mean line
  # Add standard deviation error bars
  geom_errorbar(data = summary_data, aes(y = mean, ymin = mean - sd, ymax = mean + sd), 
                width = 0.2, color = "black") +  theme_minimal() +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("black", "red")) +
  stat_compare_means(
    comparisons = my_comparisons,
    label = "p.signif",
    method = "t.test", method.args = list(alternative = "two.sided"), label.y = 120000, size = 8, # Increase size of the asterisks
    vjust = 0.5, # Adjust vertical alignment of labels
    tip.length = 0.05 # Adjust line tip length
  ) + theme(legend.position = "none")  + ylab("Granularity") + xlab("") +
  theme(
    plot.title = element_text(size = 16, face = "bold"),  # Title font size
    axis.title.x = element_text(size = 14),              # X-axis title font size
    axis.title.y = element_text(size = 14),              # Y-axis title font size
    axis.text.x = element_text(size = 12),               # X-axis text font size
    axis.text.y = element_text(size = 12)                # Y-axis text font size
  ) 

#Intensity ch02 

Intensity_dat_ch02<- data.frame(Ab = c(rep("IgG",3),rep("RLRb",3)), value = c(19467.65,20586.06,21764.43,70420.64,71974.44,104731.19))
t.test(data = Intensity_dat_ch02, value~Ab)

summary_data <- Intensity_dat_ch02 %>%
  group_by(Ab) %>%
  summarise(
    mean = mean(value),
    sd = sd(value)
  )

P3<- ggplot(Intensity_dat_ch02, aes(x = Ab, y = value, color = Ab)) +
  geom_point(size = 4, position = position_jitter(width = 0.2), alpha = 0.8) +  # Show individual points
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "black", aes(ymin = ..y.., ymax = ..y..)) +  # Add mean line
  # Add standard deviation error bars
  geom_errorbar(data = summary_data, aes(y = mean, ymin = mean - sd, ymax = mean + sd), 
                width = 0.2, color = "black") +  theme_minimal() +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("black", "red")) +
  stat_compare_means(
    comparisons = my_comparisons,
    label = "p.signif",
    method = "t.test", method.args = list(alternative = "two.sided"), label.y = 110000, size = 8, # Increase size of the asterisks
    vjust = 0.5, # Adjust vertical alignment of labels
    tip.length = 0.05 # Adjust line tip length
  ) + theme(legend.position = "none")  + ylab("Alexa fluor 488") + xlab("") +
  theme(
    plot.title = element_text(size = 16, face = "bold"),  # Title font size
    axis.title.x = element_text(size = 14),              # X-axis title font size
    axis.title.y = element_text(size = 14),              # Y-axis title font size
    axis.text.x = element_text(size = 12),               # X-axis text font size
    axis.text.y = element_text(size = 12)                # Y-axis text font size
  ) 


library(patchwork)

# Combined plot 
# Combine the plots into a 2x2 grid with labels
combined_plot <- (P + P1) / (P2 + P3) +
  plot_annotation(tag_levels = 'A') # Automatically label plots with A, B, C, D

# Combine the plots into a single row with labels
combined_plot <- P + P1 + P2 + P3 +
  plot_layout(ncol = 4) +          # Arrange in 1 row (4 columns)
  plot_annotation(tag_levels = 'A') # Automatically label plots with A, B, C, D

print(combined_plot)

sessionInfo()

