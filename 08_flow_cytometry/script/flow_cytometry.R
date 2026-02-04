# ============================================================
# Phagocytosis analysis + plots
# ============================================================

# ---- Packages ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
})


# ---- Read input ----

df1 <- read.csv(
  "~/immune_cells/cnidarian_immune_cells/08_flow_cytometry/input/phagocytosis_results.csv",
  header = TRUE
)

# ============================================================
# Part 1: Untreated (ANOVA + Tukey + barplot)
# ============================================================

df_untreated <- df1 %>% filter(Assay == "Untreated")

anova_result <- aov(Percent ~ Fraction, data = df_untreated)
summary(anova_result)

posthoc <- TukeyHSD(anova_result)
print(posthoc)

pairwise_results <- as.data.frame(posthoc$Fraction) %>%
  mutate(Significance = cut(
    `p adj`,
    breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
    labels = c("***", "**", "*", "ns")
  ))

comparisons_untreated <- list(
  c("RLRb::mCherry NaCl", "RLRb::mCherry pIC"),
  c("RLRb::mCherry NaCl", "RLRb::mCherry uninjected"),
  c("RLRb::mCherry NaCl", "WT"),
  c("RLRb::mCherry pIC", "RLRb::mCherry uninjected"),
  c("RLRb::mCherry pIC", "WT"),
  c("RLRb::mCherry uninjected", "WT")
)

p_untreated <- ggbarplot(
  df_untreated,
  x = "Fraction",
  y = "Percent",
  add = c("mean_sd", "jitter"),
  add.params = list(shape = "Fraction"),
  fill = "Fraction",
  palette = c("black", "grey", "yellow", "red"),
  position = position_dodge(0.8),
  ylab = "% of cells",
  xlab = "",
  legend.title = ""
) +
  stat_compare_means(comparisons = comparisons_untreated,
                     label = "p.signif",
                     method = "t.test") +
  theme_minimal() +
  theme(legend.position = "none")

p_untreated


# ============================================================
# Part 2: Assay-specific scatter + mean ± SD + t-test
#   - E.coli, Ovalbumin, S.aureus
# ============================================================

my_comparisons <- list(c("mCherry_negative", "mCherry_positive"))

make_assay_plot <- function(df_all, assay, label_y) {
  df_assay <- df_all %>% filter(Assay == assay)
  
  summary_data <- df_assay %>%
    group_by(Assay, Fraction) %>%
    summarise(
      mean = mean(Percent, na.rm = TRUE),
      sd   = sd(Percent, na.rm = TRUE),
      .groups = "drop"
    )
  
  ggplot(df_assay, aes(x = Fraction, y = Percent, color = Fraction)) +
    geom_point(size = 4,
               position = position_jitter(width = 0.2),
               alpha = 0.8) +
    stat_summary(
      fun = mean,
      geom = "crossbar",
      width = 0.5,
      color = "black",
      aes(ymin = after_stat(y), ymax = after_stat(y))
    ) +
    geom_errorbar(
      data = summary_data,
      aes(
        y = mean,
        ymin = mean - sd,
        ymax = mean + sd
      ),
      width = 0.2,
      color = "black"
    ) +
    scale_color_manual(values = c("black", "red")) +
    stat_compare_means(
      comparisons = my_comparisons,
      label = "p.signif",
      method = "t.test",
      method.args = list(alternative = "two.sided"),
      label.y = label_y,
      size = 8,
      vjust = 0.5,
      tip.length = 0.05
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 16, face = "bold"),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x  = element_text(size = 12),
      axis.text.y  = element_text(size = 12)
    ) +
    ylab("Percent") +
    xlab("") +
    ggtitle(assay)
}

P1 <- make_assay_plot(df1, assay = "E.coli", label_y = 45)
P2 <- make_assay_plot(df1, assay = "Ovalbumin", label_y = 60)
P3 <- make_assay_plot(df1, assay = "S.aureus", label_y = 60)

P1
P2
P3


# Apoptosis (Apotracker green) assay --------------------------------------

library(tidyverse)
library(ggpubr)

# ---------- Generate a table with the results ----------
df <- tribble(
  ~`mCherry+`, ~`mCherry-`, ~type,
  11.2, 48.6, "percent",
  11.0, 47.8, "percent",
  9.52, 46.5, "percent",
  6.99, 38.3, "percent",
  26.8, 111.0, "MFI",
  26.8, 92.6,  "MFI",
  24.7, 73.4,  "MFI",
  25.9, 35.6,  "MFI"
) %>% group_by(type) %>% mutate(replicate = row_number()) %>% ungroup()

# ---------- Plotting function ----------
paired_plot <- function(dat, ylab, ylim = NULL,
                                col_plus = "#5BA4A4",   
                                col_minus = "#D9A441",  
                                save_pdf = NULL, width_mm = 55, height_mm = 75) {
  
  # long format (mCherry+ left, mCherry- right)
  long <- dat %>%
    pivot_longer(c(`mCherry+`, `mCherry-`),
                 names_to = "group", values_to = "value") %>%
    mutate(group = factor(group, levels = c("mCherry+", "mCherry-")))
  
  # paired t-test
  wide <- dat %>% select(replicate, `mCherry+`, `mCherry-`)
  tt   <- t.test(wide$`mCherry+`, wide$`mCherry-`, paired = TRUE)
  pval <- tt$p.value
  star <- if      (pval < 0.001) "***"
  else if (pval < 0.01)  "**"
  else if (pval < 0.05)  "*"
  else                   "n.s."
  
  # bracket position
  y_max <- if (is.null(ylim)) max(long$value) else max(ylim)
  p_ann <- tibble(group1 = "mCherry+", group2 = "mCherry-",
                  y.position = y_max * 1.03, p = pval, label = star)
  
  # plot 
  g <- ggplot(long, aes(group, value)) +
    geom_line(aes(group = replicate),
              color = "black", linewidth = 1.0, alpha = 0.85) +
    geom_point(aes(color = group),
               size = 3.8, stroke = 0.9, shape = 21, fill = "white") +
    scale_color_manual(values = c("mCherry+" = col_plus, "mCherry-" = col_minus), guide = "none") +
    labs(x = NULL, y = ylab) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x = element_text(face = "bold"),
      plot.margin = margin(3, 3, 3, 3)
    ) +
    coord_cartesian(ylim = ylim, clip = "off") +
    stat_pvalue_manual(p_ann, xmin = "group1", xmax = "group2",
                       label = "label", tip.length = 0.01, size = 4)
  
  if (!is.null(save_pdf)) {
    ggsave(save_pdf, g, device = cairo_pdf,
           width = width_mm/25.4, height = height_mm/25.4, units = "in")
  }
  
  list(plot = g, ttest = tt)
}

# View and print exact p-values
p_percent$plot; p_percent$ttest
p_mfi$plot;     p_mfi$ttest



# Proliferation 

# EdU incorporation assay -------------------------------------------------

## ===============================
##  EdU+ fraction: RLRb+ vs RLRb-
## ===============================

library(tidyverse)
library(ggpubr)

## ---- 1. Input data: Edu- and Edu+ counts ----

# RLRb+ cells
plus_edu_minus <- c(32.8, 34.5, 42.0, 40.4)  # Edu-
plus_edu_plus  <- c(12.9, 14.3, 15.4, 14.0)  # Edu+

# RLRb- cells
minus_edu_minus <- c(45.4, 41.3, 33.8, 37.7)  # Edu-
minus_edu_plus  <- c(8.94, 9.87, 8.82, 7.87)  # Edu+

## ---- 2. Compute EdU positive % out of total----

perc_plus  <- plus_edu_plus  / (plus_edu_minus  + plus_edu_plus)  * 100
perc_minus <- minus_edu_plus / (minus_edu_minus + minus_edu_plus) * 100

# put into a tidy data frame
df <- tibble(
  replicate = 1:4,
  `RLRb+` = perc_plus,
  `RLRb-` = perc_minus
)

df  

## ---- 3. Paired t-test (RLRb+ vs RLRb-) ----

tt <- t.test(df$`RLRb+`, df$`RLRb-`, paired = TRUE)
tt  

pval <- tt$p.value
star <- if      (pval < 0.001) "***" else
  if      (pval < 0.01)  "**"  else
    if      (pval < 0.05)  "*"   else "n.s."

## ---- 4. Long format for plotting ----
long <- df %>%
  pivot_longer(cols = c(`RLRb-`, `RLRb+`),
               names_to = "group",
               values_to = "value") %>%
  mutate(group = factor(group, levels = c("RLRb+", "RLRb-")))  # RLRb+ on the LEFT

# tighter y-limits based on data
y_min <- min(long$value)
y_max <- max(long$value)
y_range <- y_max - y_min
ylim_use <- c(y_min - 0.15 * y_range,   
              y_max + 0.25 * y_range)   

p_ann <- tibble(
  group1 = "RLRb+",
  group2 = "RLRb-",
  y.position = y_max + 0.18 * y_range,  
  label = star
)

## ---- 5. Compact paired plot ----
col_plus  <- "#5BA4A4"  # teal
col_minus <- "#D9A441"  # mustard

g <- ggplot(long, aes(group, value)) +
  geom_line(aes(group = replicate),
            color = "black", linewidth = 1.0, alpha = 0.85) +
  geom_point(aes(color = group),
             size = 3.6, shape = 21, fill = "white", stroke = 1) +
  scale_color_manual(values = c("RLRb+" = col_plus, "RLRb-" = col_minus),
                     guide = "none") +
  labs(x = NULL, y = "EdU⁺ cells (% of total)") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x  = element_text(face = "bold"),
    axis.title.y = element_text(margin = margin(r = 4)),
    plot.margin  = margin(3, 3, 3, 3)
  ) +
  coord_cartesian(ylim = ylim_use, clip = "off") +
  stat_pvalue_manual(
    p_ann,
    xmin = "group1", xmax = "group2",
    label = "label",
    size = 4,
    tip.length = 0.01
  )

g


# DRAQ5 -------------------------------------------------------------------


## ---- Input data ----
draq <- c(89.6, 83.6, 88.8, 82.9, 83.3, 71.1)

df <- tibble(
  condition = "DRAQ5+",
  value = draq
)

## ---- Plot ----
g_draq <- ggplot(df, aes(x = condition, y = value)) +
  # raw points
  geom_point(size = 3, shape = 21, fill = "white", stroke = 1) +
  # mean ± SEM
  stat_summary(fun = mean, geom = "point",
               shape = 23, size = 3.6, fill = "grey90") +
  stat_summary(fun.data = mean_se, geom = "errorbar",
               width = 0.08, linewidth = 0.9) +
  labs(x = NULL, y = "DRAQ5⁺ cells (%)") +
  coord_cartesian(ylim = c(0, 100), clip = "off") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(face = "bold"),
    plot.margin  = margin(3, 3, 3, 3)
  )

g_draq



# Viability -------------------------------------------------------------

#––––– Packages –––––#
library(ggplot2)
library(ggpubr)
library(dplyr)

#––––– Data –––––#
df <- data.frame(
  viability = c(79.8, 81, 82.7, 90.4, 89.5, 83.9,
                81.3, 76.3, 88.2, 81.4, 84.5, 86),
  treatment = c("pIC", "NaCl", "pIC", "NaCl", "pIC", "NaCl",
                "pIC", "NaCl", "pIC", "NaCl", "pIC", "NaCl")
)

df$treatment <- factor(df$treatment, levels = c("NaCl", "pIC"))

#––––– Stats –––––#
t_test <- t.test(viability ~ treatment, data = df)
t_test

#––––– Boxplot –––––#
ggplot(df, aes(x = treatment, y = viability, fill = treatment)) +
  geom_boxplot(
    width = 0.45,
    alpha = 0.75,
    color = "black",
    linewidth = 0.6,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    size = 3.5,           
    shape = 21,
    stroke = 0.5,
    color = "black",
    aes(fill = treatment)
  ) +
  stat_compare_means(
    method = "t.test",
    label = "p.format",
    size = 4.5,
    label.y = 96
  ) +
  scale_fill_manual(
    values = c("NaCl" = "#999999",     
               "pIC"  = "#0072B2")     
  ) +
  scale_y_continuous(
    name = "Viability (%)",
    limits = c(70, 97),
    breaks = seq(70, 95, 5),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(x = "") +
  theme_classic(base_size = 13) +
  theme(
    axis.text  = element_text(color = "black"),
    axis.title.y = element_text(margin = margin(r = 8)),
    axis.line  = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    legend.position = "none"
  )



# RLRb::mCherry reporter line  --------------------------------------------
# Create the data frame
# Original data
data <- data.frame(
  treatment = c("NaCl", "NaCl", "NaCl", "NaCl", "pIC", "pIC", "pIC", "pIC", "WT", "WT", "WT", "WT"),
  cells = c(8.15, 9.43, 9.20, 5.16, 28.00, 35.50, 32.30, 24.70, 1.07, 0.49, 2.16, 0.37)
)

# Reorder the treatment levels
data$treatment <- factor(data$treatment, levels = c("WT", "NaCl", "pIC"))

# Plot with ggpubr
plot <- ggbarplot(
  data, 
  x = "treatment", 
  y = "cells", 
  add = c("mean_sd", "jitter"), 
  fill = "gray70", 
  error.plot = "errorbar", 
  width = 0.7,
  position = position_dodge(0.8)
) +
  stat_compare_means(
    method = "t.test",
    comparisons = list(c("WT", "NaCl"), c("WT", "pIC"), c("NaCl", "pIC")),
    label = "p.signif",
    exact = FALSE
  ) +
  labs(
    title = "",
    x = "Treatment",
    y = "% of mCherry positive cells"
  ) +
  theme_minimal(base_size = 16) + 
  theme(
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 14),  
    legend.position = "none",
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold")
  )

# Print the plot
print(plot)


# Granularity 
data <- data.frame(
  mCherry = c("mCherry negative", "mCherry negative", "mCherry negative", "mCherry positive", "mCherry positive", "mCherry positive"),
  Percent = c(39.9,57.8,19.0, 84.0,89.7,68.8)
)

# Reorder the treatment levels
data$mCherry <- factor(data$mCherry, levels = c("mCherry negative", "mCherry positive"))

# Plot with ggpubr
plot <- ggbarplot(
  data, 
  x = "mCherry", 
  y = "Percent", 
  add = c("mean_sd", "jitter"), 
  fill = "gray70", 
  error.plot = "errorbar", 
  width = 0.7,
  position = position_dodge(0.8)
) +
  stat_compare_means(
    method = "t.test",
    comparisons = list(c("mCherry negative", "mCherry positive")),
    label = "p.signif",
    exact = FALSE
  ) +
  labs(
    title = "",
    x = "Treatment",
    y = "% of large granular cells"
  ) +
  theme_minimal(base_size = 16) + 
  theme(
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 14),  
    legend.position = "none",
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold") 
  )

print(plot)

sessionInfo()
