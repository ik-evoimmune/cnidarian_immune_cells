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
  "~/immune_cells/cnidarian_immune_cells/08_phagocytosis/input/phagocytosis_results.csv",
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

sessionInfo()
