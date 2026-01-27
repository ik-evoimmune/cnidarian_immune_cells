---
title: "Mass-spec analysis"
author: "Itamar Kozlovski"
date: "2026-01-20"
---

```         
# Load libraries
library(readr)
library(ggplot2)
library(dplyr)

mass_spec<- readxl::read_excel("~/immune_cells/cnidarian_immune_cells/06_Mass_spec_analysis/input_files/MaxQuant-Perseus-stats-Rubi-20210201.xlsx", sheet = 2)

## New names:
## • `` -> `...10`

glimpse(mass_spec)

## Rows: 567
## Columns: 12
## $ `p-value`                             [3m[38;5;246m<dbl>[39m[23m 0.001412908…
## $ `ratio of LFQ intensity: RLRb vs IgG` [3m[38;5;246m<dbl>[39m[23m 2.002604e+0…
## $ `Fasta headers`                       [3m[38;5;246m<chr>[39m[23m "NVE16598",…
## $ `LFQ intensity alpha-RLRb-1`          [3m[38;5;246m<dbl>[39m[23m 35.97522, 2…
## $ `LFQ intensity alpha-RLRb-2`          [3m[38;5;246m<dbl>[39m[23m 30.99920, 2…
## $ `LFQ intensity alpha-RLRb-3`          [3m[38;5;246m<dbl>[39m[23m 36.37318, 2…
## $ `LFQ intensity Rabbit-IgG-1`          [3m[38;5;246m<dbl>[39m[23m 19.44059, 2…
## $ `LFQ intensity Rabbit-IgG-2`          [3m[38;5;246m<dbl>[39m[23m 19.77698, 2…
## $ `LFQ intensity Rabbit-IgG-3`          [3m[38;5;246m<dbl>[39m[23m 21.26127, 2…
## $ ...10                                 [3m[38;5;246m<lgl>[39m[23m NA, NA, NA,…
## $ `-LOG(P-value)`                       [3m[38;5;246m<dbl>[39m[23m 2.849886, 2…
## $ `log Difference`                      [3m[38;5;246m<dbl>[39m[23m 14.2895896,…

#Read dictionary file
dictionary<- read.csv("~/immune_cells/cnidarian_immune_cells/06_Mass_spec_analysis/input_files/Dictionary_Kozlovski_et.al.csv")

# Rank by p-value
df <- mass_spec %>%
  arrange(desc(`log Difference`)) %>%
  mutate(rank = row_number())

df$`Fasta headers`[1] <- "RLRb"
df$`Fasta headers`[2] <- "TRIP11"
df$`Fasta headers`[3] <- "ALDH1A2"
df$`Fasta headers`[4] <- "LRP1B"
df$`Fasta headers`[5] <- "ADGRL3"

top_hits <- df %>% filter(rank <= 5)


merged_df<- left_join(df,dictionary, by = c('Fasta headers' = 'NVE_gene_model'))

top_hits <- merged_df %>% filter(rank <= 5)


library(ggplot2)
library(ggrepel)

ggplot(df, aes(x = rank, y = `log Difference`)) +
  geom_point(size = 1, alpha = 0.7, color = "black") +
  geom_text_repel(
    data = top_hits,
    aes(label = `Fasta headers`),
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
```

![](mass_spec_analysis_files/figure-markdown_strict/unnamed-chunk-1-1.png)

```         
ratio <- top_hits$`ratio of LFQ intensity: RLRb vs IgG`[1]
log10(ratio)

## [1] 4.301595
```
