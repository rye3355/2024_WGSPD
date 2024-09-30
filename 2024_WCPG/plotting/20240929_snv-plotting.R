library(ggplot2)
library(data.table)
library(tidyr)
library(dplyr)
library(RColorBrewer)

ptv <- fread("inputs/20240929_ptv-ac5_old-analysis_rate-ratio.tsv")
syn <- fread("inputs/20240929_syn-ac5_old-analysis_rate-ratio.tsv")


ptv$Category[ptv$Category == "All genes (18356)"] <- "All genes"
syn$Category[syn$Category == "All genes (19369)"] <- "All genes"

ptv$type <- "PTV"
syn$type <- "SYN"

df <- rbind(ptv, syn)
df$Category <- factor(df$Category, 
                      levels = c("Constrained Genes (3570)", "Schizophrenia (32)",
                                 "Autism Spectrum Disorder (72)", "Neurodevelopmental Disorder (373)",
                                 "All genes"))
df$type <- factor(df$type, levels = c("PTV", "SYN"))



odds_ratio_plot <- ggplot(df, aes(x = Category, y = OR, color = Category)) +
  geom_point(size = 4) +  # Points for estimates
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2, size = 1) +  # Error bars
  labs(
    title = "Odds Ratios by Category",
    x = "Category",
    y = "Odds Ratio",
    color = "Type"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "gray80"),
    panel.grid.minor = element_blank()  # Optional: remove minor grid lines for clarity
  ) +
  facet_wrap(~ type)  # Facet by type
odds_ratio_plot

ggsave("outputs/ancestry-pca-plots/20240929_ancestry-pca_legend-right.jpg", plot = pca_plot, width = 11, height = 8)
