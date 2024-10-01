library(ggplot2)
library(data.table)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(scales) 

ptv <- fread("inputs/SNV/20240929_ptv-ac5_old-analysis_rate-ratio.tsv")
syn <- fread("inputs/SNV/20240929_syn-ac5_old-analysis_rate-ratio.tsv")


ptv$Category[ptv$Category == "All genes (18356)"] <- "All genes"
syn$Category[syn$Category == "All genes (19369)"] <- "All genes"

ptv$type <- "PTV"
syn$type <- "SYN"

df <- rbind(ptv, syn)
df$Category <- factor(df$Category, 
                      levels = rev(c("Schizophrenia (32)",
                                 "Autism Spectrum Disorder (72)", "Neurodevelopmental Disorder (373)",
                                 "Constrained Genes (3570)", "All genes")))
df$type <- factor(df$type, levels = c("PTV", "SYN"))



odds_ratio_plot <- ggplot(df, aes(x = OR, y = Category, color = Category)) +
  geom_point(size = 4.5, position=position_dodge(width = 3)) +  # Points for estimates
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), 
                 height = 0.3, size = 1, 
                 position=position_dodge(width = 3)) +  # Error bars
  facet_wrap(~type, ncol = 1) +
  labs(
    #title = "",
    #x = "",
    y = "Odds Ratio",
    color = "Gene set"
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#999999", size = 0.75) +
  scale_color_manual(breaks = c("Schizophrenia (32)", 
                                "Autism Spectrum Disorder (72)", "Neurodevelopmental Disorder (373)",
                                "Constrained Genes (3570)", "All genes"),
                     values = c("#0072B2", 
                                "#009E73", "#D55E00", 
                                "#CC79A7", "#999999")) +
  scale_x_continuous(limits = c(0.85, 2), 
                     breaks = seq(1, 2, 0.2),
                     labels = label_number(accuracy = 0.01)) +
  theme_classic(base_size = 18) +
  theme(
    legend.position = "right", legend.text.align = 0,
    strip.text = element_blank(),
    plot.title = element_blank(),
    axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(),
    axis.text.x = element_text(size = 20),
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
    panel.grid.major = element_line(color = "grey80", size = 0.5),  
    panel.grid.minor = element_blank(),
    text = element_text(family = "Arial")
  ) +
  guides(color = guide_legend(override.aes = list(shape = 15, linetype = "blank"))) 
odds_ratio_plot

ggsave("outputs/SNV-plots/20240929_snv-plotting_ptv-syn-ac5_old-analysis.jpg", plot = odds_ratio_plot, width = 11, height = 6.5)
