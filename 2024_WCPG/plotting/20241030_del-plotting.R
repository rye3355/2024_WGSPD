library(ggplot2)
library(data.table)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(scales) 

del <- fread("inputs/SV/del/fisher_del-ac5_table.tsv")
colnames(del) <- c("Gene_set", "Category", "OR", "CI_low", "CI_high", "pval", 
                   "rarevar_case", "rarevar_control", "novar_case", "novar_control",
                   "length_thresh", "corrected", "both_corrected")
del$Gene_set <- case_when(del$Gene_set == "ASD" ~ "Autism Spectrum Disorder (72)",
                          del$Gene_set == "NDD" ~ "Neurodevelopmental Disorder (373)",
                          del$Gene_set == "pLI-constrained" ~ "Constrained Genes (3570)",
                          del$Gene_set == "SCHEMA" ~ "Schizophrenia (32)",)
del$Category <- case_when(del$Category == "PREDICTED_LOF" ~ "Exonic",
                          del$Category == "PREDICTED_INTRONIC" ~ "Intronic",
                          del$Category == "PREDICTED_PROMOTER" ~ "Promoter",
                          del$Category == "PREDICTED_UTR" ~ "5'UTR",)
df <- del
df$Gene_set <- factor(df$Gene_set, 
                      levels = rev(c("Schizophrenia (32)",
                                     "Autism Spectrum Disorder (72)", "Neurodevelopmental Disorder (373)",
                                     "Constrained Genes (3570)")))
df$Category <- factor(df$Category, 
                      levels = c("Exonic", "Intronic", 
                                 "Promoter","5'UTR"))


odds_ratio_plot <- ggplot(df, aes(x = OR, y = Category, color = Gene_set)) +
  geom_point(size = 4.5, position=position_dodge(width = 200)) +  # Points for estimates
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), 
                 height = 0.3, size = 1, 
                 position=position_dodge(width = 200)) +  # Error bars
  facet_wrap(~Category, ncol = 2, ) +
  labs(
    #title = "",
    #x = "",
    y = "Odds Ratio",
    color = "Gene set"
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#999999", size = 0.75) +
  scale_color_manual(breaks = c("Schizophrenia (32)", 
                                "Autism Spectrum Disorder (72)", "Neurodevelopmental Disorder (373)",
                                "Constrained Genes (3570)"),
                     values = c("#0072B2", 
                                "#009E73", "#D55E00", 
                                "#CC79A7")) +
  scale_x_continuous(limits = c(0, 300), 
                     breaks = seq(0, 8, 1),
                     labels = label_number(accuracy = 1)) +
  coord_cartesian(xlim = c(0.5, 7)) +
  theme_classic(base_size = 18) +
  theme(
    legend.position = "bottom", legend.text.align = 0,
    plot.title = element_blank(),
    axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(),
    axis.text.x = element_text(size = 20),
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
    panel.grid.major = element_line(color = "grey80", size = 0.5),  
    panel.grid.minor = element_blank(),
    text = element_text(family = "Arial")
  ) +
  guides(color = guide_legend(override.aes = list(shape = 15, linetype = "blank"), nrow = 2)) 
odds_ratio_plot

ggsave("outputs/SV-plots/20241030_sv-plotting_del-ac5_old-analysis.jpg", plot = odds_ratio_plot, width = 8, height = 8)






