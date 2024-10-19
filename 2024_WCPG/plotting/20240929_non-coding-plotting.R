library(ggplot2)
library(data.table)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(scales) 

sv <- fread("inputs/non-coding/fisher_non-coding_table.tsv")
snv <- fread("inputs/non-coding/rate-ratio_SNV_non-coding_table.tsv")
colnames(snv) <- colnames(sv)

df <- rbind(sv, snv)
df$Category <- factor(df$Category, 
                      levels = rev(c("SNV",
                                     "DEL", "DUP")))



odds_ratio_plot <- ggplot(df, aes(x = OR, y = Category, color = Category)) +
  geom_point(size = 4.5) +  # Points for estimates
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), 
                 height = 0.3, size = 1, 
                 ) +  # Error bars
  labs(
    #title = "",
    #x = "",
    y = "Odds Ratio",
    color = "Gene set"
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#999999", size = 0.75) +
  scale_color_manual(breaks = c("SNV",
                                "DEL", "DUP"),
                     values = c("#999999",
                                "#C34631", "#3D74AD")) +
  scale_x_continuous(limits = c(0, 5), 
                     breaks = seq(0.5, 3, 0.5),
                     labels = label_number(accuracy = 0.01)) +
  coord_cartesian(xlim = c(0.5, 2.5)) +
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

ggsave("outputs/non-coding/20241004_non-coding_old-sv-new-snv.jpg", plot = odds_ratio_plot, width = 8, height = 4)
