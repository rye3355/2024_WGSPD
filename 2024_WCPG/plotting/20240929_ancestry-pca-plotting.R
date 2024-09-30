library(ggplot2)
library(data.table)
library(tidyr)
library(dplyr)
library(RColorBrewer)


meta <- fread("/Users/rye/Projects/2024_WGSPD/subsetting/files/gnomad_v3.1_subset-metadata.tsv")
manifest <- fread("/Users/rye/Projects/2024_WGSPD/Analysis/QC/20240905_WGSPD_final-qcd-manifest.tsv")
all(manifest$s %in% meta$s)

# Pull over PC info
merged <- merge(manifest, meta[, c("s", "population_inference.pca_scores")],
                by = "s", 
                all.x = T, all.y = F)

# Convert PC info to fields
merged <- merged %>%
  separate(population_inference.pca_scores, into = paste0("PC", c(1:16)), sep = ",", convert = T) %>%
  mutate(across(everything(), ~gsub("\\[|\\]", "", .)))
merged <- merged %>%
  mutate(across(starts_with("PC"), as.numeric))


# Subset to wanted populations and reformat others
merged <- merged[!(merged$POP %in% c("ami", "mid", "oth")),]
merged %>%
  count(POP) %>%
  mutate(percentage = n / sum(n) * 100)
#   POP     n percentage
# 1 afr 11632 40.8498683
# 2 amr  2279  8.0035119
# 3 asj   377  1.3239684
# 4 eas    40  0.1404741
# 5 fin  4785 16.8042142
# 6 nfe  9329 32.7620720
# 7 sas    33  0.1158911

merged$"Genetic Ancestry" <- case_when(merged$POP == "afr" ~ "African American (39.88%)",
                                       merged$POP == "amr" ~ "Admixed (8.40%)",
                                       merged$POP == "asj" ~ "Ashkenazi Jewish (1.35%)",
                                       merged$POP %in% c("eas", "sas") ~ "Asian (0.20%)",
                                       merged$POP == "fin" ~ "Finnish (14.51%)",
                                       merged$POP == "nfe" ~ "Non-Finnish European (33.20%)")
merged$`Genetic Ancestry` <- factor(merged$`Genetic Ancestry`, 
                                    levels = c("African American (39.88%)", "Non-Finnish European (33.20%)",
                                               "Finnish (14.51%)", "Admixed (8.40%)", "Ashkenazi Jewish (1.35%)",
                                               "Asian (0.20%)"))

colors <- c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", "#b15928")


pca_plot <- ggplot(merged, aes(x = PC1, y = PC2, color = `Genetic Ancestry`)) +
  geom_point(size = 1.5) +
  scale_color_manual(values = colors) +  # Apply professional color palette
  labs(
    x = "Principal Component 1",
    y = "Principal Component 2",
    color = "Genetic Ancestry"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, margin = margin(t = 5, r = 5, b = 5, l = 5)),  # Adjust axis titles' margins
    axis.text = element_text(size = 10, margin = margin(t = 0, r = 0, b = 0, l = 0)),  # Remove extra margin around tick labels
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "bottom",  # Move legend to the bottom for better presentation
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add a black border around the plot
    plot.margin = margin(10, 10, 10, 10)  # Adjust the plot margin for overall balance
  )
pca_plot

ggsave("outputs/ancestry-pca-plots/20240929_ancestry-pca_legend-bottom.jpg", plot = pca_plot, width = 9, height = 9)

pca_plot <- ggplot(merged, aes(x = PC1, y = PC2, color = `Genetic Ancestry`)) +
  geom_point(size = 2) +
  scale_color_manual(values = colors) +  # Apply professional color palette
  labs(
    x = "Principal Component 1",
    y = "Principal Component 2",
    color = "Genetic Ancestry"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, margin = margin(t = 5, r = 5, b = 5, l = 5)),  # Adjust axis titles' margins
    axis.text = element_text(size = 10, margin = margin(t = 0, r = 0, b = 0, l = 0)),  # Remove extra margin around tick labels
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "right",  # Move legend to the bottom for better presentation
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add a black border around the plot
    plot.margin = margin(10, 10, 10, 10)  # Adjust the plot margin for overall balance
  )
ggsave("outputs/ancestry-pca-plots/20240929_ancestry-pca_legend-right.jpg", plot = pca_plot, width = 11, height = 8)
