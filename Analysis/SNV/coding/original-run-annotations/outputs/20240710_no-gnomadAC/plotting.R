# gs://fc-712bf694-df47-4018-8788-bfdc120cdd67/202404_GD-analysis/outputs/20240501_Finish_gd-case-con_pass-samples_fisher_BP.tsv
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)

# Setup
d <- fread("/Users/rye/Projects/2024_WGSPD/Analysis/SNV/coding/original-run-annotations/outputs/20240710_no-gnomadAC/20240710_BDSCZ_synonymous_AC5_case-control_CMH_individual.tsv")
d <- d[d$RES.test_statistic != Inf,]


# Convert string representation of array to array and sum
d$a <- str_replace_all(d$case_carriers, c("\\["="", "\\]"=""))
d$a <- lapply(strsplit(d$a, ","), as.integer)
d$case_carriers_total <- sapply(d$a, sum)

d$a <- str_replace_all(d$control_carriers, c("\\["="", "\\]"=""))
d$a <- lapply(strsplit(d$a, ","), as.integer)
d$control_carriers_total <- sapply(d$a, sum)

# Remove less than 10 carriers
d <- d[d$case_carriers_total + d$control_carriers_total >= 100, ]

# Remove non-autosomal and p = 1
d$p.value <- d$RES.p_value
d <- d[order(d$p.value, decreasing = F),]
d$o <- -log10(d$p.value)
d$e <- -log10(ppoints(d$p.value))
d$clower <- -log10(qbeta(p = (1 - 0.95) / 2, shape1 = 1:dim(d)[1], shape2 = dim(d)[1]:1))
d$cupper <- -log10(qbeta(p = (1 + 0.95) / 2, shape1 = 1:dim(d)[1], shape2 = dim(d)[1]:1))

d$estimate <- d$OR
d$col <- case_when(d$estimate < 1 ~ "OR < 1",
                   d$estimate >= 1 & d$estimate < 2 ~ "1 ≤ OR < 2",
                   d$estimate >= 2 ~ "OR ≥ 2")
d$col <- factor(d$col, levels = rev(c("OR < 1", "1 ≤ OR < 2","OR ≥ 2")))
d$GD_clean <- d$gene_symbol

p <- ggplot(d, aes(x = e, y = o, label = GD_clean)) +
  geom_point(aes(color = col), size = 2) +
  scale_x_continuous(limits = c(0-0.05, max(d$e) + 0.05), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0-0.25, max(d$o) + 2), expand = c(0, 0)) +
  scale_color_manual(name = "", values=c("#CD323D", "#32CDC2", "black")) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  geom_ribbon(aes(ymin=clower, ymax=cupper), alpha=0.2) +
  theme_bw() +
  labs(x = expression(paste("Expected -log"[10], plain(P))), y = expression(paste("Observed -log"[10], plain(P))), title = "Synonymous Singletons (Gene carrier count >= 100, 2557)") +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    panel.border = element_blank(),
    axis.title.x = element_text(hjust=0.5), 
    axis.title.y = element_text(hjust=0.5),
    plot.title = element_text(hjust = 0.5),
    axis.line = element_line(colour = "black")) +
  geom_label_repel(data = d[c(1:5), ], 
                   nudge_x = -0.25, nudge_y = 0.25,
                   size=4, segment.size=0.25)
  # 
  # geom_label_repel(data = d[c(2), ], 
  #                  nudge_x = -0.25, nudge_y = 0.25,
  #                  size=4, segment.size=0.25) +
  # geom_label_repel(data = d[c(3:4), ], 
  #                  nudge_x = 0.25, nudge_y = -0.25,
  #                  size=4, segment.size=0.25) + 
  # geom_label_repel(data = d[c(5), ], 
  #                  nudge_x = 0.3, nudge_y = -0.3,
  #                  size=4, segment.size=0.25) + 
  # geom_label_repel(data = d[c(6:10), ], 
  #                  nudge_x = -0.4, nudge_y = 0.3,
  #                  size=4, segment.size=0.25)  
p


# library(ggridges)
# ggplot(d, aes(x=p.value)) + 
#   geom_histogram(color="black", fill="white")


