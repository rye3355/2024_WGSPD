# /Users/rye/Projects/WGSPD/gnomad-v3-subset/VEP/BDSCZvsCONT/06_BDSCZvsCONT_VEP_association_ptv-ac5_double-counting.R
library(tidyverse)
library(tidyr)
library(dplyr)
library(broom)
library(data.table)
library(ggrepel)
library(gridExtra)
library(fmsb)

setwd("/Users/rye/Projects/WGSPD/gnomad-v3-subset/VEP/BDSCZvsCONT/")
df = data.frame(matrix(vector(), 0, 9),
                stringsAsFactors=F)
df_rr = data.frame(matrix(vector(), 0, 9),
                   stringsAsFactors=F)



manifest_file <- "../../variant-qc/BDSCZvsCONT/DS-MANIFEST-WGSPD-WGS_BDSCZvsCONT_FIXED-PHENOTYPES_POPULATION.tsv"
manifest <- fread(manifest_file, sep = '\t', header=T)

sample_to_keep_file <- "../../variant-qc/BDSCZvsCONT/12_BDSCZvsCONT_ancestry-4-MAD_passing-samples.tsv"
samples_to_keep <- fread(sample_to_keep_file, sep = '\t', header = T)
table(samples_to_keep$AFFECTED_STATUS_FIXED, samples_to_keep$Population, useNA = "always")
table(samples_to_keep$AFFECTED_STATUS_FIXED, useNA = "always")
#    CASE CONTROL    <NA> 
#.  7749    4025       0 
table(samples_to_keep$Population, useNA = "always")
#  afr  amr  fin  nfe  sas <NA> 
#  7119  314 1746 2593    2    0 

dim(samples_to_keep) #11774    2

# Read in large counts table
gene_count <- fread("20230522_BDSCZvsCONT_pLoF_ac5.tsv", header=T, sep = "\t")
gene_count <- subset(gene_count, select = c("gene_symbol", samples_to_keep$s))
rownames(gene_count) <- gene_count$gene_symbol

# GnomAD LEOUF
## Extract gnomad genes
gnomad <- fread("gnomad.v2.1.1.lof_metrics.by_gene.txt.gz", data.table = F,h=T)
gnomad[gnomad$gene == "DOPEY1",]$gene <- "DOP1A"
gnomad_genes <- gnomad[gnomad$pLI>.9,]$gene


## Subset gene_count down to gnomad_genes
gnomad_genes_count <- as.data.frame(t(gene_count[gene_count$"gene_symbol" %in% gnomad_genes, ]))
gnomad_genes_count <- gnomad_genes_count[-1, ]
gnomad_genes_count <- mutate_all(gnomad_genes_count, function(x) as.numeric(as.character(x)))

## Split into case/control
gnomad_genes_count_case <- gnomad_genes_count[samples_to_keep[samples_to_keep$AFFECTED_STATUS_FIXED == "CASE"]$s, ]
gnomad_genes_count_control <- gnomad_genes_count[samples_to_keep[samples_to_keep$AFFECTED_STATUS_FIXED == "CONTROL"]$s, ]

gnomad_genes_count_case <- gnomad_genes_count_case %>% mutate(Total = select(., V1:V2027) %>% rowSums(na.rm = TRUE))
gnomad_genes_count_control <- gnomad_genes_count_control %>% mutate(Total = select(., V1:V2027) %>% rowSums(na.rm = TRUE))

gnomad_df <- data.frame("rarevar" = c(sum(gnomad_genes_count_case$Total > 0),  sum(gnomad_genes_count_control$Total > 0)),
                        "novar" =   c(sum(gnomad_genes_count_case$Total == 0), sum(gnomad_genes_count_control$Total == 0)),
                        row.names = c("case", "control"))
# rarevar novar
# case       2591  5158
# control    1184  2841

a1 <- fisher.test(gnomad_df)
a1
# 	Fisher's Exact Test for Count Data
# 
# data:  gnomad_df
# p-value = 9.221e-06
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.108918 1.310407
# sample estimates:
# odds ratio 
#   1.205309 

df <- rbind(df, c("Constrained Genes (3570)", a1$estimate, a1$conf.int[1], a1$conf.int[2], a1$p.value,
                  gnomad_df[1, 1], gnomad_df[2, 1], gnomad_df[1, 2], gnomad_df[2, 2]))

a1_rr <- rateratio(sum(gnomad_genes_count_case$Total), sum(gnomad_genes_count_control$Total), 7749, 4025)
#           Cases Person-time
# Exposed    3272        7749
# Unexposed  1474        4025
# Total      4746       11774
a1_rr
# 	Incidence rate ratio estimate and its significance probability
# 
# data:  sum(gnomad_genes_count_case$Total) sum(gnomad_genes_count_control$Total) 7749 4025
# p-value = 5.553e-06
# 95 percent confidence interval:
#  1.084262 1.226134
# sample estimates:
# [1] 1.153018
df_rr <- rbind(df_rr, c("Constrained Genes (3570)", a1_rr$estimate, a1_rr$conf.int[1], a1_rr$conf.int[2], a1_rr$p.value,
                        sum(gnomad_genes_count_case$Total), sum(gnomad_genes_count_control$Total), 7749, 4025))






# SCHEMA
## Extract SCHEMA genes
schemagenes <- fread("SCHEMA_gene_results.tsv", data.table = F, h = T)
schemagenes <- schemagenes[complete.cases(schemagenes$`P meta`), ]
schema_genes <- schemagenes[schemagenes$`P meta`<1.30e-04,]$gene_id


## Subset gene_count down to gnomad_genes
schema_gene_names <- gnomad$gene[gnomad$gene_id %in% schema_genes]

schema_genes_count <- as.data.frame(t(gene_count[gene_count$"gene_symbol" %in% schema_gene_names, ]))
schema_genes_count <- schema_genes_count[-1, ]
schema_genes_count <- mutate_all(schema_genes_count, function(x) as.numeric(as.character(x)))

## Split into case/control
schema_genes_count_case <- schema_genes_count[samples_to_keep[samples_to_keep$AFFECTED_STATUS_FIXED == "CASE"]$s, ]
schema_genes_count_control <- schema_genes_count[samples_to_keep[samples_to_keep$AFFECTED_STATUS_FIXED == "CONTROL"]$s, ]

schema_genes_count_case <- schema_genes_count_case %>% mutate(Total = select(., V1:V27) %>% rowSums(na.rm = TRUE))
schema_genes_count_control <- schema_genes_count_control %>% mutate(Total = select(., V1:V27) %>% rowSums(na.rm = TRUE))

schema_df <- data.frame("rarevar" = c(sum(schema_genes_count_case$Total > 0),  sum(schema_genes_count_control$Total > 0)),
                        "novar" =   c(sum(schema_genes_count_case$Total == 0), sum(schema_genes_count_control$Total == 0)),
                        row.names = c("case", "control"))
schema_df
# rarevar novar
# case        115  7634
# control      39  3986

a2 <- fisher.test(schema_df)
a2
# 	Fisher's Exact Test for Count Data
# 
# data:  schema_df
# p-value = 0.02071
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.059841 2.279114
# sample estimates:
# odds ratio 
#   1.539585 

df <- rbind(df, c("Schizophrenia (32)", a2$estimate, a2$conf.int[1], a2$conf.int[2], a2$p.value,
                  schema_df[1, 1], schema_df[2, 1], schema_df[1, 2], schema_df[2, 2]))

a2_rr <- rateratio(sum(schema_genes_count_case$Total), sum(schema_genes_count_control$Total), 7749, 4025)
#           Cases Person-time
# Exposed     127        7749
# Unexposed    48        4025
# Total       175       11774
a2_rr
# 	Incidence rate ratio estimate and its significance probability
# 
# data:  sum(schema_genes_count_case$Total) sum(schema_genes_count_control$Total) 7749 4025
# p-value = 0.0595
# 95 percent confidence interval:
#  0.9859651 1.9155957
# sample estimates:
# [1] 1.374304
df_rr <- rbind(df_rr, c("Schizophrenia (32)", a2_rr$estimate, a2_rr$conf.int[1], a2_rr$conf.int[2], a2_rr$p.value,
                        sum(schema_genes_count_case$Total), sum(schema_genes_count_control$Total), 7749, 4025))










# ASD
## Extract ASD genes
asd <- read.table("asd.tada.exome.results.txt", header=T)
asd_genes <- asd[asd$FDR_TADA_ASD<0.001,]$gene

## Subset gene_count down to asd_genes
asd_genes_count <- as.data.frame(t(gene_count[gene_count$"gene_symbol" %in% asd_genes, ]))
asd_genes_count <- asd_genes_count[-1, ]
asd_genes_count <- mutate_all(asd_genes_count, function(x) as.numeric(as.character(x)))

## Split into case/control
asd_genes_count_case <- asd_genes_count[samples_to_keep[samples_to_keep$AFFECTED_STATUS_FIXED == "CASE"]$s, ]
asd_genes_count_control <- asd_genes_count[samples_to_keep[samples_to_keep$AFFECTED_STATUS_FIXED == "CONTROL"]$s, ]

asd_genes_count_case <- asd_genes_count_case %>% mutate(Total = select(., V1:V54) %>% rowSums(na.rm = TRUE))
asd_genes_count_control <- asd_genes_count_control %>% mutate(Total = select(., V1:V54) %>% rowSums(na.rm = TRUE))

asd_df <- data.frame("rarevar" = c(sum(asd_genes_count_case$Total > 0),  sum(asd_genes_count_control$Total > 0)),
                        "novar" =   c(sum(asd_genes_count_case$Total == 0), sum(asd_genes_count_control$Total == 0)),
                        row.names = c("case", "control"))
# rarevar novar
# case        115  7634
# control      45  3980

a3 <- fisher.test(asd_df)
a3
# 	Fisher's Exact Test for Count Data
# 
# data:  asd_df
# p-value = 0.111
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.9340239 1.9291022
# sample estimates:
# odds ratio 
#   1.332327 

df <- rbind(df, c("Autism Spectrum Disorder (72)", a3$estimate, a3$conf.int[1], a3$conf.int[2], a3$p.value,
                  asd_df[1, 1], asd_df[2, 1], asd_df[1, 2], asd_df[2, 2]))


a3_rr <- rateratio(sum(asd_genes_count_case$Total), sum(asd_genes_count_control$Total), 7749, 4025)
#           Cases Person-time
# Exposed     119        7749
# Unexposed    48        4025
# Total       167       11774
a3_rr
# 	Incidence rate ratio estimate and its significance probability
# 
# data:  sum(asd_genes_count_case$Total) sum(asd_genes_count_control$Total) 7749 4025
# p-value = 0.1381
# 95 percent confidence interval:
#  0.9210457 1.8004070
# sample estimates:
# [1] 1.287733
df_rr <- rbind(df_rr, c("Autism Spectrum Disorder (72)", a3_rr$estimate, a3_rr$conf.int[1], a3_rr$conf.int[2], a3_rr$p.value,
                        sum(asd_genes_count_case$Total), sum(asd_genes_count_control$Total), 7749, 4025))









# NDD
## Extract NDD genes
ndd_genes <- asd[asd$FDR_TADA_NDD<0.001,]$gene

## Subset gene_count down to asd_genes
ndd_genes_count <- as.data.frame(t(gene_count[gene_count$"gene_symbol" %in% ndd_genes, ]))
ndd_genes_count <- ndd_genes_count[-1, ]
ndd_genes_count <- mutate_all(ndd_genes_count, function(x) as.numeric(as.character(x)))

## Split into case/control
ndd_genes_count_case <- ndd_genes_count[samples_to_keep[samples_to_keep$AFFECTED_STATUS_FIXED == "CASE"]$s, ]
ndd_genes_count_control <- ndd_genes_count[samples_to_keep[samples_to_keep$AFFECTED_STATUS_FIXED == "CONTROL"]$s, ]

ndd_genes_count_case <- ndd_genes_count_case %>% mutate(Total = select(., V1:V274) %>% rowSums(na.rm = TRUE))
ndd_genes_count_control <- ndd_genes_count_control %>% mutate(Total = select(., V1:V274) %>% rowSums(na.rm = TRUE))

ndd_df <- data.frame("rarevar" = c(sum(ndd_genes_count_case$Total > 0),  sum(ndd_genes_count_control$Total > 0)),
                     "novar" =   c(sum(ndd_genes_count_case$Total == 0), sum(ndd_genes_count_control$Total == 0)),
                     row.names = c("case", "control"))
# rarevar novar
# case        516  7233
# control     201  3824

a4 <- fisher.test(ndd_df)
a4
# 	Fisher's Exact Test for Count Data
# 
# data:  ndd_df
# p-value = 0.0002964
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.145267 1.613087
# sample estimates:
# odds ratio 
#   1.357195 

df <- rbind(df, c("Neurodevelopmental Disorder (373)", a4$estimate, a4$conf.int[1], a4$conf.int[2], a4$p.value,
                  ndd_df[1, 1], ndd_df[2, 1], ndd_df[1, 2], ndd_df[2, 2]))


a4_rr <- rateratio(sum(ndd_genes_count_case$Total), sum(ndd_genes_count_control$Total), 7749, 4025)
#           Cases Person-time
# Exposed     560        7749
# Unexposed   213        4025
# Total       773       11774
a4_rr
# 	Incidence rate ratio estimate and its significance probability
# 
# data:  sum(ndd_genes_count_case$Total) sum(ndd_genes_count_control$Total) 7749 4025
# p-value = 0.0001017
# 95 percent confidence interval:
#  1.166287 1.599013
# sample estimates:
# [1] 1.365616
df_rr <- rbind(df_rr, c("Neurodevelopmental Disorder (373)", a4_rr$estimate, a4_rr$conf.int[1], a4_rr$conf.int[2], a4_rr$p.value,
                        sum(ndd_genes_count_case$Total), sum(ndd_genes_count_control$Total), 7749, 4025))









# All genes
all_count <- as.data.frame(t(gene_count))
all_count <- all_count[-1,]
all_count <- mutate_all(all_count, function(x) as.numeric(as.character(x)))
all_count_case <- all_count[samples_to_keep[samples_to_keep$AFFECTED_STATUS_FIXED == "CASE"]$s, ]
all_count_control <- all_count[samples_to_keep[samples_to_keep$AFFECTED_STATUS_FIXED == "CONTROL"]$s, ]

all_count_case <- all_count_case %>% mutate(Total = select(., V1:V18356) %>% rowSums(na.rm = TRUE))
all_count_control <- all_count_control %>% mutate(Total = select(., V1:V18356) %>% rowSums(na.rm = TRUE))

all_df <- data.frame("rarevar" = c(sum(all_count_case$Total > 0),  sum(all_count_control$Total > 0)),
                     "novar" =   c(sum(all_count_case$Total == 0), sum(all_count_control$Total == 0)),
                     row.names = c("case", "control"))
# rarevar novar
# case        7729  20
# control     4005  20

a5 <- fisher.test(all_df)
a5
# 	Fisher's Exact Test for Count Data
# 
# data:  all_df
# p-value = 0.04405
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.9840734 3.7841875
# sample estimates:
# odds ratio 
#   1.929788 

a5_rr <- rateratio(sum(all_count_case$Total), sum(all_count_control$Total), 7749, 4025)
#            Cases Person-time
# Exposed    71851        7749
# Unexposed  36544        4025
# Total     108395       11774
a5_rr
# 	Incidence rate ratio estimate and its significance probability
# 
# data:  sum(all_count_case$Total) sum(all_count_control$Total) 7749 4025
# p-value = 0.001058
# 95 percent confidence interval:
#  1.008481 1.034204
# sample estimates:
# [1] 1.021261
df_rr <- rbind(df_rr, c("All genes (18356)", a5_rr$estimate, a5_rr$conf.int[1], a5_rr$conf.int[2], a5_rr$p.value,
                        sum(all_count_case$Total), sum(all_count_control$Total), 7749, 4025))

df <- rbind(df, c("All genes (18356)", a5_rr$estimate, a5_rr$conf.int[1], a5_rr$conf.int[2], a5_rr$p.value,
                     sum(all_count_case$Total), sum(all_count_control$Total), 7749, 4025))
colnames(df) <- c("Category", 
                  "OR", "CI_low", "CI_high", "pval", 
                  "rarevar_case", "rarevar_control", "novar_case", "novar_control")
colnames(df_rr) <- c("Category", 
                  "OR", "CI_low", "CI_high", "pval", 
                  "rarevar_case", "rarevar_control", "novar_case", "novar_control")

write.table(df, "/Users/rye/Projects/2024_WGSPD/2024_WCPG/plotting/inputs/20240929_ptv-ac5_old-analysis_Fisher.tsv",
            sep = "\t", col.names = T, row.names = F, quote = F)
write.table(df_rr, "/Users/rye/Projects/2024_WGSPD/2024_WCPG/plotting/inputs/20240929_ptv-ac5_old-analysis_rate-ratio.tsv",
            sep = "\t", col.names = T, row.names = F, quote = F)
