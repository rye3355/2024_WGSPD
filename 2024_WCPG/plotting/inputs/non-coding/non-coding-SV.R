#/Users/rye/Projects/WGSPD/gnomad-v3-subset/20231004_SV/02_fin-ctrl_noncoding
library(data.table)
library(fmsb)

out_root = "~/tmp/20240930_WGSPD/"
system2("mkdir", paste("-p", out_root))
#system2("gsutil", paste("cp gs://fc-54cd2a03-28fe-43ab-9142-1d265515b386/WGSPD/20231004_SV/gnomAD_SV_v3.releasable.WGSDP-Fin-CTRL.with_annotations_PASSING.bed.gz", out_root))
#system2("gsutil", paste("cp gs://fc-54cd2a03-28fe-43ab-9142-1d265515b386/WGSPD/20230822_SNV-non-coding/gene-sets/asd.tada.exome.results.txt", out_root))
#system2("gsutil", paste("cp gs://fc-54cd2a03-28fe-43ab-9142-1d265515b386/WGSPD/20230822_SNV-non-coding/gene-sets/gnomad.v2.1.1.lof_metrics.by_gene.txt.gz", out_root))
#system2("gsutil", paste("cp gs://fc-54cd2a03-28fe-43ab-9142-1d265515b386/WGSPD/20230822_SNV-non-coding/gene-sets/SCHEMA_gene_results.tsv", out_root))
#system2("gsutil", paste("cp gs://fc-54cd2a03-28fe-43ab-9142-1d265515b386/WGSPD/20231004_SV/gencode.v44.annotation.gtf.gz", out_root))
#system2("gsutil", paste("cp gs://fc-54cd2a03-28fe-43ab-9142-1d265515b386/WGSPD/20231005_noncoding-constraint/constraint_z_genome_1kb_filtered.ft17-4-9.nc.sorted.copy.bed.gz", out_root))



# Process input files
## Bed
merged <- fread(paste0(out_root, "gnomAD_SV_v3.releasable.WGSDP-Fin-CTRL.with_annotations_PASSING.bed.gz")) # 555528



# Process noncoding file
constrained <- fread(paste0(out_root, "constraint_z_genome_1kb_filtered.ft17-4-9.nc.sorted.copy.bed.gz"))
colnames(constrained) <- c("Chr", "Start", "End", "z")
constrained$ID <- seq(c(1:dim(constrained)[1]))
constrained <- constrained[constrained$z >= 4,] # 19471
write.table(constrained, paste0(out_root, "noncoding-constraint_z4.input"),
            quote = F, sep = "\t", row.names = F, col.names = F)



system2("mkdir", paste0("-p ", out_root, "outputs/"))

df_fisher = data.frame(matrix(vector(), 0, 9),
                       stringsAsFactors=F)
df_rateratio = data.frame(matrix(vector(), 0, 9),
                          stringsAsFactors=F)
case_total = 6336
control_total = 2793 + 2484
############################################################
## DEL
############################################################
p <- merged[merged$FILTER == "PASS" & 
              complete.cases(merged$AC_withFin) & complete.cases(merged$CASE_AC) &
              merged$AC_withFin <= 2 & 
              merged$svtype == "DEL" &
              merged$SVLEN > 5000, ] 
# AC <= 10: 206650
# AC <= 5: 174927
# AC <= 10 and > 100bp: 175655
# AC <= 10 and > 500bp: 106832
# AC <= 5 and >500bp: 93132
# AC <= 5 and >1kb: 66336
# AC <= 5 and >5kb: 36015
# AC <= 2 and >5kb: 28428
write.table(p[, c("#chrom", "start", "end", "name", 
                  "CASE_N_HET", "CASE_N_HOMALT", "CONTROL_N_HET_withFin", "CONTROL_N_HOMALT_withFin")],  
            paste0(out_root, "gnomAD_SV_del-ac2-5kb.input"),
            quote = F, sep = "\t", row.names = F, col.names = F)
system2("./bedtools", paste("intersect -F 1 -wa -wb -a",
                            paste0(out_root, "gnomAD_SV_del-ac2-5kb.input"),
                            "-b", 
                            paste0(out_root, "noncoding-constraint_z4.input"),
                            ">",
                            paste0(out_root, "outputs/gnomAD_SV_del-ac2-5kb_intersect-noncoding-constraint-z4.output")))


constrained_res <- fread(paste0(out_root, "outputs/gnomAD_SV_del-ac2-5kb_intersect-noncoding-constraint-z4.output"))
colnames(constrained_res) <- c(c("#chrom", "start", "end", "name", 
                                 "CASE_N_HET", "CASE_N_HOMALT", "CONTROL_N_HET_withFin", "CONTROL_N_HOMALT_withFin"), colnames(constrained))
length(unique(constrained_res$ID)) # 2974 / 19471 constrained regions represented
length(unique(constrained_res$name)) # 1375 / 28428 unique SVs represented

counts <- p[p$name %in% constrained_res$name, c("CASE_N_HET", "CASE_N_HOMALT", "CONTROL_N_HET_withFin", "CONTROL_N_HOMALT_withFin")]

a <- sum(counts$CASE_N_HET, counts$CASE_N_HOMALT)
b <- sum(counts$CONTROL_N_HET_withFin, counts$CONTROL_N_HOMALT_withFin)
dat <- data.frame(
  "rarevar" = c(a, b),
  "novar" = c(case_total - a, control_total - b),
  row.names = c("case", "control"),
  stringsAsFactors = FALSE
)
#         rarevar novar
# case       1144  5192
# control     533  4744
f <- fisher.test(dat)
f
# 	Fisher's Exact Test for Count Data
# 
# data:  dat
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.754644 2.193576
# sample estimates:
# odds ratio 
#   1.961015 

df_fisher <- rbind(df_fisher, c("DEL", f$estimate, f$conf.int[1], f$conf.int[2], f$p.value,
                                a, b, case_total - a, control_total - b))

a2_rr <- rateratio(a, b, case_total, control_total)
#           Cases Person-time
# Exposed    1144        6336
# Unexposed   533        5277
# Total      1677       11613
a2_rr
# 	Incidence rate ratio estimate and its significance probability
# 
# data:  a b case_total control_total
# p-value < 2.2e-16
# 95 percent confidence interval:
#  1.612987 1.981119
# sample estimates:
# [1] 1.787602
df_rateratio <- rbind(df_rateratio, c("DEL", a2_rr$estimate, a2_rr$conf.int[1], a2_rr$conf.int[2], a2_rr$p.value,
                                      a, b, case_total, control_total))


############################################################
## DUP
############################################################
p <- merged[merged$FILTER == "PASS" & 
              complete.cases(merged$AC_withFin) & complete.cases(merged$CASE_AC) &
              merged$AC_withFin <= 2 & 
              merged$svtype == "DUP" &
              merged$SVLEN > 5000, ] 
# AC <= 5 and >5kb: 16038
# AC <= 2 and >5kb: 13536
write.table(p[, c("#chrom", "start", "end", "name", 
                  "CASE_N_HET", "CASE_N_HOMALT", "CONTROL_N_HET_withFin", "CONTROL_N_HOMALT_withFin")],  
            paste0(out_root, "gnomAD_SV_dup-ac2-5kb.input"),
            quote = F, sep = "\t", row.names = F, col.names = F)
system2("./bedtools", paste("intersect -F 1 -wa -wb -a",
                            paste0(out_root, "gnomAD_SV_dup-ac2-5kb.input"),
                            "-b", 
                            paste0(out_root, "noncoding-constraint_z4.input"),
                            ">",
                            paste0(out_root, "outputs/gnomAD_SV_dup-ac2-5kb_intersect-noncoding-constraint-z4.output")))


constrained_res <- fread(paste0(out_root, "outputs/gnomAD_SV_dup-ac2-5kb_intersect-noncoding-constraint-z4.output"))
colnames(constrained_res) <- c(c("#chrom", "start", "end", "name", 
                                 "CASE_N_HET", "CASE_N_HOMALT", "CONTROL_N_HET_withFin", "CONTROL_N_HOMALT_withFin"), colnames(constrained))
length(unique(constrained_res$ID)) # 6005 / 19471 constrained regions represented
length(unique(constrained_res$name)) # 1564 / 13536 unique SVs represented

counts <- p[p$name %in% constrained_res$name, c("CASE_N_HET", "CASE_N_HOMALT", "CONTROL_N_HET_withFin", "CONTROL_N_HOMALT_withFin")]

a <- sum(counts$CASE_N_HET, counts$CASE_N_HOMALT)
b <- sum(counts$CONTROL_N_HET_withFin, counts$CONTROL_N_HOMALT_withFin)
dat <- data.frame(
  "rarevar" = c(a, b),
  "novar" = c(case_total - a, control_total - b),
  row.names = c("case", "control"),
  stringsAsFactors = FALSE
)
#        rarevar novar
# case       1246  5090
# control     570  4707
f <- fisher.test(dat)
f
# 	Fisher's Exact Test for Count Data
# 
# data:  dat
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.814808 2.253411
# sample estimates:
# odds ratio 
#   2.021376 

df_fisher <- rbind(df_fisher, c("DUP", f$estimate, f$conf.int[1], f$conf.int[2], f$p.value,
                                a, b, case_total, control_total))

a2_rr <- rateratio(a, b, case_total, control_total)
#           Cases Person-time
# Exposed    1246        6336
# Unexposed   570        5277
# Total      1816       11613
a2_rr
# 	Incidence rate ratio estimate and its significance probability
# 
# data:  a b case_total control_total
# p-value < 2.2e-16
# 95 percent confidence interval:
#  1.648819 2.010283
# sample estimates:
# [1] 1.820602
df_rateratio <- rbind(df_rateratio, c("DUP", a2_rr$estimate, a2_rr$conf.int[1], a2_rr$conf.int[2], a2_rr$p.value,
                                      a, b, case_total - a, control_total - b))




colnames(df_fisher) <- c("Category", 
                         "OR", "CI_low", "CI_high", "pval", 
                         "rarevar_case", "rarevar_control", "novar_case", "novar_control")
colnames(df_rateratio) <- c("Category", 
                            "RR", "CI_low", "CI_high", "pval", 
                            "rarevar_case", "rarevar_control", "case_total", "control_total")

write.table(df_fisher, file = paste0(out_root, "outputs/fisher_non-coding_table.tsv"), quote = F, 
            sep = "\t", col.names = T, row.names = F)
write.table(df_rateratio, file = paste0(out_root, "outputs/rate-ratio_non-coding_table.tsv"),quote = F, 
            sep = "\t", col.names = T, row.names = F)


system2("gsutil", paste("cp",
                        paste0(out_root, "outputs/*.tsv"),
                        "gs://fc-5571d2c9-c15e-466a-9ee3-d0c0fee335aa/2024_WCPG/non-coding/"))
