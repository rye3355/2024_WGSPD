#/Users/rye/Projects/WGSPD/gnomad-v3-subset/20231005_noncoding-constraint/04_noncoding-constraint_counts-export.py
library(data.table)
library(fmsb)


df_fisher = data.frame(matrix(vector(), 0, 9),
                       stringsAsFactors=F)
df_rateratio = data.frame(matrix(vector(), 0, 9),
                          stringsAsFactors=F)
case_total = 7749
control_total = 4025


a <- 7749
b <- 4063
dat <- data.frame(
  "rarevar" = c(a, b),
  "novar" = c(5, 2),
  row.names = c("case", "control"),
  stringsAsFactors = FALSE
)
#         rarevar novar
# case       7844     5
# control    4063     2
f <- fisher.test(dat)
f
# 	Fisher's Exact Test for Count Data
# 
# data:  dat
# p-value = 1
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.07260709 4.66238938
# sample estimates:
# odds ratio 
#  0.7628839 

df_fisher <- rbind(df_fisher, c("SNV", f$estimate, f$conf.int[1], f$conf.int[2], f$p.value,
                                a, b, 5, 2))

a2_rr <- rateratio(a, b, case_total, control_total)
#           Cases Person-time
# Exposed    7749        7749
# Unexposed  4063        4025
# Total     11812       11774
a2_rr
# 	Incidence rate ratio estimate and its significance probability
# 
# data:  a b case_total control_total
# p-value = 0.6276
# 95 percent confidence interval:
#  0.953744 1.028978
# sample estimates:
# [1] 0.9906473
df_rateratio <- rbind(df_rateratio, c("SNV", a2_rr$estimate, a2_rr$conf.int[1], a2_rr$conf.int[2], a2_rr$p.value,
                                      a, b, case_total, control_total))




colnames(df_fisher) <- c("Category", 
                         "OR", "CI_low", "CI_high", "pval", 
                         "rarevar_case", "rarevar_control", "novar_case", "novar_control")
colnames(df_rateratio) <- c("Category", 
                            "RR", "CI_low", "CI_high", "pval", 
                            "rarevar_case", "rarevar_control", "case_total", "control_total")

write.table(df_fisher, file = "fisher_SNV_non-coding_table.tsv", quote = F, 
            sep = "\t", col.names = T, row.names = F)
write.table(df_rateratio, file = "rate-ratio_SNV_non-coding_table.tsv",quote = F, 
            sep = "\t", col.names = T, row.names = F)

