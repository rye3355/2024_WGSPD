# Look for enrichment of ASD genes associated SVs now with Fin Controls
library(data.table)
library(fmsb)

out_root = "~/tmp/20240930_WGSPD/"
system2("mkdir", paste("-p", out_root))
#system2("gsutil", paste("cp gs://fc-54cd2a03-28fe-43ab-9142-1d265515b386/WGSPD/20231004_SV/gnomAD_SV_v3.releasable.WGSDP-Fin-CTRL.with_annotations_PASSING.bed.gz", out_root))
#system2("gsutil", paste("cp gs://fc-54cd2a03-28fe-43ab-9142-1d265515b386/WGSPD/20230822_SNV-non-coding/gene-sets/asd.tada.exome.results.txt", out_root))
#system2("gsutil", paste("cp gs://fc-54cd2a03-28fe-43ab-9142-1d265515b386/WGSPD/20230822_SNV-non-coding/gene-sets/gnomad.v2.1.1.lof_metrics.by_gene.txt.gz", out_root))
#system2("gsutil", paste("cp gs://fc-54cd2a03-28fe-43ab-9142-1d265515b386/WGSPD/20230822_SNV-non-coding/gene-sets/SCHEMA_gene_results.tsv", out_root))
#system2("gsutil", paste("cp gs://fc-54cd2a03-28fe-43ab-9142-1d265515b386/WGSPD/20231004_SV/gencode.v44.annotation.gtf.gz", out_root))



# Process input files
## Bed
merged <- fread(paste0(out_root, "gnomAD_SV_v3.releasable.WGSDP-Fin-CTRL.with_annotations_PASSING.bed.gz")) # 555528

### Filter to DELETIONS
p <- merged[merged$FILTER == "PASS" & 
              complete.cases(merged$AC_withFin) & 
              merged$AC_withFin <= 5 & 
              merged$svtype %in% c("INS:ME:ALU", "INS:ME:LINE1", "INS:ME:SVA"), ] 
# AC <= 5: 555528

### Find all PREDICTED fields
predicted_cols <- colnames(p)[grepl("PREDICTED", colnames(p), fixed = TRUE)]
### But we only want some
#predicted_cols <- c("PREDICTED_INTERGENIC", "PREDICTED_INTRONIC", "PREDICTED_LOF", "PREDICTED_NEAREST_TSS", "PREDICTED_PROMOTER", "PREDICTED_UTR")

predicted_cols <- c("PREDICTED_INTERGENIC", "PREDICTED_INTRONIC", "PREDICTED_LOF", "PREDICTED_PROMOTER", "PREDICTED_UTR")

# Helper function to see if at least one gene of interest shows up in entry
contains_genes <- function(entry) {
  a <- strsplit(entry, ",")[[1]]
  return(any(genes %in% a))
}

case_total = 6336
control_total = 2793 + 2484
system2("mkdir", paste0("-p ", out_root, "outputs/"))


test <- function(p, predicted_cols, genes, geneset) {
  df = data.frame(matrix(vector(), 0, 13),
                  stringsAsFactors=F)
  df_rr = data.frame(matrix(vector(), 0, 13),
                     stringsAsFactors=F)
  for (c in predicted_cols) {
    if (c == "PREDICTED_NEAREST_TSS") {
      next
    }
    cat(c, "\n")
    
    ## Length threshold
    l = 100
    if(c == "PREDICTED_INTRONIC") {
      l = 5000
    }
    
    ## Subset to SVs above length threshold
    s <- p[p$SVLEN > l, ]
    cat(dim(s), "\n")
    ## Check for NAs
    if(all(is.na(s[[c]]))) {
      cat("PREDICTED all NA entries", "\n", "\n")
      next
    }
    
    s <- s[!is.na((s[[c]])),]
    o <- dim(s)[1]
    if(!any(sapply(s[[c]], is.character))) {
      cat("PREDICTED no gene symbol information", "\n")
      next
    }
    s <- s[sapply(s[[c]][sapply(s[[c]], is.character)], contains_genes),]
    cat(paste0("Number of >", l, "b DELs within ", length(genes), geneset, " genes: ", dim(s)[1], " (out of ", o, ")"), "\n")
    
    a <- sum(s$CASE_N_HET[!is.na(s$CASE_N_HET)], s$CASE_N_HOMALT[!is.na(s$CASE_N_HOMALT)])
    b <- sum(s$CONTROL_N_HET_withFin[!is.na(s$CONTROL_N_HET_withFin)], s$CONTROL_N_HOMALT_withFin[!is.na(s$CONTROL_N_HOMALT_withFin)])
    corrected <- F
    both <- F
    if(a == 0 | b == 0) {
      corrected <- T
      if(a == 0 & b == 0) {
        both <- T
      }
      a <- a + 1
      b <- b + 1
    }
    dat <- data.frame(
      "rarevar" = c(a, b),
      "novar" = c(case_total - a, control_total - b),
      row.names = c("case", "control"),
      stringsAsFactors = FALSE
    )
    print(dat)
    f <- fisher.test(dat)
    print(f)
    df <- rbind(df, c(geneset, c, f$estimate, f$conf.int[1], f$conf.int[2], f$p.value,
                      a, b, case_total - a, control_total - b, l, corrected, both))
    
    rr <- rateratio(a, b, case_total, control_total)
    print(rr)
    df_rr <- rbind(df_rr, c(geneset, c, rr$estimate, rr$conf.int[1], rr$conf.int[2], rr$p.value,
                            a, b, case_total, control_total, l, corrected, both))
    cat("\n", "\n", "\n")
  }
  colnames(df) <- c("Gene_set", "Category", "OR", "CI_low", "CI_high", "pval", 
                    "rarevar_case", "rarevar_control", "novar_case", "novar_control",
                    "length_thresh", "corrected", "both_corrected")
  colnames(df_rr) <- c("Gene_set", "Category", "RR", "CI_low", "CI_high", "pval", 
                       "rarevar_case", "rarevar_control", "novar_case", "novar_control",
                       "length_thresh", "corrected", "both_corrected")
  return(list(fisher = df, rr = df_rr))
}



df_fisher = data.frame(matrix(vector(), 0, 13),
                       stringsAsFactors=F)
df_rateratio = data.frame(matrix(vector(), 0, 13),
                          stringsAsFactors=F)

############################################################
## ASD
############################################################
ASD <- fread(paste0(out_root, "asd.tada.exome.results.txt"))
genes <- ASD[ASD$FDR_TADA_ASD<0.001,]$gene

sink(paste0(out_root, "outputs/ASD_PREDICTED_del-ac5.log"))
dfs <- test(p, predicted_cols, genes, "ASD")
sink()
df_fisher <- rbind(df_fisher, dfs$fisher)
df_rateratio <- rbind(df_rateratio, dfs$rr)
############################################################
## NDD
############################################################
NDD <- fread(paste0(out_root, "asd.tada.exome.results.txt"))
genes <- NDD[NDD$FDR_TADA_NDD<0.001,]$gene

sink(paste0(out_root, "outputs/NDD_PREDICTED_del-ac5.log"))
dfs <- test(p, predicted_cols, genes, "NDD")
sink()
df_fisher <- rbind(df_fisher, dfs$fisher)
df_rateratio <- rbind(df_rateratio, dfs$rr)
############################################################
## pLI-constrained
############################################################
gnomad <- fread(paste0(out_root, "gnomad.v2.1.1.lof_metrics.by_gene.txt.gz"))
gnomad[gnomad$gene == "DOPEY1",]$gene <- "DOP1A"
genes <- gnomad[gnomad$pLI>.9,]$gene

sink(paste0(out_root, "outputs/pLI-constrained_PREDICTED_del-ac5.log"))
dfs <- test(p, predicted_cols, genes, "pLI-constrained")
sink()
df_fisher <- rbind(df_fisher, dfs$fisher)
df_rateratio <- rbind(df_rateratio, dfs$rr)
############################################################
## SCHEMA
############################################################
SCHEMA <- fread(paste0(out_root, "SCHEMA_gene_results.tsv"))
SCHEMA <- SCHEMA[order(SCHEMA$`P meta`),]
SCHEMA <- SCHEMA[c(1:32), ]
SCHEMA <- merge(SCHEMA, gnomad[, c("gene_id", "gene")],
                by = "gene_id")
genes <- SCHEMA$gene

sink(paste0(out_root, "outputs/SCHEMA_PREDICTED_del-ac5.log"))
dfs <- test(p, predicted_cols, genes, "SCHEMA")
sink()
df_fisher <- rbind(df_fisher, dfs$fisher)
df_rateratio <- rbind(df_rateratio, dfs$rr)








write.table(df_fisher, file = paste0(out_root, "outputs/fisher_del-ac5_table.tsv"), append = T, quote = F, 
            sep = "\t", col.names = F, row.names = F)
write.table(df_rateratio, file = paste0(out_root, "outputs/rate-ratio_del-ac5_table.tsv"), append = T, quote = F, 
            sep = "\t", col.names = F, row.names = F)


system2("gsutil", paste("cp -r",
                        paste0(out_root, "outputs/*"),
                        "gs://fc-5571d2c9-c15e-466a-9ee3-d0c0fee335aa/2024_WCPG/del/"))
