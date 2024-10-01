# Look for enrichment of ASD genes associated SVs now with Fin Controls
library(data.table)

# Process input files
## Bed
merged <- fread("/Users/rye/Projects/WGSPD/gnomad-v3-subset/20231004_SV/gnomAD_SV_v3.releasable.WGSDP.with_annotations.bed.gz") # 555528


table(merged$svtype)
#        BND          CNV          CPX          CTX          DEL DEL:ME:HERVK DEL:ME:LINE1          DUP          INS 
#     251914          719         4987           14       795028          336         3910       113000        52442 
# INS:ME:ALU INS:ME:LINE1   INS:ME:SVA          INV 
#      66426        18050         6911          576 



## Add in length annotation
merged$length <- as.integer(merged$end) - as.integer(merged$start)

### Filter to DELETIONS
p <- merged[merged$FILTER == "PASS" & 
              complete.cases(merged$AC_withFin) & 
              merged$AC_withFin <= 5 & 
              merged$svtype == "DEL", ] 
# AC <= 10: 206650
# AC <= 5: 174927
# AC <= 10 and > 100bp: 175655
# AC <= 10 and > 500bp: 106832
# AC <= 5 and >100bp: 147453
# AC <= 5 and >500bp: 93132
# AC <= 5 and >5kb: 36015

### Find all PREDICTED fields
predicted_cols <- colnames(p)[grepl("PREDICTED", colnames(p), fixed = TRUE)]
### But we only want some
#predicted_cols <- c("PREDICTED_INTRONIC", "PREDICTED_LOF", "PREDICTED_PROMOTER", "PREDICTED_UTR")

## ASD
ASD <- fread("~/tmp/asd.tada.exome.results.txt")
genes <- ASD[ASD$FDR_TADA_ASD<0.001,]$gene


# Helper function to see if at least one gene of interest shows up in entry
contains_genes <- function(entry) {
  a <- strsplit(entry, ",")[[1]]
  return(any(genes %in% a))
}

df = data.frame(matrix(vector(), 0, 13),
                stringsAsFactors=F)
case_total = 6336
control_total = 2793 + 2484
sink(paste0(out_root, "ASD_PREDICTED_del-ac5.log"))
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
  s <- p[p$length > l, ]
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
  cat(paste0("Number of >", l, "b DELs within ", length(genes), " ASD genes: ", dim(s)[1], " (out of ", o, ")"), "\n")
  
  
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
  if(any(dat < 0)) {
    cat("Carrier heuristic failed", "\n")
    next
  }
  f <- fisher.test(dat)
  print(f)
  
  df <- rbind(df, c("ASD", c, f$estimate, f$conf.int[1], f$conf.int[2], f$p.value,
                    a, b, case_total - a, control_total - b, l, corrected, both))
  cat("\n", "\n", "\n")
  
}
sink()

colnames(df) <- c("Gene_set", "Category", "OR", "CI_low", "CI_high", "pval", 
                  "rarevar_case", "rarevar_control", "novar_case", "novar_control",
                  "length_thresh", "corrected", "both_corrected")
write.table(df, file = paste0(out_root, "gene-set_del-ac5_table.tsv"), append = T, quote = F, 
            sep = "\t", col.names = F, row.names = F)
