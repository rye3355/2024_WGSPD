# In neale-bipex-wes terra workspace
# gs://fc-54cd2a03-28fe-43ab-9142-1d265515b386/WGSPD/202411_SV/joined.tsv.gz
library(data.table)

system2("mkdir", "-p /home/rstudio/tmp/20241113_WGSPD_SV/")
system2("gsutil", "cp gs://fc-54cd2a03-28fe-43ab-9142-1d265515b386/WGSPD/202411_SV/joined.tsv.gz /home/rstudio/tmp/20241113_WGSPD_SV/")

dir <- "/home/rstudio/tmp/20241113_WGSPD_SV/"
d <- fread(paste0(dir, "joined.tsv.gz")) # 2063885
colnames(d) <- c("CHROM",	"POS",	"ID",	"REF",	"ALT",	"QUAL",	"FILTER",	"INFO")




table(d$ALT, useNA = "always")
#  <BND>          <CNV>          <CPX>          <CTX> <DEL:ME:HERVK> <DEL:ME:LINE1>          <DEL> 
# 343492            657          14464             98            647           7724        1146928 
#  <DUP>   <INS:ME:ALU> <INS:ME:LINE1>   <INS:ME:SVA>          <INS>          <INV>           <NA> 
# 253146         168462          28891          17170          80080           2126              0



# Filter to pass
table(d$FILTER, useNA = "always")
#                           FAIL_MANUAL_REVIEW                  FAIL_MANUAL_REVIEW;HIGH_NCR 
#                                           67                                            3 
#                                     HIGH_NCR                     HIGH_NCR;IGH_MHC_OVERLAP 
#                                        79097                                          514 
# HIGH_NCR;IGH_MHC_OVERLAP;LOWQUAL_WHAM_SR_DEL          HIGH_NCR;IGH_MHC_OVERLAP;UNRESOLVED 
#                                          493                                         1624 
#                 HIGH_NCR;LOWQUAL_WHAM_SR_DEL                          HIGH_NCR;UNRESOLVED 
#                                        67124                                        73518 
#                              IGH_MHC_OVERLAP          IGH_MHC_OVERLAP;LOWQUAL_WHAM_SR_DEL 
#                                         5424                                          882 
#                   IGH_MHC_OVERLAP;UNRESOLVED                          LOWQUAL_WHAM_SR_DEL 
#                                         7280                                       126104 
#  LOWQUAL_WHAM_SR_DEL;OUTLIER_SAMPLE_ENRICHED                      OUTLIER_SAMPLE_ENRICHED 
#                                       181528                                       105606 
#                                         PASS                           REFERENCE_ARTIFACT 
#                                      1144040                                           55 
#                                   UNRESOLVED                                         <NA> 
#                                       270526                                            0 
d <- d[d$FILTER == "PASS",] # 1144040

d$split <- strsplit(d$INFO, ";")
pattern <- "PREDICTED_.*=.*"
extract_matches <- function(row) {
  grep(pattern, row, value = TRUE)  # Searches and returns the values that match the pattern
}
cpx$matches <- apply(cpx, 1, function(x) extract_matches(x$split))







# Things we want from info:
# PREDICTED_* things
# SVLEN
# SVTYPE (check against ALT)
# AC




