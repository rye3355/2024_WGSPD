
library(data.table)

d <- fread("/Users/rye/Projects/2024_WGSPD/Analysis/SV/WGSDP_Releasable_4.1_202405/joined.tsv.gz") # 2063922
colnames(d) <- c("CHROM",	"POS",	"ID",	"REF",	"ALT",	"QUAL",	"FILTER",	"INFO")




table(d$ALT, useNA = "always")
#        <BND>          <CNV>          <CPX>          <CTX> <DEL:ME:HERVK> <DEL:ME:LINE1>          <DEL>          <DUP> 
#       343502            659          14464             98            647           7724        1146943         253156 
# <INS:ME:ALU> <INS:ME:LINE1>   <INS:ME:SVA>          <INS>          <INV>           <NA> 
#       168462          28891          17170          80080           2126              0 

cpx <- d[d$ALT == "<CPX>" & d$FILTER == "PASS",]

# Filter to pass
table(d$FILTER, useNA = "always")
#                           FAIL_MANUAL_REVIEW                  FAIL_MANUAL_REVIEW;HIGH_NCR 
#                                           67                                            3 
#                                     HIGH_NCR                     HIGH_NCR;IGH_MHC_OVERLAP 
#                                        79103                                          514 
# HIGH_NCR;IGH_MHC_OVERLAP;LOWQUAL_WHAM_SR_DEL          HIGH_NCR;IGH_MHC_OVERLAP;UNRESOLVED 
#                                          493                                         1624 
#                 HIGH_NCR;LOWQUAL_WHAM_SR_DEL                          HIGH_NCR;UNRESOLVED 
#                                        67125                                        73523 
#                              IGH_MHC_OVERLAP          IGH_MHC_OVERLAP;LOWQUAL_WHAM_SR_DEL 
#                                         5424                                          882 
#                   IGH_MHC_OVERLAP;UNRESOLVED                          LOWQUAL_WHAM_SR_DEL 
#                                         7280                                       126104 
#  LOWQUAL_WHAM_SR_DEL;OUTLIER_SAMPLE_ENRICHED                      OUTLIER_SAMPLE_ENRICHED 
#                                       181528                                       105606 
#                                         PASS                           REFERENCE_ARTIFACT 
#                                      1144060                                           55 
#                                   UNRESOLVED                                         <NA> 
#                                       270531                                            0 
d <- d[d$FILTER == "PASS",]

cpx$split <- strsplit(cpx$INFO, ";")
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




