library(data.table)
library(plyr)
d <- fread("/Users/rye/Projects/2024_WGSPD/Analysis/SV/WGSDP_Releasable_4.1_202405/joined.tsv.gz") # (2063922)
colnames(d) <- c("CHROM",	"POS",	"ID",	"REF",	"ALT",	"QUAL",	"FILTER",	"INFO")


ins <- d[d$ALT == "<INS>" & d$FILTER == "PASS",] # Filter to high quality INS (72821)




#######################################################################################################
# Processing INFO field
#######################################################################################################
extract_matches <- function(row, pattern) {
  grep(pattern, row, value = TRUE)  # Searches and returns the values that match the pattern
}
ins$split <- strsplit(ins$INFO, ";")

# Double check SVTYPE is as expected
ins$matches <- gsub("SVTYPE=", "", apply(ins, 1, function(x) extract_matches(x$split, "SVTYPE=.*")))
stopifnot(all(ins$matches == "INS"))


# Get SVLEN 
ins$matches <- gsub("SVLEN=", "", apply(ins, 1, function(x) extract_matches(x$split, "SVLEN=.*")))
ins$SVLEN <- as.numeric(ins$matches)
stopifnot(!any(is.na(ins$SVLEN)))
summary(ins$SVLEN)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 50.0     57.0     70.0    243.8    110.0 719427.0 


# Get AC
ins$matches <- gsub("AC=", "", apply(ins, 1, function(x) extract_matches(x$split, "AC=.*")))
ins$AC <- as.numeric(ins$matches)
stopifnot(!any(is.na(ins$AC)))
summary(ins$AC)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    0       1       3    1063      39   45347 


# Get predicted gene fields
ins$matches <- apply(ins, 1, function(x) extract_matches(x$split, "PREDICTED_.*=.*"))
summary(sapply(ins$matches, length))
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   1.000   1.000   1.002   1.000   3.000 

convert_to_named_pair <- function(x) {
  split_pair <- strsplit(x, "=")[[1]]  # Split the string by "="
  key <- split_pair[1]                 # Extract the key
  value <- split_pair[2]               # Extract the value
  return(setNames(value, key))         # Return a named value (list with the key)
} 
pairs <- sapply(ins$matches, function(x) unlist(lapply(x, convert_to_named_pair)))
predicted_fields <- unique(names(unlist(pairs)))
table(names(unlist(pairs)))
# PREDICTED_INTRONIC         PREDICTED_LOF PREDICTED_NEAREST_TSS    PREDICTED_PROMOTER         PREDICTED_UTR 
#              28050                   356                 42929                   662                   988
ins <- cbind(ins, ldply(pairs, rbind))
all((length(predicted_fields) - apply(ins, 1, function(x) sum(is.na(x)))) == (sapply(ins$matches, length))) # Check expected number of NAs

# Some commas in gene lists
v <- c("matches", predicted_fields)
View(ins[order(ins$PREDICTED_LOF), ..v])


# Get just predicted annotations
ins$matches <-  apply(ins, 1, function(x) extract_matches(x$split, "^PREDICTED_[^=]+$"))
ins$predicted <- sapply(ins$matches, "[", 1)
table(ins$predicted)
# PREDICTED_INTERGENIC 
#                43511
ins$split[5]






# Final filter to wanted SVLEN and AC
ins <- ins[ins$AC <= 10 & ins$AC > 0, ]














