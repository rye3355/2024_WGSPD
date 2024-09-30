library(data.table)
library(plyr)
d <- fread("/Users/rye/Projects/2024_WGSPD/Analysis/SV/WGSDP_Releasable_4.1_202405/joined.tsv.gz") # (2063922)
colnames(d) <- c("CHROM",	"POS",	"ID",	"REF",	"ALT",	"QUAL",	"FILTER",	"INFO")


del <- d[d$ALT == "<DEL>" & d$FILTER == "PASS",] # Filter to high quality DEL (598205)




#######################################################################################################
# Processing INFO field
#######################################################################################################
extract_matches <- function(row, pattern) {
  grep(pattern, row, value = TRUE)  # Searches and returns the values that match the pattern
}
del$split <- strsplit(del$INFO, ";")

# Double check SVTYPE is as expected
del$matches <- gsub("SVTYPE=", "", apply(del, 1, function(x) extract_matches(x$split, "SVTYPE=.*")))
stopifnot(all(del$matches == "DEL"))


# Get SVLEN 
del$matches <- gsub("SVLEN=", "", apply(del, 1, function(x) extract_matches(x$split, "SVLEN=.*")))
del$SVLEN <- as.numeric(del$matches)
stopifnot(!any(is.na(del$SVLEN)))
summary(del$SVLEN)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#   50      137      704     8621     4071 75976918 


# Get AC
del$matches <- gsub("AC=", "", apply(del, 1, function(x) extract_matches(x$split, "AC=.*")))
del$AC <- as.numeric(del$matches)
stopifnot(!any(is.na(del$AC)))
summary(del$AC)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    0       1       3    1063      39   45347 


# Get predicted gene fields
del$matches <- apply(del, 1, function(x) extract_matches(x$split, "PREDICTED_.*=.*"))
summary(sapply(del$matches, length))
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   1.000   1.000   1.002   1.000   3.000 

convert_to_named_pair <- function(x) {
  split_pair <- strsplit(x, "=")[[1]]  # Split the string by "="
  key <- split_pair[1]                 # Extract the key
  value <- split_pair[2]               # Extract the value
  return(setNames(value, key))         # Return a named value (list with the key)
} 
pairs <- sapply(del$matches, function(x) unlist(lapply(x, convert_to_named_pair)))
predicted_fields <- unique(names(unlist(pairs)))
table(names(unlist(pairs)))
# PREDICTED_INTRONIC         PREDICTED_LOF PREDICTED_NEAREST_TSS    PREDICTED_PROMOTER         PREDICTED_UTR 
#              28050                   356                 42929                   662                   988
del <- cbind(del, ldply(pairs, rbind))
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














