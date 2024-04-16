# Filtering samples based on 4 MAD for the following metrics:
# c("n_insertion", "n_deletion", "n_snp", "r_het_hom_var", "r_ti_tv", "r_insertion_deletion")
# Also call rate >= 0.9

library(data.table)
library(stringr)

source("../plotting/pretty_plotting.r")

# Read in sample-qc data
d <- fread("files/20240408_subset_sample_qc1.tsv")

names(d) <- gsub("sample_qc1\\.", "", names(d))


# Filtering samples

## Call rate >= 0.9
d$CALL_RATE <- d$call_rate >= 0.9

## MAD threshold
MAD <- 4
## Fields to filter on
fields <- c("n_insertion", "n_deletion", "n_snp", "r_het_hom_var", "r_ti_tv", "r_insertion_deletion")
thresholds_low <- c()
thresholds_high <- c()

## Iterate
for (field in fields) {
  # Compute MAD and median
  MAD <- mad(d[, ..field])
  MED <- median(d[, field])
  low <- MED - MAD_THRESH * MAD
  high <- MED + MAD_THRESH * MAD
  thresholds_low <- c(thresholds_low, low)
  thresholds_high <- c(thresholds_high, high)
}
