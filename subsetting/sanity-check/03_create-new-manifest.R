library(data.table)
library(dplyr)

# Read in subsetting manifest
gnomad_meta <- fread("../files/gnomad_v3.1_subset-metadata.tsv")
subset_manifest <- fread("../files/scz_bp_samplelist_gnomad3.tsv")
subset_manifest <- subset_manifest[subset_manifest$sample_id%in%gnomad_meta$s] # Slim to successfully subsetted

# Read in original WGSPD samples names
old_samples <- fread("WGSPD-SEQID_fixed_non-missing.txt")
stopifnot(all(old_samples$s %in% subset_manifest$sample_id)) # None missing


# Read in original WGSPD manifest
old_manifest <- fread("DS-MANIFEST-WGSPD-WGS_FIXED-PHENOTYPES.tsv")
stopifnot(all(old_samples$s %in% old_manifest$COLLABORATOR_SAMPLE_ID)) # None missing
old_manifest <- old_manifest[old_manifest$COLLABORATOR_SAMPLE_ID %in% old_samples$s,]
stopifnot(all(old_manifest$COLLABORATOR_SAMPLE_ID %in% subset_manifest$sample_id)) # All old samples in new subset
rm(old_samples)


# Read in manifest from Caroline (not all samples here)
car_manifest <- fread("../files/SCHEMA_WGS_WorkingManifest.tsv")

sum(duplicated(car_manifest$subject_id)) # 8 duplicates for some reason
# View(car_manifest[duplicated(car_manifest$subject_id) | duplicated(car_manifest$subject_id, fromLast = T),])
# "11C121536"  "09C97227"   "10C101596"  "10C113795"  "11C126482"  "10C105852"  "MH0134754"  "8007542477"
# "G84381_11C121536"   "G95284_09C97227"    "G95284_10C101596"   "G95826_10C113795"   "G95831_11C126482"   "G95835_10C105852"   "G95837_MH0134754"   "RP-1365_8007542477"

car_manifest <- car_manifest[!duplicated(car_manifest$subject_id, fromLast = T) | duplicated(car_manifest$subject_id),]


stopifnot(all(old_manifest$COLLABORATOR_SAMPLE_ID %in% car_manifest$subject_id)) # 1025 old samples not in Caroline's new manifest
table(old_manifest$COHORT[!(old_manifest$COLLABORATOR_SAMPLE_ID %in% car_manifest$subject_id)], useNA = "always") # From these cohorts:
            #   0 SUPER - PALOTIE            <NA> 
            # 173             852               0 
table(old_manifest$PRIMARY_DISEASE[!(old_manifest$COLLABORATOR_SAMPLE_ID %in% car_manifest$subject_id)], useNA = "always") # With these phenotypes:
      #     Psychosis      <NA> 
      # 173       852         0 


# See how samples in Caroline's manifest lines up with old manifest

common <- merge(car_manifest[, c("subject_id", "sex", "disease_description", "affected_status", 
                                 "primary_disease")],
                old_manifest[, c("COLLABORATOR_SAMPLE_ID", "GENDER", 
                                     "AFFECTED_STATUS", "AFFECTED_STATUS_FIXED", 
                                     "PRIMARY_DISEASE", "PRIMARY_DISEASE_FIXED")],
                    by.x = "subject_id",
                    by.y = "COLLABORATOR_SAMPLE_ID")
table(common[, c("sex", "GENDER")], useNA = "always") # Reported sex lines up pretty well -- use Caroline's
table(common[, c("affected_status", "AFFECTED_STATUS")], useNA = "always") # Lines up pretty well, use primary_disease for classifications and fall back on PRIMARY_DISEASE if missing

table(common[, c("primary_disease", "PRIMARY_DISEASE_FIXED")], useNA = "always") # Lines up pretty well, use primary_disease for classifications and fall back on PRIMARY_DISEASE if missing



# Create new manifest
d <- data.table("s" = subset_manifest$sample_id)

## Deal with sex first
d <- merge(d, car_manifest[, c("subject_id", "sex")], 
           by.x = "s", by.y = "subject_id",
           all.x = T)
d <- merge(d, old_manifest[, c("COLLABORATOR_SAMPLE_ID", "GENDER")], 
           by.x = "s", by.y = "COLLABORATOR_SAMPLE_ID",
           all.x = T)
colnames(d)[2:3] <- c("sex_new", "sex_old")
d$SEX <- coalesce(d$sex_new, d$sex_old) # 

table(d$SEX, useNA = "always")
    #        #N/A  Female    Male Unknown    <NA> 
    # 173       5    7204    8891     237   19017 
d$SEX <- case_when(d$SEX == "Female" ~ "Female",
                   d$SEX == "Male" ~ "Male",
                   !(is.na(d$SEX)) ~ "Unknown")



## Deal with primary disease next
d <- merge(d, car_manifest[, c("subject_id", "primary_disease")], 
           by.x = "s", by.y = "subject_id",
           all.x = T)
colnames(d)[5] <- "primary_disease_new"
d$primary_disease_new_fixed <- case_when(d$primary_disease_new %in% c("Bipolar", "Bipolar disorder", "Bipolar Disorder",
                                                          "Bipolar DIsorder", "bipolar disorders", "Bipolar with Psychosis",
                                                          "BP", "BP Disorder", "BP Proband", "Mania with Psychosis") 
                                ~ "BD",
                                d$primary_disease_new %in% c("Bipolar I", "Bipolar I disorder", "Bipolar I Disorder",
                                                          "Bipolar I with psychosis", "Bipolar I, Mania/Depression, multiple episodes, W/OUT psychosis", 
                                                          "Bipolar I, Mania/Depression, multiple episodes, w/psychosis", "BPI", "Manic episode (Bipolar I)",
                                                          "Manic Episode Bipolar i", "Manic episode with psychosis (bipolar I)", 
                                                          "Manic episode with psychosis (Bipolar I)") 
                                ~ "BD1",
                                d$primary_disease_new %in% c("Bipolar 2 Disorder", "Bipolar II", "Bipolar II with psychosis") ~ "BD2",
                                d$primary_disease_new %in% c("Schizophrenia", "F20", "F20.0", "F20.00", "F20.01", "F20.02", "F20.03",
                                                          "F20.05", "F20.09", "F20.1", "F20.10", "F20.11", "F20.12", "F20.14", 
                                                          "F20.2", "F20.20", "F20.22", "F20.3", "F20.30", "F20.31", "F20.32", 
                                                          "F20.39", "F20.4", "F20.5", "F20.50", "F20.55", "F20.58", "F20.6",
                                                          "F20.60", "F20.8", "F20.83", "F20.9", "SAD",
                                                          "Schizoaffective disorder, bipolar type", "Schizoaffective disorder, depressed type",
                                                          "Schizoaffective disorder, depression", "Schizoaffective disorder, depressive type",
                                                          "Schizophrenia", "schizophrenia spectrum disorders", "Schizophrenia, hebephrenic",
                                                          "Schizophrenia, Hebephrenic", "Schizophrenia, paranoid", "Schizophrenia, Paranoid",
                                                          "Schizophrenia, Undiff", "Schizophrenia, undifferentiated", 
                                                          "Schizophrenia, Undifferentiated", "Schizophrenia,undifferentiated", "scz", "SCZ") 
                                ~ "SCZ",
                                d$primary_disease_new %in% c("control", "Control", "Healthy Control") 
                                ~ "CTRL",
                                d$primary_disease_new %in% c("Case", "Schizophrenia, Depression, Bipolar Disorder") 
                                ~ "CASE",
                                !is.na(d$primary_disease_new) ~ "OTHER")
table(d[, c("primary_disease_new", "primary_disease_new_fixed")], useNA = "always")
table(d[d$primary_disease_new_fixed == "OTHER", c("primary_disease_new", "primary_disease_new_fixed")], useNA = "always")
table(d$primary_disease_new_fixed, useNA = "always")
  #  BD   BD1   BD2  CASE  CTRL OTHER   SCZ  <NA> 
  # 690  1182     6   140  4022  3102  6343 20042 


stopifnot(all(is.na(d$primary_disease_new[!(d$s %in% car_manifest$subject_id)])))
stopifnot(all(is.na(d$primary_disease_new_fixed[!(d$s %in% car_manifest$subject_id)])))

# Compare to old manifest method
d <- merge(d, old_manifest[, c("COLLABORATOR_SAMPLE_ID", "PRIMARY_DISEASE")],
           by.x = "s", by.y = "COLLABORATOR_SAMPLE_ID",
           all.x = T)
colnames(d)[7] <- "primary_disease_old"
d$primary_disease_old_fixed <- case_when(d$primary_disease_old %in% c("Bipolar", "Bipolar disorder", "Bipolar Disorder", "Bipolar DIsorder",
                                                              "bipolar disorders", "Bipolar with Psychosis", "BIPOLMAN", "BP", "BP Disorder",
                                                              "BP Proband", "Mania with Psychosis") 
                                     ~ "BD",
                                     d$primary_disease_old %in% c("Bipolar I", "Bipolar I disorder", "Bipolar I Disorder", "Bipolar I with psychosis",
                                                              "Bipolar I, Mania/Depression, multiple episodes, W/OUT psychosis",
                                                              "Bipolar I, Mania/Depression, multiple episodes, w/psychosis",
                                                              "BPI", "Manic episode (Bipolar I)", "Manic Episode Bipolar i", 
                                                              "Manic episode with psychosis (bipolar I)", "Manic episode with psychosis (Bipolar I)") 
                                     ~ "BD1",
                                     d$primary_disease_old %in% c("Bipolar 2 Disorder", "Bipolar II", "Bipolar II with psychosis") 
                                     ~ "BD2",
                                     d$primary_disease_old %in% c("anypsychosis", "Psychosis", "SAD", "Schizoaffective disorder, bipolar type",
                                                              "Schizoaffective disorder, bipolar type\xca", "Schizoaffective disorder, depressed type",
                                                              "Schizoaffective disorder, depression", "Schizoaffective disorder, depressive type",
                                                              "Schizophrenia", "schizophrenia spectrum disorders", "Schizophrenia, hebephrenic",
                                                              "Schizophrenia, Hebephrenic", "Schizophrenia, Hebephrenic\xca", "Schizophrenia, paranoid",
                                                              "Schizophrenia, Paranoid", "Schizophrenia, Undiff", "Schizophrenia, undifferentiated",
                                                              "Schizophrenia, Undifferentiated", "Schizophrenia,undifferentiated", "scz", "SCZ") 
                                     ~ "SCZ",
                                     d$primary_disease_old %in% c("control", "Control", "CONTROL", "Healthy Control") 
                                     ~ "CTRL",
                                     d$primary_disease_old %in% c("Schizophrenia, Depression, Bipolar Disorder") 
                                     ~ "CASE",
                                     !is.na(d$primary_disease_old) ~ "OTHER") 
table(d[, c("primary_disease_old", "primary_disease_old_fixed")], useNA = "always")
table(d[d$primary_disease_old_fixed == "OTHER", c("primary_disease_old", "primary_disease_old_fixed")], useNA = "always")
table(d$primary_disease_old_fixed, useNA = "always")
  #  BD   BD1   BD2  CASE  CTRL OTHER   SCZ  <NA> 
  # 807  1193     6    70  4497   874  7032 21048 

table(d[, c("primary_disease_new_fixed", "primary_disease_old_fixed")], useNA = "always")

table(d[d$primary_disease_old_fixed == "CTRL", c("primary_disease_new_fixed", "primary_disease_old_fixed")])
table(d[d$primary_disease_old_fixed == "BD" & d$primary_disease_new_fixed == "OTHER", c("primary_disease_new_fixed", "primary_disease_old_fixed")])
table(d[d$primary_disease_old_fixed == "CTRL" & d$primary_disease_new_fixed == "OTHER", c("primary_disease_new_fixed", "primary_disease_old_fixed")])
table(d[d$primary_disease_old_fixed == "CTRL" & d$primary_disease_new_fixed == "SCZ", c("primary_disease_new_fixed", "primary_disease_old_fixed")])



# Use Caroline's new manifest as default, fall back on old manifest if NA
d$PRIMARY_DISEASE <- coalesce(d$primary_disease_new_fixed, d$primary_disease_old_fixed)

# Rescue SCZ from old subset with SCZ categorization
d$PRIMARY_DISEASE[is.na(d$PRIMARY_DISEASE) &
                          is.na(d$primary_disease_new_fixed) &
                          d$primary_disease_old_fixed == "SCZ"] <- "SCZ"

# Clean up remaining unknown samples -> CONTROLS
d$PRIMARY_DISEASE[is.na(d$PRIMARY_DISEASE) &
                          is.na(d$primary_disease_new_fixed) &
                          is.na(d$primary_disease_old_fixed)] <- "CTRL" 

table(d$PRIMARY_DISEASE, useNA = "always")
  # BD   BD1   BD2  CASE  CTRL OTHER   SCZ  <NA> 
  # 690  1182     6   140 23039  3275  7195     0 

# Final Case/Con
d$CASECON <- case_when(d$PRIMARY_DISEASE %in% c("BD", "BD1", "BD2", "CASE", "SCZ") ~ "CASE",
                       d$PRIMARY_DISEASE %in% c("CTRL") ~ "CTRL",
                       T ~ "OTHER")

write.table(d, "../2024_WGSPD_merged-manifest.tsv", sep = "\t",
            quote = F, col.names = T, row.names = F)


# SUMMARY of fields:
# s = sample id (maps to gnomad id)
# sex_new = reported sex from Caroline's new manifest, NAs for samples not mapped
# sex_old = reported sex from old manifest, NAs for samples not mapped
# SEX = coalesced reported sex first from Caroline's new manifest, then falling back on old manifest if NA. Otherwise, remain NA (use later as sanity check during sex check)
# primary_disease_new = reported primary disease from Caroline's new manifest
# primary_disease_new_fixed = cleaned up primary disease from Caroline's new manifest
# primary_disease_old = reported primary disease from old manifest
# primary_disease_old_fixed = cleaned up primary disease from old manifest
# PRIMARY_DISEASE = coalesced fixed primary disease from Caroline's new manifest, then falling back on old manifest if NA
# CASECON = final case/control classifications

