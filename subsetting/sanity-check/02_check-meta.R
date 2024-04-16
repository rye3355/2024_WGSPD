library(data.table)

# Read in gnomad metadata
gnomad_meta <- fread("../files/gnomad_v3.1_subset-metadata.tsv")

# Read in subsetting manifest
subset_manifest <- fread("../files/scz_bp_samplelist_gnomad3.tsv")
stopifnot(all(gnomad_meta$s %in% subset_manifest$sample_id)) # Subsetting worked
not_subset <- subset_manifest[!(subset_manifest$sample_id %in% gnomad_meta$s),]
subset_manifest <- subset_manifest[subset_manifest$sample_id %in% gnomad_meta$s,]
table(not_subset$research_project) # 1174 samples not subset out


# Read in original WGSPD samples names
old_samples <- fread("WGSPD-SEQID_fixed_non-missing.txt")
stopifnot(all(old_samples$s %in% subset_manifest$sample_id)) # None missing


# Read in original WGSPD manifest
old_manifest <- fread("DS-MANIFEST-WGSPD-WGS_FIXED-PHENOTYPES.tsv")
stopifnot(all(old_samples$s %in% old_manifest$COLLABORATOR_SAMPLE_ID)) # None missing
old_manifest <- old_manifest[old_manifest$COLLABORATOR_SAMPLE_ID %in% old_samples$s,]
stopifnot(all(old_manifest$COLLABORATOR_SAMPLE_ID %in% subset_manifest$sample_id)) # All old samples in new subset
 


# Read in manifest from Caroline (not all samples here)
car_manifest <- fread("../files/SCHEMA_WGS_WorkingManifest.tsv")

sum(duplicated(car_manifest$subject_id)) # 8 duplicates for some reason
# View(car_manifest[duplicated(car_manifest$subject_id) | duplicated(car_manifest$subject_id, fromLast = T),])
# "11C121536"  "09C97227"   "10C101596"  "10C113795"  "11C126482"  "10C105852"  "MH0134754"  "8007542477"
# "G84381_11C121536"   "G95284_09C97227"    "G95284_10C101596"   "G95826_10C113795"   "G95831_11C126482"   "G95835_10C105852"   "G95837_MH0134754"   "RP-1365_8007542477"

car_manifest <- car_manifest[!duplicated(car_manifest$subject_id, fromLast = T) | duplicated(car_manifest$subject_id),]


stopifnot(all(old_manifest$COLLABORATOR_SAMPLE_ID %in% car_manifest$subject_id)) # 1025 old samples not in Caroline's new manifest
table(old_manifest$COHORT[!(old_manifest$COLLABORATOR_SAMPLE_ID %in% car_manifest$subject_id)], useNA = "always") # From these cohorts:
# 0 SUPER - PALOTIE            <NA> 
#   173             852               0 
table(old_manifest$PRIMARY_DISEASE[!(old_manifest$COLLABORATOR_SAMPLE_ID %in% car_manifest$subject_id)], useNA = "always") # With these phenotypes:
#       Psychosis      <NA> 
#   173       852         0 


# See how samples in Caroline's manifest lines up with old manifest

common <- merge(car_manifest[, c("subject_id", "sex", "disease_description", "affected_status", 
                                 "primary_disease")],
                old_manifest[, c("COLLABORATOR_SAMPLE_ID", "GENDER", 
                                     "AFFECTED_STATUS", "AFFECTED_STATUS_FIXED", 
                                     "PRIMARY_DISEASE", "PRIMARY_DISEASE_FIXED")],
                    by.x = "subject_id",
                    by.y = "COLLABORATOR_SAMPLE_ID")
table(common[, c("sex", "GENDER")], useNA = "always") # Reported sex lines up pretty well -- use Caroline's

table(common[, c("disease_description", "AFFECTED_STATUS")], useNA = "always") # Don't use disease_description
table(common[, c("affected_status", "AFFECTED_STATUS")], useNA = "always") # 233 disagreeing case/control, use old manifest for full?
table(common[, c("affected_status", "AFFECTED_STATUS_FIXED")], useNA = "always") # Don't use affected_status

table(common[, c("primary_disease", "PRIMARY_DISEASE_FIXED")], useNA = "always") # Lines up pretty well, use primary_disease for classifications and fall back on PRIMARY_DISEASE if missing


