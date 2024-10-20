# Quality Control

Filter out variants outside of LCRs, filter out low quality reads. 
```bash
hailctl dataproc start rye \
    --num-workers 5 \
    --packages gnomad,"git+https://github.com/broadinstitute/gnomad_qc.git@main" \
    --autoscaling-policy=test-5-200 \
    --max-idle=15m
hailctl dataproc submit rye 00_variant-qc.py
```

Run Hail sample qc to generate basic metrics on both 1) all samples and 2) samples passing gnomadv3 HQ sample QC
```bash
hailctl dataproc submit rye 01_sample-qc.py
gsutil cp gs://2024-wgspd/qc/20240903_subset_sample_qc1.tsv files/
gsutil cp gs://2024-wgspd/qc/20240903_subset_hq_sample_qc1.tsv files/
```

Double check initial sample-level filtering based on read metrics 
```bash
Rscript 02_sample-filtering.R
```


Initial sample-level filtering based on read metrics
```bash
Rscript 02_sample-filtering.R
```


LD Prune, pull resources using QoB from gnomad bucket
```python
import hail as hl
hl.init(gcs_requester_pays_configuration = 'wes-bipolar',
        tmp_dir = 'gs://wes-bipolar-tmp-4day',
        default_reference = 'GRCh38')

v2_qc_mt_liftover = hl.read_matrix_table("gs://gnomad/sample_qc/mt/gnomad.joint.high_callrate_common_biallelic_snps.pruned.grch38.mt/") # Read pruned MT, see from gnomad_qc.v2.resources.sample_qc import get_liftover_v2_qc_mt
v2_qc_mt_liftover.rows().write("gs://2024-wgspd/qc/20240423_gnomad.joint.high_callrate_common_biallelic_snps.pruned.grch38.ht", overwrite=True)
```

Using LD pruned variants from gnomadv3, run relatedness analysis (pc_relate)
```bash
hailctl dataproc submit rye 03_pc-relate.py
```


After relatedness output, get maximally independent (unrelated) set of samples
```python
import hail as hl
mt = hl.read_matrix_table("gs://2024-wgspd/qc/20240408_subset_initial-var-QC.mt") # Read in data
samples = mt.cols()
meta = hl.import_table("gs://2024-wgspd/files/gnomad_v3.1_subset-metadata.tsv", key="s") # Read in meta
samples = samples.annotate(high_quality = meta[samples.s].high_quality)
samples = samples.filter(samples.high_quality == "true") # Filter to high quality
manifest = hl.import_table("gs://2024-wgspd/files/2024_WGSPD_merged-manifest.tsv", key="s") # Read in manifest
samples = samples.annotate(case_con = manifest[samples.s].CASECON)
samples = samples.filter(hl.set(["CASE", "CTRL"]).contains(samples.case_con)) # Filter to CASE CON
samples = samples.annotate(is_case = samples.case_con == "CASE") # Define ordering (prefer CASE > CTRL)
def tie_breaker(l, r): # Define tiebreaker function 
    return hl.if_else(l.is_case & ~r.is_case, -1,
                      hl.if_else(~l.is_case & r.is_case, 1, 0))

pc_rel = hl.read_table("gs://2024-wgspd/qc/20240903_pc-relate_relatedness.ht") # Relatedness (4337 entries)
pairs = pc_rel.filter(pc_rel['kin'] > 0.05) # Slim to related, as used in gnomadQC (4337 pairs)
pairs_with_case = pairs.key_by(
    i=hl.struct(id=pairs.i, is_case=samples[pairs.i].is_case),
    j=hl.struct(id=pairs.j, is_case=samples[pairs.j].is_case)) # Annotate with ordering to use in tie_breaker

related_samples_to_remove = hl.maximal_independent_set(
   pairs_with_case.i, pairs_with_case.j, False, tie_breaker) # Find maximally indep set (samples to remove)

related_samples_to_remove.aggregate(hl.agg.counter(related_samples_to_remove.node.is_case)) 
# {False: 1328, True: 796}

result = manifest.filter(hl.is_defined(
    related_samples_to_remove.key_by(
       s = related_samples_to_remove.node.id.s)[manifest.s]))

result.export(output = "gs://2024-wgspd/qc/20240905_pc-relate_samples-to-remove.tsv")
```


Check concordance of related samples with gnomadv3 QC and write final manifest
```R
library(data.table)
system2("gsutil", "cp gs://2024-wgspd/qc/20240425_pc-relate_samples-to-remove.tsv files/")
to_remove = fread("files/20240905_pc-relate_samples-to-remove.tsv")
meta = fread("../../subsetting/files/gnomad_v3.1_subset-metadata.tsv")
manifest = fread("../../subsetting/2024_WGSPD_merged-manifest.tsv")

table(meta[meta$s %in% to_remove$s, "sample_filters.all_samples_related"], useNA = "always")
# 1458 / 2124 match up with gnomad
# FALSE  TRUE  <NA> 
#   666  1458     0 


table(to_remove$PRIMARY_DISEASE, useNA = "always")
#  BD  BD1  BD2 CASE CTRL  SCZ <NA> 
#  97   30    1   44 1328  624    0

table(to_remove$CASECON, useNA = "always")
# CASE CTRL <NA> 
# 796 1328    0

manifest = manifest[manifest$s %in% meta$s[meta$high_quality],] # Slim to high quality (32739)
manifest = manifest[manifest$CASECON %in% c("CASE", "CTRL"),] # Slim to case/control (31248)
manifest = manifest[!(manifest$s %in% to_remove$s)] # Remove relateds (29124)
manifest = merge(manifest, meta[, c("s", "population_inference.pop")], by = "s", all.x = T, all.y = F) # Get inferred populations
table(manifest$population_inference.pop)
#   afr   ami   amr   asj   eas   fin   mid   nfe   oth   sas 
# 11633    67  2279   377    40  4785     8  9329   574    32 
colnames(manifest) <- c("s", "SEX", "PRIMARY_DISEASE", "CASECON", "COHORT", "POP")
manifest$CHIP <- "WGS"

write.table(manifest, "20240905_WGSPD_final-qcd-manifest.tsv", sep = "\t",
            quote = F, col.names = T, row.names = F)

system2("gsutil", "cp 20240905_WGSPD_final-qcd-manifest.tsv gs://2024-wgspd/qc/")
```

Write out final QC'd MT
```bash
hailctl dataproc submit rye 04_remove-related.py
```

