# Quality Control

## Subsetting

First, subset from gnomadv3 using QoB using sample list (copy in *files/scz_bp_samplelist_gnomad3_FINAL.tsv*). 
```bash
python3 01_subset-v3-vds.py
```

Next, on dataproc, remove many alleles, split multi, and add call stats.
```bash
python3 02_densify-split-add-call-stats.py
```

Finally, on dataproc, filter to passing variants as determined by gnomadv4 (really v3 dataset, but analysis done under a different name) initial qc
```bash
python3 03_filter-passing-vars.py
```

# Missing samples (907 / 36434)
Used sample manifest with mappings to gnomadv3 sample IDs in *files/*
Note that not all 36434 samples from the manifest are subset out from gnomadv3. 35527 expected samples were pulled out (426 with no mapping to gnomadv3, 481 removed from gnomadv3 hard filtering; 907 total).
See *00_check-samples.ipynb* for more details



# Create new manifest with subset metadata
Use gnomadv3 sample metadata HT to create new, merged manifest.
```bash
python3 04_get-sample-info.py
gsutil cp gs://2024-wgspd/gnomad_v3.1_subset-metadata.tsv files/
#Rscript 05_check-meta.R
```
