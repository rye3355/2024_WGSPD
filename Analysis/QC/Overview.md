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
gsutil cp gs://2024-wgspd/qc/20240408_subset_sample_qc1.tsv files/
gsutil cp gs://2024-wgspd/qc/20240408_subset_hq_sample_qc1.tsv files/
```

Double check initial sample-level filtering based on read metrics 
```bash
Rscript 02_sample-filtering.R
```


Initial sample-level filtering based on read metrics
```bash
Rscript 02_sample-filtering.R
```


LD Prune
```bash

```


# More QC
You can also refer to specific scripts: [script1.sh](script1.sh)

```bash
bash script1.sh
```
