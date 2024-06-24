# Single Nucleotide Variant analysis


First, update coding SNV results

Subset to coding regions
```bash
hailctl dataproc start rye \
    --num-workers 5 \
    --packages gnomad \
    --autoscaling-policy=test-5-200 \
    --vep GRCh38 \
    --max-idle=15m

hailctl dataproc submit rye coding/00_filter-to-coding.py
```

```bash

```





Using additional annotatations from RGC (https://rgc-research.regeneron.com/me/resources)