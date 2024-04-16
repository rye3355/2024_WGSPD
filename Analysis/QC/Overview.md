# Quality Control

Filter out variants outside of LCRs, filter out low quality reads. 
```bash
hailctl dataproc start rye \
    --num-workers 5 \
    --packages gnomad \
    --autoscaling-policy=test-5-200 \
    --max-idle=15m
hailctl dataproc submit rye 00_variant-qc.py
```

Run Hail sample qc to generate basic metrics.
```bash
hailctl dataproc submit rye 01_sample-qc.py
gsutil cp gs://2024-wgspd/qc/20240408_subset_sample_qc1.tsv files/
```

Initial sample-level filtering based on read metrics
```bash
Rscript 02_sample-qc_MAD-4.R
```

or C++

```cpp
#include <iostream>

int main() {
    std::cout << "Hello World!";
    return 0;
}
```

# More QC
You can also refer to specific scripts: [script1.sh](script1.sh)

```bash
bash script1.sh
```
