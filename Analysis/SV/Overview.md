# Structural Variation

Grab VCFs (from Xuefang)
```bash
gsutil -m cp -r gs://talkowski-sv-gnomad-v3-release/WGSDP_Releasable_4.1_202405/annotated_by_pop.vcfs/ gs://2024-wgspd/WGSDP_Releasable_4.1_202405/
```


Merge VCFs
```bash
git clone https://github.com/talkowski-lab/svtk.git
cd svtk
pip install -e .

hailctl dataproc start rye \
    --num-workers 5 \
    --packages gnomad,"git+https://github.com/broadinstitute/gnomad_qc.git@main" \
    --autoscaling-policy=test-5-200 \
    --max-idle=15m
hailctl dataproc submit rye 00_variant-qc.py
```