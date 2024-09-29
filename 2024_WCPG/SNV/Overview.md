# Single Nucleotide Variant analysis


Using outputs from normal SNV analysis (2024_WGSPD/Analysis/SNV/Overview.md)
Picking up after 2024_WGSPD/Analysis/SNV/coding/01_vep-annotate.py
Using:
gs://2024-wgspd/snv/coding/202240928_subset_post-qc_protein-coding.mt
gs://2024-wgspd/snv/coding/20240928_subset_post-qc_protein-coding_VEP-annotated_original.ht


Generate counts in a similar way as initial round of analysis


hailctl dataproc start rye \
    --num-workers 5 \
    --autoscaling-policy=test-5-200 \
    --max-idle=15m

```bash
# pLoF
# writes to: gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/20240928_WGSPD_pLoF-MAC5-gnomadAC10
hailctl dataproc submit rye 20240928_generate-counts_pLoF-MAC5-gnomadAC10.py
# writes to: gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/20240928_WGSPD_pLoF-MAC5
hailctl dataproc submit rye 20240928_generate-counts_pLoF-MAC5.py


# Synonymous
# writes to: gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/20240928_WGSPD_synonymous-MAC5-gnomadAC10
hailctl dataproc submit rye 20240928_generate-counts_synonymous-MAC5-gnomadAC10.py
# writes to: gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/20240928_WGSPD_synonymous-MAC5
hailctl dataproc submit rye 20240928_generate-counts_synonymous-MAC5.py
```




Fisher exact tests for the above counts
```bash
# pLoF
python3 coding/03_fisher-test.py \
    --mt gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/20240928_WGSPD_pLoF-MAC5-gnomadAC10 \
    --gene_lists individual,gnomAD-constrained,SCHEMA,NDD,ASD,all \
    --out gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/ \
    --file_prefix 20240928_WGSPD_pLoF-MAC5-gnomadAC10 \
    --tmp gs://wes-bipolar-tmp-4day/20240928/

python3 coding/03_fisher-test.py \
    --mt gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/20240928_WGSPD_pLoF-MAC5\
    --gene_lists individual,gnomAD-constrained,SCHEMA,NDD,ASD,all \
    --out gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/ \
    --file_prefix 20240928_WGSPD_pLoF-MAC5 \
    --tmp gs://wes-bipolar-tmp-4day/20240928/
# {'CASE': 7927, 'CTRL': 21197}

# Synonymous
# writes to: 
python3 coding/03_fisher-test.py \
    --mt gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/20240928_WGSPD_synonymous-MAC5-gnomadAC10 \
    --gene_lists individual,gnomAD-constrained,SCHEMA,NDD,ASD,all \
    --out gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/ \
    --file_prefix 20240928_WGSPD_synonymous-MAC5-gnomadAC10 \
    --tmp gs://wes-bipolar-tmp-4day/20240928/

python3 coding/03_fisher-test.py \
    --mt gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/20240928_WGSPD_synonymous-MAC5 \
    --gene_lists individual,gnomAD-constrained,SCHEMA,NDD,ASD,all \
    --out gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/ \
    --file_prefix 20240928_WGSPD_synonymous-MAC5 \
    --tmp gs://wes-bipolar-tmp-4day/20240928/

# {'CASE': 7927, 'CTRL': 21197}
```



Fisher looks weird, so rate ratio instead
```bash
# pLoF
python3 coding/06_rate-ratio-test.py \
    --mt gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/20240928_WGSPD_pLoF-MAC5-gnomadAC10 \
    --gene_lists gnomAD-constrained,SCHEMA,NDD,ASD,all \
    --out gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/ \
    --file_prefix 20240929_WGSPD_pLoF-MAC5-gnomadAC10 \
    --tmp gs://wes-bipolar-tmp-4day/20240929/

python3 coding/06_rate-ratio-test.py \
    --mt gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/20240928_WGSPD_pLoF-MAC5 \
    --gene_lists gnomAD-constrained,SCHEMA,NDD,ASD,all \
    --out gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/ \
    --file_prefix 20240929_WGSPD_pLoF-MAC5 \
    --tmp gs://wes-bipolar-tmp-4day/20240929/


# synonymous
python3 coding/06_rate-ratio-test.py \
    --mt gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/20240928_WGSPD_synonymous-MAC5-gnomadAC10 \
    --gene_lists gnomAD-constrained,SCHEMA,NDD,ASD,all \
    --out gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/ \
    --file_prefix 20240929_WGSPD_synonymous-MAC5-gnomadAC10 \
    --tmp gs://wes-bipolar-tmp-4day/20240929/

python3 coding/06_rate-ratio-test.py \
    --mt gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/20240928_WGSPD_synonymous-MAC5 \
    --gene_lists gnomAD-constrained,SCHEMA,NDD,ASD,all \
    --out gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/ \
    --file_prefix 20240929_WGSPD_synonymous-MAC5 \
    --tmp gs://wes-bipolar-tmp-4day/20240929/
```







Still looks weird, try raising AC threshold
```bash
# pLoF
# writes to: gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/20240929_WGSPD_pLoF-MAC10
hailctl dataproc submit rye 20240929_generate-counts_pLoF-MAC10.py
# writes to: gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/20240929_WGSPD_pLoF-MAC10-gnomadAC10
hailctl dataproc submit rye 20240929_generate-counts_pLoF-MAC10-gnomadAC10.py

# Without gnomadAC10
hailctl dataproc submit rye coding/06_rate-ratio-test.py \
    --mt gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/20240929_WGSPD_pLoF-MAC10 \
    --gene_lists gnomAD-constrained,SCHEMA,NDD,ASD,all \
    --out gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/ \
    --file_prefix 20240929_WGSPD_pLoF-MAC10 \
    --tmp gs://wes-bipolar-tmp-4day/20240929/
hailctl dataproc submit rye coding/03_fisher-test.py \
    --mt gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/20240929_WGSPD_pLoF-MAC10 \
    --gene_lists gnomAD-constrained,SCHEMA,NDD,ASD,all \
    --out gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/ \
    --file_prefix 20240929_WGSPD_pLoF-MAC10 \
    --tmp gs://wes-bipolar-tmp-4day/20240929/

# With gnomadAC10
hailctl dataproc submit rye coding/06_rate-ratio-test.py \
    --mt gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/20240929_WGSPD_pLoF-MAC10-gnomadAC10 \
    --gene_lists gnomAD-constrained,SCHEMA,NDD,ASD,all \
    --out gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/ \
    --file_prefix 20240929_WGSPD_pLoF-MAC10-gnomadAC10 \
    --tmp gs://wes-bipolar-tmp-4day/20240929/
hailctl dataproc submit rye coding/03_fisher-test.py \
    --mt gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/20240929_WGSPD_pLoF-MAC10-gnomadAC10 \
    --gene_lists gnomAD-constrained,SCHEMA,NDD,ASD,all \
    --out gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/ \
    --file_prefix 20240929_WGSPD_pLoF-MAC10-gnomadAC10 \
    --tmp gs://wes-bipolar-tmp-4day/20240929/
```

# Synonymous
```bash
# writes to: gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/20240929_WGSPD_synonymous-MAC10
hailctl dataproc submit rye2 20240929_generate-counts_synonymous-MAC10.py
# writes to: gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/20240929_WGSPD_synonymous-MAC10-gnomadAC10
hailctl dataproc submit rye2 20240929_generate-counts_synonymous-MAC10-gnomadAC10.py

# Without gnomadAC10
hailctl dataproc submit rye2 coding/06_rate-ratio-test.py \
    --mt gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/20240929_WGSPD_synonymous-MAC10 \
    --gene_lists gnomAD-constrained,SCHEMA,NDD,ASD,all \
    --out gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/ \
    --file_prefix 20240929_WGSPD_synonymous-MAC10 \
    --tmp gs://wes-bipolar-tmp-4day/20240929/
hailctl dataproc submit rye2 coding/03_fisher-test.py \
    --mt gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/20240929_WGSPD_synonymous-MAC10 \
    --gene_lists gnomAD-constrained,SCHEMA,NDD,ASD,all \
    --out gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/ \
    --file_prefix 20240929_WGSPD_synonymous-MAC10 \
    --tmp gs://wes-bipolar-tmp-4day/20240929/

# With gnomadAC10
hailctl dataproc submit rye2 coding/06_rate-ratio-test.py \
    --mt gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/20240929_WGSPD_synonymous-MAC10-gnomadAC10 \
    --gene_lists gnomAD-constrained,SCHEMA,NDD,ASD,all \
    --out gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/ \
    --file_prefix 20240929_WGSPD_synonymous-MAC10-gnomadAC10 \
    --tmp gs://wes-bipolar-tmp-4day/20240929/
hailctl dataproc submit rye2 coding/03_fisher-test.py \
    --mt gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/20240929_WGSPD_synonymous-MAC10-gnomadAC10 \
    --gene_lists gnomAD-constrained,SCHEMA,NDD,ASD,all \
    --out gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/ \
    --file_prefix 20240929_WGSPD_synonymous-MAC10-gnomadAC10 \
    --tmp gs://wes-bipolar-tmp-4day/20240929/
```
