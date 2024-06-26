# Single Nucleotide Variant analysis


First, update coding SNV results

Subset to coding regions and write in some helpful annotations.
Coding regions defined Gencode-v46 (https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz)
MPC, gnomAD-nonPsych, DiscovEHR, liftover OS (see below), AM, score_ml, MisFit_S, MisFit_D 
```bash
hailctl dataproc start rye \
    --num-workers 5 \
    --packages gnomad \
    --autoscaling-policy=test-5-200 \
    --vep GRCh38 \
    --requester-pays-allow-all \
    --max-idle=15m

hailctl dataproc submit rye coding/00_filter-to-coding.py
```

Lifting over OS annotation GRCh37 -> GRCh38
```bash
hailctl dataproc submit rye 20240617_convert-OS.py
```


VEP annotation and processing consequence categories
Using the same consequence category annotation system as initial analysis (for now)
{'damaging_missense': 354982, 'non_coding': 131179272, 'other_missense': 1592839, 'pLoF': 137917, 'synonymous': 882295, None: 10212733}
```bash
hailctl dataproc submit rye coding/01_vep-annotate.py
```


TODO: try different consequence category annotation methods



Generate HT for convenience from additional VCF annotatations from RGC (https://rgc-research.regeneron.com/me/resources)
```python
import hail as hl
hl.init(default_reference = 'GRCh38',
        tmp_dir = 'gs://wes-bipolar-tmp-4day/')

recode = {f"{i}":f"chr{i}" for i in (list(range(1, 23)) + ['X', 'Y'])} # Since improperly formatted
mt = hl.import_vcf('gs://bipex2/annotations/RGC/rgc_me_variant_frequencies_chr*_20231004.vcf.gz', force = True, reference_genome='GRCh38', contig_recoding=recode) # Read in
ht_sites = mt.rows() # Sites only, so take only rows
ht_sites = ht_sites.transmute(**ht_sites.info)
ht_sites.n_partitions()
ht_sites = ht_sites.repartition(50)
ht_sites.write('gs://bipex2/annotations/RGC/rgc_me_variant_frequencies.ht')
```




Generate counts in a similar way as initial round of analysis
```bash
# pLoF
hailctl dataproc submit rye coding/02_VEP-counts-export.py \
    --mt gs://2024-wgspd/snv/coding/202240618_subset_post-qc_protein-coding.mt \
    --vep_ht gs://2024-wgspd/snv/coding/20240625_subset_post-qc_protein-coding_VEP-annotated_original.ht \
    --gnomAD_AC_thresh 10 \
    --cons_cat pLoF \
    --ac 5 \
    --non_gnomAD_psych \
    --out gs://2024-wgspd/snv/coding/outputs/original-run-annotations/ \
    --file_prefix 20240626_BDSCZ \
    --tmp gs://wes-bipolar-tmp-4day/20240626/



#### TO RUN
# Missense
hailctl dataproc submit rye coding/02_VEP-counts-export.py \
    --mt gs://2024-wgspd/snv/coding/202240618_subset_post-qc_protein-coding.mt \
    --vep_ht gs://2024-wgspd/snv/coding/20240625_subset_post-qc_protein-coding_VEP-annotated_original.ht \
    --gnomAD_AC_thresh 10 \
    --cons_cat other_missense,damaging_missense \
    --ac 5 \
    --non_gnomAD_psych \
    --out gs://2024-wgspd/snv/coding/outputs/original-run-annotations/ \
    --file_prefix 20240626_BDSCZ \
    --tmp gs://wes-bipolar-tmp-4day/20240626/

# Synonymous
hailctl dataproc submit rye coding/02_VEP-counts-export.py \
    --mt gs://2024-wgspd/snv/coding/202240618_subset_post-qc_protein-coding.mt \
    --vep_ht gs://2024-wgspd/snv/coding/20240625_subset_post-qc_protein-coding_VEP-annotated_original.ht \
    --gnomAD_AC_thresh 10 \
    --cons_cat synonymous \
    --ac 5 \
    --non_gnomAD_psych \
    --out gs://2024-wgspd/snv/coding/outputs/original-run-annotations/ \
    --file_prefix 20240626_BDSCZ \
    --tmp gs://wes-bipolar-tmp-4day/20240626/
```


Exact tests for the above counts
```bash
# pLoF
hailctl dataproc submit rye coding/03_fisher-test.py \
    --mt gs://2024-wgspd/snv/coding/outputs/original-run-annotations/20240626_BDSCZ_gnomadAC10_pLoF_AC5_counts.mt \
    --individual \
    --gene_lists gnomAD-constrained,SCHEMA,NDD,ASD \
    --out gs://2024-wgspd/snv/coding/outputs/original-run-annotations/ \
    --file_prefix 20240626_BDSCZ_gnomadAC10_pLoF_AC5_non-gnomAD-psych \
    --tmp gs://wes-bipolar-tmp-4day/20240626/
```


hailctl dataproc submit rye coding/03_fisher-test.py \
    --mt gs://2024-wgspd/snv/coding/outputs/original-run-annotations/20240626_BDSCZ_gnomadAC10_pLoF_AC5_counts.mt \
    --manifest gs://2024-wgspd/files/20240523_WGSPD_final-qcd-manifest.tsv \
    --annotate_pop \
    --annotate_chip \
    --minimum_group_size 200 \
    --minimum_cases 50 \
    --gene_lists gnomAD-constrained,SCHEMA,NDD,ASD \
    --out gs://2024-wgspd/snv/coding/outputs/original-run-annotations/ \
    --file_prefix 20240626_BDSCZ_gnomadAC10_pLoF_AC5_non-gnomAD-psych \
    --tmp gs://wes-bipolar-tmp-4day/20240626/

hailctl dataproc start rye \
    --num-workers 5 \
    --packages gnomad \
    --autoscaling-policy=test-5-200 \
    --requester-pays-allow-all \
    --max-idle=30m