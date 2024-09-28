# Structural Variation

Grab VCFs (from Xuefang)
```bash
gsutil -m cp -r gs://talkowski-sv-gnomad-v3-release/WGSDP_Releasable_4.1_202405/annotated_by_pop.vcfs/ gs://2024-wgspd/WGSDP_Releasable_4.1_202405/
gsutil -m cp -r gs://2024-wgspd/WGSDP_Releasable_4.1_202405/ .
```


Merge VCFs
```bash
git clone https://github.com/talkowski-lab/svtk.git
cd svtk
pip install -e .
svtk vcf2bed --no-header WGSDP_Releasable_4.1_202405/annotated_by_pop.vcfs/gnomAD_SV_v3.release_4_1.WGSDP.sampleID_pcr_status.chr18.annotated.GD_Evi_fixed.sites.gz test.bed
```