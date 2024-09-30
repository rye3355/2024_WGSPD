# Structural Variation

Grab VCFs (from Xuefang)
```bash
gsutil -m cp -r gs://talkowski-sv-gnomad-v3-release/WGSDP_Releasable_4.1_202405/annotated_by_pop.vcfs/ gs://2024-wgspd/WGSDP_Releasable_4.1_202405/
gsutil -m cp -r gs://2024-wgspd/WGSDP_Releasable_4.1_202405/ .
```


Merge VCFs
```bash
# Unzip
gunzip WGSDP_Releasable_4.1_202405/annotated_by_pop.vcfs/gnomAD_SV_v3.release_4_1.WGSDP.sampleID_pcr_status.chr*.gz 
# Get line of first occurence of #CHROM (296)
grep -n -m 1 '#CHROM' gnomAD_SV_v3.release_4_1.WGSDP.sampleID_pcr_status.chr1.annotated.GD_Evi_fixed.sites | sed  's/\([0-9]*\).*/\1/'
# Remove all but header

for i in $(seq 1 22);
do
    sed '1,296d' WGSDP_Releasable_4.1_202405/annotated_by_pop.vcfs/gnomAD_SV_v3.release_4_1.WGSDP.sampleID_pcr_status.chr"${i}".annotated.GD_Evi_fixed.sites > test/chr"${i}".sites
done
cat test/* > WGSDP_Releasable_4.1_202405/joined.tsv
rm -r test
gzip WGSDP_Releasable_4.1_202405/joined.tsv

gsutil cp WGSDP_Releasable_4.1_202405/joined.tsv.gz gs://2024-wgspd/WGSDP_Releasable_4.1_202405/
```