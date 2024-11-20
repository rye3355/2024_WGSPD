# Structural Variation

Grab VCFs (from Xuefang)
See 20241104_WGSPD_Xuefang_get-files
```bash
gsutil -m cp -r gs://2024-wgspd/sv/20241104_WGSPD_Xuefang .
```


Merge VCFs
```bash
# Unzip
gunzip gnomAD_SV_v3.release_4_1.bipex.chr*.gz 
# Get line of first occurence of #CHROM (790)
grep -n -m 1 '#CHROM' gnomAD_SV_v3.release_4_1.bipex.chr1.sites | sed  's/\([0-9]*\).*/\1/'# Remove all but header
# Remove all but header
mkdir no_header
for i in $(seq 1 22);
do
    sed '1,790d' gnomAD_SV_v3.release_4_1.bipex.chr"${i}".sites > no_header/chr"${i}".sites
done
cat no_header/* > joined.tsv
rm -r no_header
gzip joined.tsv

gsutil cp joined.tsv.gz gs://2024-wgspd/sv/20241104_WGSPD_Xuefang/
gsutil cp gs://2024-wgspd/sv/20241104_WGSPD_Xuefang/joined.tsv.gz gs://fc-54cd2a03-28fe-43ab-9142-1d265515b386/WGSPD/202411_SV/ 
```


hailctl dataproc start rye --num-workers 5 --autoscaling-policy=test-5-200 --max-idle=10m

Merge GT VCFs into a VDS
```python
import hail as hl

gvcfs = [f"gs://2024-wgspd/sv/20241104_WGSPD_sample-GT/gnomAD_SV_v3.release_4_1.1KGP.chr{i}.annotated.vcf.gz" for i in list(range(1,23)) + ['X', 'Y']]

combiner = hl.vds.new_combiner(
    output_path='gs://2024-wgspd/sv/20241104_WGSPD_sample-GT/20241104_WGSPD_sample-GT.vds',
    temp_path='gs://wes-bipolar-tmp-4day/',
    gvcf_paths=gvcfs,
    reference_genome="GRCh38",
    use_genome_default_intervals=True,
)

combiner.run()

```



