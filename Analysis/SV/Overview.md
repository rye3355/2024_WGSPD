# Structural Variation

Grab VCFs (from Xuefang)
For aggregateed info, copy over from gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/bipex_anno/
(copies to: gs://2024-wgspd/sv/20241104_WGSPD_Xuefang/)

Merge aggregated VCFs
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

Reformat merged version
```python
import hail as hl
hl.init(default_reference = "GRCh38")

# Read
ht = hl.import_table("gs://2024-wgspd/sv/20241104_WGSPD_Xuefang/joined.tsv.gz", force = True, no_header=True)
# Format/rename some fields
ht = ht.transmute(locus = hl.locus(ht.f0, hl.int(ht.f1)),
                  sv_id = ht.f2,
                  sv_type = ht.f4.replace("<", "").replace("[<>]", ""),
                  qual = ht.f5,
                  filters = ht.f6.split(";")
                  )
ht = ht.drop(ht.f3)
# Convert info field to a dictionary
ht = ht.transmute(info = hl.dict(
        hl.map(
            lambda d: hl.if_else(
                d.contains('='),
                (d.split('=')[0], d.split('=')[1]),
                (d, "")
            ),
            ht.f7.replace('"', '').split(';')
        )
    ))
# Also get keys
ht = ht.annotate(keys = ht.info.keys())
unique_elements = ht.aggregate(hl.agg.explode(lambda x: hl.agg.collect_as_set(x), ht.keys))
ht = ht.annotate_globals(all_keys = unique_elements)

# Write
ht = ht.repartition(100)
ht.write("gs://2024-wgspd/sv/20241104_WGSPD_Xuefang/gnomAD_SV_v3_joined.ht")
```



For GT-level info, see 20241104_WGSPD_Xuefang_get-files 
(copies to: gs://2024-wgspd/sv/20241104_WGSPD_sample-GT/)

Merge GT VCFs
```bash

hailctl dataproc start rye --num-workers 5 --autoscaling-policy=test-5-200 --max-idle=10m

# Merge GT VCFs into a MT
hailctl dataproc submit rye 00_merge-SV-VCFs.py
# Writes to gs://2024-wgspd/sv/20241104_WGSPD_sample-GT/gnomAD_SV_v3.release_4_1.1KGP_merged.mt
```
NOTE: not all samples from final-qcd manifest are present in MT, some are not releasable
 CASECON releasable unreleasable
   CASE       6866         1061
   CTRL      15969         5228







hailctl dataproc start rye --num-workers 5 --num-secondary-workers 10 --pkgs="ipython<8.22" --max-idle=10m
hailctl dataproc start rye --num-workers 5 --autoscaling-policy=test-5-200 --max-idle=10m