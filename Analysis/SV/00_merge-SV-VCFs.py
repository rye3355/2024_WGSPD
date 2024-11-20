"""
Merge VCFs
"""
import hail as hl

hl.init(default_reference = 'GRCh38',
                tmp_dir = "gs://wes-bipolar-tmp-4day/")

# Paths
gvcfs = [f"gs://2024-wgspd/sv/20241104_WGSPD_sample-GT/gnomAD_SV_v3.release_4_1.1KGP.chr{i}.annotated.vcf.gz" for i in list(range(1,23)) + ['X', 'Y']]
# Read
mts = [hl.import_vcf(gvcf, force_bgz = True) for gvcf in gvcfs]
"""
# Make sure all are compatible
r = mts[0]
assert(all([r.describe() == mt.describe() for mt in mts]))

# Checked sample names match up for only X and Y bc smaller...
assert(mts[-2].s.collect() == mts[-1].s.collect())
"""


# Union rows
merged = hl.MatrixTable.union_rows(*mts)
merged.write("gs://2024-wgspd/sv/20241104_WGSPD_sample-GT/gnomAD_SV_v3.release_4_1.1KGP_merged.mt", overwrite = True)


"""
counts = [mt.count() for mt in mts]
merge = hl.read_matrix_table("gs://2024-wgspd/sv/20241104_WGSPD_sample-GT/gnomAD_SV_v3.release_4_1.1KGP_merged.mt")
assert(sum([a[0] for a in counts]) == merge.count()[0])
"""