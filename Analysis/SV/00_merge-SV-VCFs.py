"""
Merge VCFs
"""
import hail as hl

hl.init(default_reference = 'GRCh38',
                tmp_dir = "gs://wes-bipolar-tmp-4day/")


# Read original data MT
mt = hl.import_vcf("gs://2024-wgspd/WGSDP_Releasable_4.1_202405/annotated_by_pop.vcfs/gnomAD_SV_v3.release_4_1.WGSDP.sampleID_pcr_status.chr22.annotated.GD_Evi_fixed.sites.gz",
                   force_bgz=True)

print(mt.describe())
#print(mt.info.show())

s = mt.info.PREDICTED_UTR.collect()
print(s)