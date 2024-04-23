"""
First pass of hail sample qc before and after filtering to gnomad-defined HQ samples
"""
import hail as hl

hl.init(default_reference = 'GRCh38',
                tmp_dir = "gs://wes-bipolar-tmp-4day/")


# Read original data MT
MT = "gs://2024-wgspd/qc/20240408_subset_initial-var-QC.mt"
mt = hl.read_matrix_table(MT)


"""
# Run sample_qc
mt = hl.methods.sample_qc(mt, name = 'sample_qc1')

# Save 
mt.cols().select('sample_qc1').flatten().export(output='gs://2024-wgspd/qc/20240408_subset_sample_qc1.tsv')
"""


# Slim to high quality samples
print(f"Starting: {mt.count()}")
# Starting: (332778683, 35527)
meta = hl.import_table("gs://2024-wgspd/gnomad_v3.1_subset-metadata.tsv", key="s")
mt = mt.annotate_cols(high_quality = meta[mt.s].high_quality)
mt = mt.filter_cols(mt.high_quality == "true")
print(f"HQ Filtered: {mt.count()}")
# HQ filtered: (332778683, 32739)

# Slim to wanted case/control samples
manifest = hl.import_table("gs://2024-wgspd/2024_WGSPD_merged-manifest.tsv", key = "s")
mt = mt.annotate_cols(case_con = manifest[mt.s].CASECON)
mt = mt.filter_cols(hl.set(["CASE", "CTRL"]).contains(mt.case_con))
print(f"Case/Control filtered: {mt.count()}")
# Case/Control filtered: (332778683, 30659)

# Run sample_qc again
mt = hl.methods.sample_qc(mt, name = 'hq_sample_qc1')

# Save 
mt.cols().select('hq_sample_qc1').flatten().export(output='gs://2024-wgspd/qc/20240408_subset_hq_sample_qc1.tsv')
