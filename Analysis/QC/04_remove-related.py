"""
Remove related samples and checkpoint to final, qc'd subset
"""
import hail as hl

hl.init(default_reference = 'GRCh38',
                tmp_dir = "gs://wes-bipolar-tmp-4day/")

# Read original data MT
MT = "gs://2024-wgspd/qc/20240408_subset_initial-var-QC.mt"
mt = hl.read_matrix_table(MT)


# Slim to high quality samples
print(f"Starting: {mt.count()}")
# Starting: (332778683, 35527)
meta = hl.import_table("gs://2024-wgspd/files/gnomad_v3.1_subset-metadata.tsv", key="s")
mt = mt.filter_cols(meta[mt.s].high_quality == "true")
print(f"HQ Filtered: {mt.count()}")
# HQ filtered: (332778683, 32739)

# Slim to wanted case/control samples
manifest = hl.import_table("gs://2024-wgspd/files/2024_WGSPD_merged-manifest.tsv", key = "s")
mt = mt.annotate_cols(case_con = manifest[mt.s].CASECON)
mt = mt.filter_cols(hl.set(["CASE", "CTRL"]).contains(mt.case_con))
print(f"Case/Control filtered: {mt.count()}")
# Case/Control filtered: (332778683, 30659)

# Remove related samples
related = hl.import_table("gs://2024-wgspd/qc/20240425_pc-relate_samples-to-remove.tsv", key = "s")
mt = mt.filter_cols(hl.is_defined(related[mt.s]), keep = False)
print(f"Relateds filtered: {mt.count()}")
# Relateds filtered: (332778683, 30659)


mt.write("gs://2024-wgspd/qc/20240426_subset_final-qcd.mt", overwrite = True)

