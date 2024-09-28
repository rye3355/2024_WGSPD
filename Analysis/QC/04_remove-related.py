"""
Remove related samples and checkpoint to final, qc'd subset
"""
import hail as hl

hl.init(default_reference = 'GRCh38',
                tmp_dir = "gs://wes-bipolar-tmp-4day/")

# Read original data MT
MT = "gs://2024-wgspd/qc/20240408_subset_initial-var-QC.mt"
mt = hl.read_matrix_table(MT)


# Slim to wanted samples
manifest = hl.import_table("gs://2024-wgspd/qc/20240905_WGSPD_final-qcd-manifest.tsv", key = "s")
mt = mt.filter_cols(hl.is_defined(manifest[mt.s]), keep = True)
print(f"Final kept samples: {mt.count()}")
# Final kept samples: (332778683, 29124)


mt.repartition(50000, shuffle = False).write("gs://2024-wgspd/qc/20240905_subset_final-qcd.mt", overwrite = True)

