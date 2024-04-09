"""
First pass of hail sample qc
"""
import hail as hl

hl.init(default_reference = 'GRCh38',
                tmp_dir = "gs://wes-bipolar-tmp-4day/")


# Read original data MT
MT = "gs://2024-wgspd/qc/20240408_subset_initial-var-QC.mt"
mt = hl.read_matrix_table(MT)

# Run sample_qc
mt = hl.methods.sample_qc(mt, name = 'sample_qc1')

# Save 
mt.cols().select('sample_qc1').flatten().export(output='gs://2024-wgspd/qc/20240408_subset_sample_qc1.tsv')