"""
hailctl dataproc start rye \
    --num-workers 10 \
    --num-secondary-workers 200 \
    --master-machine-type n1-highmem-16 \
    --worker-machine-type n1-highmem-16 \
    --packages gnomad,"git+https://github.com/broadinstitute/gnomad_qc.git@main" \
    --max-idle=15m
"""
import hail as hl
hl.init(gcs_requester_pays_configuration = 'wes-bipolar',
        tmp_dir = 'gs://wes-bipolar-tmp-4day',
        default_reference = 'GRCh38')


# Read subset
mt = hl.read_matrix_table("gs://gnomad-subsets-2024/gnomad-v3/202403/20240328_subset_dense-callstats.mt")
print(f"Starting: {mt.count()}")
# Starting: (588713326, 35527)

# Not done here, but after filtering out AS_lowqual and chrM: 437659765
# There are also some other qc steps done before the next filtering table by gnomad (see message from Qin)

# HT from Qin
v = hl.read_table("gs://gnomad-subsets-2024/gnomad-v3/gnomad.genomes.v4.0.final_filter.ht")
# 990067456

# Annotate in filters
mt = mt.annotate_rows(filters = v[mt.locus, mt.alleles].filters)
# Filter to passing
mt = mt.filter_rows(hl.len(mt.filters) == 0)
print(f"Passing: {mt.count()}")
# Passing: (348351706, 35527)

# Write
MT = "gs://gnomad-subsets-2024/gnomad-v3/202403/20240402_subset_passing-vars.mt"
print(f"Writing to: {MT}")
mt.write(MT)