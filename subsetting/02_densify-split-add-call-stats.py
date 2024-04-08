"""
Doing it this way to use dataproc (batch corruption seems to be happening again)
Remove many alleles, split multi, 


hailctl dataproc start rye \
    --num-workers 10 \
    --num-secondary-workers 200 \
    --master-machine-type n1-highmem-16 \
    --worker-machine-type n1-highmem-16 \
    --packages gnomad,"git+https://github.com/broadinstitute/gnomad_qc.git@main" \
    --max-idle=15m
"""
import hail as hl

from gnomad.utils.annotations import annotate_adj


hl.init(gcs_requester_pays_configuration = 'wes-bipolar',
        tmp_dir = 'gs://wes-bipolar-tmp-4day',
        default_reference = 'GRCh38')


# Read subset VDS
vds = hl.vds.read_vds("gs://gnomad-subsets-2024/gnomad-v3/v3_filtered-samples.vds")
#print(vds.validate())

# Remove many alleles
print(f"Removing many alleles")
many_alleles = vds.variant_data.filter_rows(vds.variant_data.alleles.length() > 6).rows()
vds = hl.vds.filter_variants(vds, many_alleles, keep=False)

# Split multi
print(f"Splitting multi")
vds = hl.vds.split_multi(vds, filter_changed_loci = True)

# Densify
print(f"Densifying VDS")
mt = hl.vds.to_dense_mt(vds)
mt.describe()

"""
# Get info
info_ht = get_info().ht()# Note: AC and AC_raw are computed over all gnomAD samples
info_ht = info_ht.annotate(info=info_ht.info.drop('AC', 'AC_raw'))

mt = mt.drop(mt.gvcf_info) # Note: gvcf_info is sparse-specific
mt = mt.annotate_rows(info=info_ht[mt.row_key].info)

# Add callstats
print("Adding subset callstats")
mt = mt.annotate_rows(subset_callstats_raw=hl.agg.call_stats(mt.GT, mt.alleles))
mt = mt.annotate_rows(
        info=mt.info.annotate(
        AC_raw=mt.subset_callstats_raw.AC[1],
        AN_raw=mt.subset_callstats_raw.AN,
        AF_raw=hl.float32(mt.subset_callstats_raw.AF[1])))
mt = mt.drop('subset_callstats_raw')
mt = annotate_adj(mt)
mt = mt.annotate_rows(
        subset_callstats_adj=hl.agg.filter(mt.adj, hl.agg.call_stats(mt.GT, mt.alleles)))
mt = mt.annotate_rows(
        info=mt.info.annotate(
                AC=mt.subset_callstats_adj.AC[1],
                AN=mt.subset_callstats_adj.AN,
                AF=hl.float32(mt.subset_callstats_adj.AF[1])))
mt = mt.drop('subset_callstats_adj')
"""

# Add callstats
print("Adding subset callstats")
mt = mt.annotate_rows(subset_callstats_raw=hl.agg.call_stats(mt.GT, mt.alleles))
mt = mt.annotate_rows(
        AC_raw=mt.subset_callstats_raw.AC[1],
        AN_raw=mt.subset_callstats_raw.AN,
        AF_raw=hl.float32(mt.subset_callstats_raw.AF[1]))
mt = mt.drop('subset_callstats_raw')
mt = annotate_adj(mt)
mt = mt.annotate_rows(
        subset_callstats_adj=hl.agg.filter(mt.adj, hl.agg.call_stats(mt.GT, mt.alleles)))
mt = mt.annotate_rows(
    AC=mt.subset_callstats_adj.AC[1],
    AN=mt.subset_callstats_adj.AN,
    AF=hl.float32(mt.subset_callstats_adj.AF[1]))
mt = mt.drop('subset_callstats_adj')


# Write
MT = "gs://gnomad-subsets-2024/gnomad-v3/202403/20240328_subset_dense-callstats.mt"
print(f"Writing to: {MT}")
mt.write(MT)