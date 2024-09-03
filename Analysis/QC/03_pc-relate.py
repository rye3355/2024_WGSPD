"""
Using LD pruned variants from gnomadv3, run relatedness analysis (pc_relate)

Downsampling number of variants from pruned to 50% -- still too large

Need more memory workers
hailctl dataproc start rye \
    --num-workers 10 \
    --num-secondary-workers 10 \
    --max-idle=15m
"""
import hail as hl

hl.init(default_reference = 'GRCh38',
                tmp_dir = "gs://wes-bipolar-tmp-4day/",
                driver_cores=8, worker_memory='highmem')

# Read original data MT
MT = "gs://2024-wgspd/qc/20240408_subset_initial-var-QC.mt"
mt = hl.read_matrix_table(MT)


# Slim to high quality samples
print(f"Starting: {mt.count()}")
# Starting: (332778683, 35527)
meta = hl.import_table("gs://2024-wgspd/files/gnomad_v3.1_subset-metadata.tsv", key="s")
mt = mt.annotate_cols(high_quality = meta[mt.s].high_quality)
mt = mt.filter_cols(mt.high_quality == "true")
print(f"HQ Filtered: {mt.count()}")
# HQ Filtered: (332778683, 32739)

# Slim to wanted case/control samples
manifest = hl.import_table("gs://2024-wgspd/files/2024_WGSPD_merged-manifest.tsv", key = "s")
mt = mt.annotate_cols(case_con = manifest[mt.s].CASECON)
mt = mt.filter_cols(hl.set(["CASE", "CTRL"]).contains(mt.case_con))
print(f"Case/Control filtered: {mt.count()}")
# Case/Control filtered: (332778683, 31248)




# Slim to pruned variants from gnomad
gnomad_pruned_vars = hl.read_table("gs://2024-wgspd/qc/20240423_gnomad.joint.high_callrate_common_biallelic_snps.pruned.grch38.ht")
print(f"gnomad pruned: {gnomad_pruned_vars.count()}")
# gnomad pruned: 94148

pruned = mt.filter_rows(hl.is_defined(gnomad_pruned_vars[mt.locus, mt.alleles]))
pruned.write("gs://wes-bipolar-tmp-4day/20240903_pc-relate_pruned-checkpoint.mt", overwrite = True)
pruned = hl.read_matrix_table("gs://wes-bipolar-tmp-4day/20240903_pc-relate_pruned-checkpoint.mt")
print(f"Pruned subset: {pruned.count()}")
# Pruned subset: (90588, 30659)
pruned = pruned.repartition(500)
#pruned = pruned.sample_rows(0.5)


eig, scores, _ = hl.hwe_normalized_pca(
            pruned.GT, k=10, compute_loadings=False
        )
scores = scores.checkpoint("gs://wes-bipolar-tmp-4day/20240903_pc-relate_scores-checkpoint.ht", overwrite = True)
scores = hl.read_table("gs://wes-bipolar-tmp-4day/20240903_pc-relate_scores-checkpoint.ht")

print("Starting pc_relate")
relatedness_ht = hl.pc_relate(
            pruned.GT,
            min_individual_maf=0.01,
            scores_expr=scores[pruned.col_key].scores,
            block_size=4096,
            min_kinship=0.1,
            statistics="all",
        )

print("Starting write")
relatedness_ht.write("gs://2024-wgspd/qc/20240903_pc-relate_relatedness.ht", overwrite = True)



