"""
Filter to variants outside of LCRs.

Filter out low quality reads satisfying:
 - Defined and
 - If homref: GQ < 20 | DP < 10
 - If het:  (AD[0] + AD[1]) / DP) < 0.8 | 
            AD[1] / DP < 0.2 |
            PL[0] < 20 |
            DP < 10
 - If homvar:   (AD[1] / DP) < 0.8 |
                PL[0] < 20 |
                DP < 10 

Write to gs://2024-wgspd/qc/20240408_subset_initial-var-QC.mt

hailctl dataproc start rye \
    --num-workers 5 \
    --packages gnomad \
    --autoscaling-policy=test-5-200 \
    --max-idle=15m
"""
import hail as hl

hl.init(default_reference = 'GRCh38',
                tmp_dir = "gs://wes-bipolar-tmp-4day/")


# Read original data MT
MT = "gs://gnomad-subsets-2024/gnomad-v3/202403/20240402_subset_passing-vars.mt"
mt = hl.read_matrix_table(MT)
print(f"Original data count: {mt.count()}")
# Original data count: (348351706, 35527)

# LCRs
## Read in LCRs
LCRs = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/inputs/LCR-hs38.bed'
LCR_intervals = hl.import_locus_intervals(LCRs, reference_genome = 'GRCh38')

## Filter out variants in LCRs
mt = mt.filter_rows(~hl.is_defined(LCR_intervals[mt.locus]))
print(f"Number of variants outside of LCRs: {mt.count()}")
# Number of variants outside of LCRs: (332778683, 35527)


# Filter out low quality reads
mt = mt.filter_entries(
    hl.is_defined(mt.GT) &
    (
        (mt.GT.is_hom_ref() & 
            (
                # ((mt.AD[0] / mt.DP) < 0.8) | # Has to be removed because allele depth no longer defined for hom ref calls.
                (mt.GQ < 20) |
                (mt.DP < 10)
        	)
        ) |
        (mt.GT.is_het() & 
        	( 
                (((mt.AD[0] + mt.AD[1]) / mt.DP) < 0.8) | 
                ((mt.AD[1] / mt.DP) < 0.2) | 
                (mt.PL[0] < 20) |
                (mt.DP < 10)
        	)
        ) |
        (mt.GT.is_hom_var() & 
        	(
                ((mt.AD[1] / mt.DP) < 0.8) |
                (mt.PL[0] < 20) |
                (mt.DP < 10)
        	)
        )
    ),
    keep = False
)

# Write
mt.write("gs://2024-wgspd/qc/20240408_subset_initial-var-QC.mt", overwrite = True)