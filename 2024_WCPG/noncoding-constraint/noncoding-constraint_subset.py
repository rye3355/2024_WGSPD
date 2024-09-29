"""
Subset down to noncoding constraint regions defined in bed file from Siwei:
https://drive.google.com/file/d/1EpyOsVloQ_hYzbb4VhmXCUYY9is_GjiF/view?usp=drive_link

Filtering to z > 4 (99th percentile)

gsutil cp /Users/rye/Projects/WGSPD/gnomad-v3-subset/20231005_noncoding-constraint/constraint_z_genome_1kb_filtered.ft17-4-9.nc.sorted.copy.bed.gz gs://gnomad-wgspd-subset/noncoding-constraint/
"""
import hail as hl
hl.init(default_reference = 'GRCh38',
                tmp_dir = "gs://wes-bipolar-tmp-4day/")


# Read in data MT
MT= 'gs://2024-wgspd/snv/coding/202240928_subset_post-qc_protein-coding.mt'
mt = hl.read_matrix_table(MT)





# Read in bed file with intervals
BED = 'gs://gnomad-wgspd-subset/noncoding-constraint/constraint_z_genome_1kb_filtered.ft17-4-9.nc.sorted.copy.bed.gz'
bed = hl.import_bed(BED, reference_genome = 'GRCh38', force = True)
bed = bed.annotate(target_int = hl.float(bed.target)) # Fix target 
#bed.count() # 1843559
bed_4 = bed.filter(bed.target_int > 4) # Filter to z > 4
#bed_4.count() # 19471

# Filter
mt_filtered = hl.filter_intervals(mt, bed_4.interval.collect(), keep = True)
mt_filtered.persist() # wes-bipolar-tmp-4day/persist_MatrixTable3DMIi2lUGh
#mt_filtered.count() # (1240621, 29124)

mt_filtered = hl.variant_qc(mt_filtered)
mt_1 = mt_filtered.filter_rows((mt_filtered.variant_qc.n_non_ref==1), keep = True)
mt_5 = mt_filtered.filter_rows((hl.min(mt_filtered.variant_qc.AC) <= 5), keep = True)
mt_10 = mt_filtered.filter_rows((hl.min(mt_filtered.variant_qc.AC) <= 10), keep = True)

mt_1 = mt_1.checkpoint("gs://wes-bipolar-tmp-4day/noncoding-constriant/target-int4_MAC1")
mt_5 = mt_5.checkpoint("gs://wes-bipolar-tmp-4day/noncoding-constriant/target-int4_MAC5")
mt_10 = mt_10.checkpoint("gs://wes-bipolar-tmp-4day/noncoding-constriant/target-int4_MAC10")


#{'CASE': 7927, 'CTRL': 21197}
def compute_rate_ratio(m: hl.MatrixTable):
    m = m.annotate_cols(total_obs = hl.agg.sum(m.GT.is_non_ref())) # Label samples as carriers (any count > 0)
    res = m.cols()
    a = res.aggregate(hl.agg.filter(res.case_con == "CASE", hl.agg.sum(res.total_obs))) # Count carriers
    b = res.aggregate(hl.agg.filter(res.case_con == "CTRL", hl.agg.sum(res.total_obs))) # Count non-carriers
    A = 7927
    B = 21197
    # {fmsb} rateratio
    rate_ratio = (a/A) / (b/B)
    norm_p = hl.qnorm(1 - (1-0.95)/2)
    c = (a - (A/(A+B)) * (a+b))/hl.sqrt((a+b) * (A/(A+B)) * (B/(A+B)))
    p_value = 2 * (1 - hl.pnorm(hl.abs(c)))
    se = hl.sqrt((1/a) + (1/b))
    ci_95_lower = rate_ratio * hl.exp(-norm_p*se)
    ci_95_upper = rate_ratio * hl.exp(norm_p*se)
    p_value = 2*(1-hl.pnorm(hl.abs(hl.log(rate_ratio) / se)))
    return [rate_ratio, hl.eval(ci_95_lower), hl.eval(ci_95_upper), hl.eval(p_value), hl.eval(se), 
            a, b, A, B]




compute_rate_ratio(mt_1)
# [0.9521359521943918, 0.9467602861217311, 0.957542141078527, 0.0, 0.0028887776402703193, 162500, 456373, 7927, 21197]
compute_rate_ratio(mt_5)
# [1.011664001356293, 1.0092853678162914, 1.0140482407414777, 0.0, 0.0012010309809038033, 955531, 2525655, 7927, 21197]
compute_rate_ratio(mt_10)
# [1.0194122256031068, 1.0173001355027165, 1.0215287007659153, 0.0, 0.001058192799652238, 1233490, 3235574, 7927, 21197]




#{'CASE': 7927, 'CTRL': 21197}
def compute_fisher(m: hl.MatrixTable):
    m = m.annotate_cols(carrier = hl.agg.any(m.GT.is_non_ref())) # Label samples as carriers (any count > 0)
    res = m.cols()
    case_carriers = res.aggregate(hl.agg.filter(res.case_con == "CASE", hl.agg.sum(res.carrier))) # Count carriers
    case_non_carriers =  7927 - case_carriers # Use dictionary to get non-carrier counts (faster than re-counting)
    control_carriers = res.aggregate(hl.agg.filter(res.case_con == "CTRL", hl.agg.sum(res.carrier))) 
    control_non_carriers = 21197 - control_carriers 

    res = hl.eval(hl.fisher_exact_test(case_carriers, case_non_carriers,
                                        control_carriers, control_non_carriers))
    return [hl.eval(res.p_value), hl.eval(res.odds_ratio), hl.eval(res.ci_95_lower), hl.eval(res.ci_95_upper),
            case_carriers, case_non_carriers, control_carriers, control_non_carriers]

compute_fisher(mt_1)
# [0.0010783373565330164, 2.592593878986912, 1.4070333356232791, 5.222051099755805, 7915, 12, 21114, 83]
compute_fisher(mt_5)
# [nan, nan, nan, nan, 7927, 0, 21197, 0]
compute_fisher(mt_10)
# [nan, nan, nan, nan, 7927, 0, 21197, 0]