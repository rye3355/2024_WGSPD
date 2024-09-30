import hail as hl
MANIFEST = 'gs://2024-wgspd/files/20240905_WGSPD_final-qcd-manifest.tsv'
manifest = hl.import_table(MANIFEST, delimiter='\t',
                          key = "s", impute = True)

#mt_1 = hl.read_matrix_table("gs://wes-bipolar-tmp-4day/noncoding-constriant/target-int4_MAC1")
mt_5 = hl.read_matrix_table("gs://wes-bipolar-tmp-4day/noncoding-constriant/target-int4_MAC5")
mt_10 = hl.read_matrix_table("gs://wes-bipolar-tmp-4day/noncoding-constriant/target-int4_MAC10")

mt_5 = mt_5.annotate_cols(case_con = manifest[mt_5.s].CASECON)
mt_10 = mt_10.annotate_cols(case_con = manifest[mt_10.s].CASECON)
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


print(compute_rate_ratio(mt_5))
print(compute_fisher(mt_5))

print(compute_rate_ratio(mt_10))
print(compute_fisher(mt_10))
