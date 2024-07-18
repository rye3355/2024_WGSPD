"""
Conduct MCH test on VEP-counts table (generated in 01_VEP-counts-export.py) 
Currently only supports stratification by ancestry x chip
"""

import hail as hl
import pandas as pd
import argparse
from datetime import date



gene_lists_info =   {"gnomAD-constrained":  {"Description": "Constrained genes as determined by gnomAD pLI > 0.9",
                                             "List path": "gs://2024-wgspd/files/gene-sets/gnomAD_pLI-constrained.tsv"},
                     "SCHEMA":              {"Description": "Top SCHEMA genes (P meta < 1.30e-04)",
                                             "List path": "gs://2024-wgspd/files/gene-sets/SCHEMA.tsv"},
                     "ASD":                 {"Description": "Top ASD genes (TADA FDR < 0.001)",
                                             "List path": "gs://2024-wgspd/files/gene-sets/ASD.tsv"},
                     "NDD":                 {"Description": "Top NDD genes (TADA FDR < 0.001)",
                                             "List path": "gs://2024-wgspd/files/gene-sets/NDD.tsv"},  
                    }

# Hail function doesn't report OR for some reason... so here's a custom function
def compute_MCH_OR(case_carriers, case_non_carriers,
                   control_carriers, control_non_carriers):
    def numerator_term(a, d, t):
        return a*d/t
    def denominator_term(b, c, t):
        return b*c/t
    n1 = hl.zip(case_carriers, case_non_carriers).map(lambda ab: ab[0] + ab[1])
    n2 = hl.zip(control_carriers, control_non_carriers).map(lambda cd: cd[0] + cd[1])
    t = hl.zip(n1, n2).map(lambda nn: nn[0] + nn[1])
    numerator = hl.sum(hl.zip(case_carriers, control_non_carriers, t).map(lambda tup: numerator_term(*tup)))
    denominator = hl.sum(hl.zip(case_non_carriers, control_carriers, t).map(lambda tup: denominator_term(*tup)))
    return numerator / denominator



def run_cmh_on_gene_lists(m: hl.MatrixTable, kept_groups: list[str], counts: hl.DictExpression):
    m = m.annotate_cols(carrier = hl.agg.any(m.agg > 0)) # Label samples as carriers (any count > 0)
    res = m.cols()
    res = res.group_by(res.group2).aggregate(n_carriers = hl.agg.sum(res.carrier)) # Count carriers by ancestry x kit x CASE/CTRL stratifications
    res = res.annotate(n_non_carriers = counts.get(res.group2) - res.n_carriers) # Use dictionary to get non-carrier counts (faster than re-counting)
    res = res.annotate(case_con = res.group2.split("_")[-1]) # Recover CASE/CTRL assignments
    case_carriers = res.aggregate(hl.agg.filter(res.case_con == "CASE", hl.agg.collect(res.n_carriers)))
    case_non_carriers = res.aggregate(hl.agg.filter(res.case_con == "CASE", hl.agg.collect(res.n_non_carriers)))
    control_carriers = res.aggregate(hl.agg.filter(res.case_con == "CTRL", hl.agg.collect(res.n_carriers)))
    control_non_carriers = res.aggregate(hl.agg.filter(res.case_con == "CTRL", hl.agg.collect(res.n_non_carriers))) # Collect counts into arrays (collect all ancestry x kit groups), slim to only rows now

    res = hl.eval(hl.cochran_mantel_haenszel_test(case_carriers, case_non_carriers,
                                                  control_carriers, control_non_carriers))
    
    OR = compute_MCH_OR(case_carriers, case_non_carriers,
                        control_carriers, control_non_carriers)

    return [OR, res.p_value, res.test_statistic]


def main(args):
    # Initialize
    hl.init(default_reference = 'GRCh38',
            tmp_dir = args.tmp)
        
    # Read in MT
    print(f"\nReading in MT ({args.mt})...\n")
    mt = hl.read_matrix_table(args.mt) 
    print(f"\nStarting (genes, samples): {mt.count()}\n")


    # Check annotation requirements
    if (args.annotate_casecon | args.annotate_pop | args.annotate_chip | args.filter_pass) and not args.manifest:
        print(f"No manifest provided to pull annotations from...")
        return
    
    # Annotate in fields from manifest
    # Requires key by "s" and specific fields "CASECON", "POP", and "CHIP" as needed
    if args.manifest:
        manifest = hl.import_table(args.manifest,
                                   delimiter = "\t", impute = True, key = "s")
        if args.filter_pass:
            mt = mt.filter_cols(hl.if_else(hl.is_defined(manifest[mt.s]),
                                           manifest[mt.s].FILTER == "PASS",
                                           False), 
                                keep = True)
        if args.annotate_casecon:
            mt = mt.annotate_cols(case_con = manifest[mt.s].CASECON)
        if args.annotate_pop:
            mt = mt.annotate_cols(pop = manifest[mt.s].POP)
        if args.annotate_chip:
            mt = mt.annotate_cols(chip = manifest[mt.s].CHIP)
        

        
    
    # Create population x chip stratification
    mt = mt.annotate_cols(group = mt.pop + "_" + mt.chip,
                          group2 = mt.pop + "_" + mt.chip + "_" + mt.case_con)
    groups = mt.aggregate_cols(hl.agg.collect_as_set(mt.group))
    counts = mt.aggregate_cols(hl.agg.counter(mt.group2))
    print(hl.eval(counts))

    for g in groups:
        if g + "_CASE" not in counts:
            counts[g + "_CASE" ] = 0
        if g + "_CTRL"  not in counts:
            counts[g + "_CTRL"] = 0
    


    # Exclude groups that are too small as desired
    excluded_groups = []
    if args.minimum_group_size or args.minimum_cases:
        if args.minimum_group_size:
            excluded_groups += [x for x in groups if ((counts[x + "_CASE"] + counts[x + "_CTRL"]) < args.minimum_group_size)]
        if args.minimum_cases:
            excluded_groups += [x for x in groups if ((counts[x + "_CASE"]) < args.minimum_cases)]   
        if args.minimum_controls:
            excluded_groups += [x for x in groups if ((counts[x + "_CTRL"]) < args.minimum_controls)]   

        print(f"\nGroups that are too small: {set(excluded_groups)}\n")
        print(f"\nTotal number of samples filtered out: {mt.filter_cols(hl.set(excluded_groups).contains(mt.group)).count()[1]}\n")
        mt = mt.filter_cols(hl.set(excluded_groups).contains(mt.group), keep = False)
        print(f"\nTotal number of samples after filtering: {mt.count()[1]}\n")

    kept_groups = [x for x in groups if x not in excluded_groups]

    print(counts)
    counts = hl.literal(counts) # Convert counts dictionary to hail format

    # Repartition if needed
    if mt.n_partitions() > 200:
        mt = mt.repartition(200, shuffle = False)

    # Persist for speed
    mt = mt.persist()


    ## Next, conduct individual gene-set burden analyses
    if args.gene_lists:
        gene_lists = args.gene_lists.replace(" ", "").split(",") # Parse command-line input


        # Running CMH on each gene individually
        # TODO: Make this more efficient with row aggregations?
        if "individual" in gene_lists:
            genes = mt.gene_symbol.collect() 
            print(f"\nAnalyzing {len(genes)} genes individually\n")
            m = mt
            res = m.group_cols_by(m.group2).aggregate(n_carriers = hl.agg.sum(m.agg > 0)).repartition(200, shuffle = False).persist() # Count carriers by ancestry x kit x CASE/CTRL stratifications
            res = res.annotate_entries(n_non_carriers = counts.get(res.group2) - res.n_carriers) # Use dictionary to get non-carrier counts (faster than re-counting)
            res = res.annotate_cols(case_con = res.group2.split("_")[-1]) # Recover CASE/CTRL assignments
            res_counts = res.annotate_rows(case_carriers = hl.agg.filter(res.case_con == "CASE", hl.agg.collect(res.n_carriers)), 
                                           case_non_carriers = hl.agg.filter(res.case_con == "CASE", hl.agg.collect(res.n_non_carriers)), 
                                           control_carriers = hl.agg.filter(res.case_con == "CTRL", hl.agg.collect(res.n_carriers)), 
                                           control_non_carriers = hl.agg.filter(res.case_con == "CTRL", hl.agg.collect(res.n_non_carriers))).rows() # Collect counts into arrays (collect all ancestry x kit groups), slim to only rows now
            res_counts = res_counts.annotate(RES = hl.cochran_mantel_haenszel_test(res_counts.case_carriers, res_counts.case_non_carriers,
                                                                                   res_counts.control_carriers, res_counts.control_non_carriers),
                                             OR = compute_MCH_OR(res_counts.case_carriers, res_counts.case_non_carriers,
                                                                 res_counts.control_carriers, res_counts.control_non_carriers)) # Annotate in MCH results and OR estimates
            res_counts = res_counts.flatten()

            #print(f"\nWriting result to: {args.out + args.file_prefix + '_case-control_CMH_individual.ht'}\n") 
            #res_counts = res_counts.checkpoint(args.out + args.file_prefix + '_case-control_CMH_individual.ht')

            print(f"\nWriting result to: {args.out + args.file_prefix + '_case-control_CMH_individual.tsv'}\n") 
            res_counts.export(args.out + args.file_prefix + '_case-control_CMH_individual.tsv', delimiter = "\t")
            
            
            gene_lists.remove("individual")


        results = {}
        for gene_list in gene_lists: # Iterate through every desired gene list analysis
            if gene_list == "individual":
                continue
            print(f"\nAnalyzing {gene_list}: {gene_lists_info[gene_list]['Description']}\n")
            print(f"Across {len(kept_groups)} groups: {kept_groups}\n")
            m = mt
            g = hl.import_table(gene_lists_info[gene_list]["List path"], 
                                delimiter = "\t", key = "gene_symbol", impute = True) # Read in list of gene symbols
            m = m.filter_rows(hl.is_defined(g[m.gene_symbol]), keep = True) # Filter down to genes of interest
            
            print(f"\nTotal number of genes: {m.rows().count()}\n")    
            results[gene_list] = [gene_list] + run_cmh_on_gene_lists(m, kept_groups, counts)

        # Only write if something to be written at this point
        if results:
            df = pd.DataFrame.from_dict(results, orient = "index",
                                        columns = ["gene_set", "OR", "p_value", "test_statistic"])
            df_ht = hl.Table.from_pandas(df, key = "gene_set")
            print(f"\nWriting result to: {args.out + args.file_prefix + '_case-control_CMH_' + '-'.join(gene_lists) + '.tsv'}\n") 
            df_ht.export(args.out + args.file_prefix + '_case-control_CMH_' + '-'.join(gene_lists) + '.tsv', delimiter='\t')



if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mt",
        help = "Path to starting VEP-counts MT (computed in 02_VEP-counts-export.py)",
        type = str,
        required = True
    )
    parser.add_argument(
        "--manifest",
        help = "Path to tab-delimited manifest keyed by samples ('s') with relevant annotations: 'FILTER', 'CASECON', 'POP', 'CHIP'",
        type = str,
        required = False
    )
    parser.add_argument(
        "--filter_pass",
        help = "Flag to filter to samples with FILTER == PASS in manifest",
        action = 'store_true'
    )
    parser.add_argument(
        "--annotate_casecon",
        help = "Flag to annotate in case control status for each sample",
        action = 'store_true'
    )
    parser.add_argument(
        "--annotate_pop",
        help = "Flag to annotate in population for each sample",
        action = 'store_true'
    )
    parser.add_argument(
        "--annotate_chip",
        help = "Flag to annotate in chip for each sample",
        action = 'store_true'
    )
    parser.add_argument(
        "--minimum_group_size",
        help = "Minimum ancestry-chip stratification group size (otherwise, drop it from analysis)",
        type = int,
        required = False
    )
    parser.add_argument(
        "--minimum_cases",
        help = "Minimum number of cases in ancestry-chip stratification group (otherwise, drop it from analysis)",
        type = int,
        required = False
    )
    parser.add_argument(
        "--minimum_controls",
        help = "Minimum number of controls in ancestry-chip stratification group (otherwise, drop it from analysis)",
        type = int,
        required = False
    )
    parser.add_argument(
        "--gene_lists",
        help = "Comma-separated gene-lists to individually filter to and analyze. Possible values: individual, gnomAD-constrained, SCHEMA, NDD, ASD.",
        type = str,
        required = False
    )
    parser.add_argument(
        "--out",
        help = "Path to output directory",
        type = str,
        required = True
    )
    parser.add_argument(
        "--file_prefix",
        help = "File output prefix",
        type = str,
        required = True
    )
    parser.add_argument(
        "--tmp",
        help = "Path to temp bucket",
        type = str,
        required = True
    )


    args = parser.parse_args()

    main(args)