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
                     "all":                 {"Description": "All provided genes"},
                    }

def run_cmh(mt: hl.MatrixTable, gene_list: str):
    print(f"\nAnalyzing {gene_list}: {gene_lists_info[gene_list]['Description']}...\n")
    
    g = hl.import_table(gene_lists_info[gene_list]["List path"], 
                        delimiter = "\t", key = "gene_symbol", impute = True) # Read in list of gene symbols
    
    m = mt.filter_rows(hl.is_defined(g[mt.gene_symbol]), keep = True) # Filter down to genes of interest
    print(f"\nTotal number of genes: {m.rows().count()}\n")
    m = m.annotate_cols(carrier = hl.agg.count_where(m.agg > 0))






    # TODO: Implement CMH stratified counts
    case_carriers = m.aggregate_cols(hl.agg.filter(m.case_con == "CASE", hl.agg.sum(m.carrier)))
    control_carriers = m.aggregate_cols(hl.agg.filter(m.case_con == "CTRL", hl.agg.sum(m.carrier)))
    case_non_carriers = m.aggregate_cols(hl.agg.filter(m.case_con == "CASE", hl.agg.sum(m.carrier == 0)))
    control_non_carriers = m.aggregate_cols(hl.agg.filter(m.case_con == "CTRL", hl.agg.sum(m.carrier == 0)))







    res = hl.eval(hl.cochran_mantel_haenszel_test(case_carriers, case_non_carriers,
                                                  control_carriers, control_non_carriers))

    return [gene_list, res.p_value, res.test_statistic]


def main(args):
    # Initialize
    hl.init(default_reference = 'GRCh38',
            tmp_dir = args.tmp)
        
    # Read in MT
    print(f"\nReading in MT ({args.mt})...\n")
    mt = hl.read_matrix_table(args.mt) 
    print(f"\nStarting (genes, samples): {mt.count()}\n")


    # Check annotation requirements
    if (args.annotate_casecon | args.annotate_pop | args.annotate_chip) and not args.manifest:
        print(f"No manifest provided to pull annotations from...")
        return
    
    # Annotate in fields from manifest
    # Requires key by "s" and specific fields "CASECON", "POP", and "CHIP" as needed
    if args.manifest:
        manifest = hl.import_table(args.manifest,
                                   delimiter = "\t", impute = True, key = "s")
        
        if args.annotate_casecon:
            mt = mt.annotate_cols(case_con = manifest[mt.s].CASECON)
        if args.annotate_pop:
            mt = mt.annotate_cols(pop = manifest[mt.s].POP)
        if args.annotate_chip:
            mt = mt.annotate_cols(chip = manifest[mt.s].CHIP)
        
    
    # Create population x chip stratification
    mt = mt.annotate_cols(group = mt.pop + "_" + mt.chip)
    counts = mt.aggregate_cols(hl.agg.counter(mt.group))

    # Exclude groups that are too small as desired
    excluded_groups = []
    if args.minimum_group_size | args.minimum_cases:
        if args.minimum_group_size:
            excluded_groups += [x for x in counts if counts[x] < args.minimum_group_size]
        if args.minimum_cases:
            case_counts = mt.aggregate_cols(hl.agg.filter(mt.case_con == "CASE", hl.agg.counter(mt.group)))
            excluded_groups += [x for x in case_counts if case_counts[x] < args.minimum_cases]   
        
        print(f"\nGroups that are too small: {set(excluded_groups)}\n")
        print(f"\nTotal number of samples filtered out: {mt.filter_cols(hl.set(excluded_groups).contains(mt.group)).count()[1]}\n")
        mt = mt.filter_cols(hl.set(excluded_groups).contains(mt.group), keep = False)
        print(f"\nTotal number of samples after filtering: {mt.count()[1]}\n")


    # Repartition if needed
    if mt.n_partitions() > 200:
        mt = mt.repartition(200, shuffle = False)

    # Persist for speed
    mt = mt.persist()



    ## Next, conduct individual gene-set burden analyses
    if args.gene_lists:
        results = {}
        gene_lists = args.gene_lists.replace(" ", "").split(",") # Parse command-line input
        for gene_list in gene_lists: # Iterate through every desired gene list analysis
            results[gene_list] = run_cmh(mt, gene_list)
            

        df = pd.DataFrame.from_dict(results, orient = "index",
                                    columns = ["gene_set", "p_value", "test_statistic"])
        df_ht = hl.Table.from_pandas(df, key = "gene_set")

        print(f"\nWriting result to: {args.out + args.file_prefix + '_case-control_Fisher-exact_' + '-'.join(gene_lists) + '.tsv'}\n") 
        df_ht.export(args.out + args.file_prefix + '_case-control_Fisher-exact_' + '-'.join(gene_lists) + '.tsv', delimiter='\t')

            


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
        help = "Path to tab-delimited manifest keyed by samples ('s') with relevant annotations: 'CASECON', 'POP', 'CHIP'",
        type = str,
        required = False
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
        "--gene_lists",
        help = "Comma-separated gene-lists to individually filter to and analyze. Possible values: all, gnomAD-constrained, SCHEMA, NDD, ASD.",
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