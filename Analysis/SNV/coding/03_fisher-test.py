"""
Conduct Fisher's exact test on VEP-counts table (generated in 01_VEP-counts-export.py) for case/control burden

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

def main(args):
    # Initialize
    hl.init(default_reference = 'GRCh38',
            tmp_dir = args.tmp)
        
    # Read in MT
    print(f"Reading in MT ({args.mt})...")
    mt = hl.read_matrix_table(args.mt) 
    mt = mt.repartition(500, shuffle = False)
    print(f"Starting (variants, samples): {mt.count()}")
    counts = mt.aggregate_cols(hl.agg.counter(mt.case_con))
    print(counts)


    if args.individual:
        print(f"Conducting case/control Fishers exact tests for each individual gene...")
        ## Output Fisher's exact tests for each individual gene
        counts = mt.annotate_rows(case_carriers = hl.agg.filter(mt.case_con == "CASE", hl.agg.count_where(mt.agg > 0)),
                                  control_carriers = hl.agg.filter(mt.case_con == "CTRL", hl.agg.count_where(mt.agg > 0)),
                                  case_non_carriers = hl.agg.filter(mt.case_con == "CASE", hl.agg.count_where(mt.agg == 0)),
                                  control_non_carriers = hl.agg.filter(mt.case_con == "CTRL", hl.agg.count_where(mt.agg == 0))) # Annotate in counts
        
        counts = counts.annotate_rows(fisher = (hl.expr.functions.fisher_exact_test(hl.int(counts.case_carriers), hl.int(counts.case_non_carriers),
                                                                                    hl.int(counts.control_carriers), hl.int(counts.control_non_carriers)))) # Test and annotate
        
        print(f"Writing result to: {args.out + args.file_prefix + '_case-control_Fisher-exact_individual.tsv'}") 
        counts.rows().flatten().export(args.out + args.file_prefix + '_case-control_Fisher-exact_individual.tsv') # Export



    ## Next, conduct individual gene-set burden analyses
    if args.gene_lists:
        results = {}
        gene_lists = args.gene_lists.replace(" ", "").split(",") # Parse command-line input
        for gene_list in gene_lists: # Iterate through every desired gene list analysis
            m = mt

            print(f"Analyzing {gene_list}: {gene_lists_info[gene_list]['Description']}...")
            if gene_list != "all": # All just considers all, so no need to filter
                g = hl.import_table(gene_lists_info[gene_list]["List path"], 
                                    delimiter = "\t", key = "gene_symbol", impute = True) # Read in list of gene symbols
                m = m.filter_rows(hl.is_defined(g[m.gene_symbol]), keep = True) # Filter down to genes of interest
            
            print(f"Total number of genes: {m.rows().count()}")
            m = m.annotate_cols(carrier = hl.agg.count_where(m.agg > 0))

            case_carriers = m.aggregate_cols(hl.agg.filter(m.case_con == "CASE", hl.agg.sum(m.carrier)))
            control_carriers = m.aggregate_cols(hl.agg.filter(m.case_con == "CTRL", hl.agg.sum(m.carrier)))
            case_non_carriers = m.aggregate_cols(hl.agg.filter(m.case_con == "CASE", hl.agg.sum(m.carrier == 0)))
            control_non_carriers = m.aggregate_cols(hl.agg.filter(m.case_con == "CTRL", hl.agg.sum(m.carrier == 0)))
            res = hl.eval(hl.fisher_exact_test(case_carriers, case_non_carriers,
                                                control_carriers, control_non_carriers))

            results[gene_list] = [res.p_value, res.odds_ratio, res.ci_95_lower, res.ci_95_upper,
                                  case_carriers, case_non_carriers, control_carriers, control_non_carriers]
            
        df = pd.DataFrame.from_dict(results, orient = "index",
                                    columns = ["p_value", "odds_ratio", "ci_95_lower", "ci_95_upper",
                                               "case_carriers", "case_non_carriers", "control_carriers", "control_non_carriers"])
        print(df)
        print(f"Writing result to: {args.out + args.file_prefix + '_case-control_Fisher-exact_' + '-'.join(gene_lists) + '.tsv'}") 
        df.to_csv(args.out + args.file_prefix + '_case-control_Fisher-exact_' + '-'.join(gene_lists) + '.tsv',
                  sep = "\t", index = False)

            


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mt",
        help = "Path to starting VEP-counts MT (computed in 02_VEP-counts-export.py)",
        type = str,
        required = True
    )
    """
    TODO: accommodate MT's without case/control annotated in already, also maybe ancestry/chip

    parser.add_argument(
        "--manifest",
        help = "Path to manifest keyed by samples ('s') with relevant annotations: 'CASECON', 'population_inference.pop', 'chip'",
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
    """
    parser.add_argument(
        "--individual",
        help = "Whether to conduct case/control Fishers exact tests for each individual gene",
        action = 'store_true'
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