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

def run_fisher_on_gene_lists(m: hl.MatrixTable, counts: hl.DictExpression):
    m = m.annotate_cols(carrier = hl.agg.any(m.agg > 0)) # Label samples as carriers (any count > 0)
    res = m.cols()
    case_carriers = res.aggregate(hl.agg.filter(res.case_con == "CASE", hl.agg.sum(res.carrier))) # Count carriers
    case_non_carriers = hl.eval(counts.get("CASE")) - case_carriers # Use dictionary to get non-carrier counts (faster than re-counting)
    control_carriers = res.aggregate(hl.agg.filter(res.case_con == "CTRL", hl.agg.sum(res.carrier))) 
    control_non_carriers = hl.eval(counts.get("CTRL")) - control_carriers 

    res = hl.eval(hl.fisher_exact_test(case_carriers, case_non_carriers,
                                        control_carriers, control_non_carriers))

    return [res.p_value, res.odds_ratio, res.ci_95_lower, res.ci_95_upper,
            case_carriers, case_non_carriers, control_carriers, control_non_carriers]




def main(args):
    # Initialize
    hl.init(default_reference = 'GRCh38',
            tmp_dir = args.tmp)
        
    # Read in MT
    print(f"\nReading in MT ({args.mt})...\n")
    mt = hl.read_matrix_table(args.mt) 

     # Check annotation requirements
    if args.annotate_casecon and not args.manifest:
        print(f"No manifest provided to pull annotations from...")
        return
    
    # Annotate in fields from manifest
    # Requires key by "s" and specific fields "CASECON" as needed
    if args.manifest:
        manifest = hl.import_table(args.manifest,
                                   delimiter = "\t", impute = True, key = "s")
        
        if args.annotate_casecon:
            mt = mt.annotate_cols(case_con = manifest[mt.s].CASECON)

    counts = mt.aggregate_cols(hl.agg.counter(mt.case_con))
    print(counts)
    counts = hl.literal(counts) # Convert counts dictionary to hail format



    if mt.n_partitions() > 200:
        mt = mt.repartition(200, shuffle = False).persist()
    print(f"\nStarting (genes, samples): {mt.count()}\n")
    
    ## Next, conduct individual gene-set burden analyses
    if args.gene_lists:
        gene_lists = args.gene_lists.replace(" ", "").split(",") # Parse command-line input


        # Running CMH on each gene individually
        # TODO: Make this more efficient with row aggregations?
        if "individual" in gene_lists:
            genes = mt.gene_symbol.collect() 
            print(f"\nAnalyzing {len(genes)} genes individually\n")
            m = mt
            res = m.group_cols_by(m.case_con).aggregate(n_carriers = hl.agg.sum(m.agg > 0)).repartition(200, shuffle = False).persist() # Count carriers by CASE/CTRL stratifications
            res = res.annotate_entries(n_non_carriers = counts.get(res.case_con) - res.n_carriers) # Use dictionary to get non-carrier counts (faster than re-counting)
            res_counts = res.annotate_rows(case_carriers = hl.agg.filter(res.case_con == "CASE", hl.agg.collect(res.n_carriers)[0]), 
                                           case_non_carriers = hl.agg.filter(res.case_con == "CASE", hl.agg.collect(res.n_non_carriers)[0]), 
                                           control_carriers = hl.agg.filter(res.case_con == "CTRL", hl.agg.collect(res.n_carriers)[0]), 
                                           control_non_carriers = hl.agg.filter(res.case_con == "CTRL", hl.agg.collect(res.n_non_carriers)[0])).rows() # Collect counts into arrays (collect all ancestry x kit groups), slim to only rows now
            res_counts = res_counts.annotate(fisher = hl.expr.functions.fisher_exact_test(hl.int(res_counts.case_carriers), hl.int(res_counts.case_non_carriers),
                                                                                          hl.int(res_counts.control_carriers), hl.int(res_counts.control_non_carriers))) # Annotate in Fisher results
            res_counts = res_counts.flatten()

            #print(f"\nWriting result to: {args.out + args.file_prefix + '_case-control_Fisher-exact_individual.ht'}\n") 
            #res_counts = res_counts.checkpoint(args.out + args.file_prefix + '_case-control_Fisher-exact_individual.ht')

            print(f"\nWriting result to: {args.out + args.file_prefix + '_case-control_Fisher-exact_individual.tsv'}\n") 
            res_counts.export(args.out + args.file_prefix + '_case-control_Fisher-exact_individual.tsv', delimiter = "\t")
            
            
            gene_lists.remove("individual")


        results = {}
        for gene_list in gene_lists: # Iterate through every desired gene list analysis
            if gene_list == "individual":
                continue
            print(f"\nAnalyzing {gene_list}: {gene_lists_info[gene_list]['Description']}\n")
            m = mt
            g = hl.import_table(gene_lists_info[gene_list]["List path"], 
                                delimiter = "\t", key = "gene_symbol", impute = True) # Read in list of gene symbols
            m = m.filter_rows(hl.is_defined(g[m.gene_symbol]), keep = True) # Filter down to genes of interest
            
            print(f"\nTotal number of genes: {m.rows().count()}\n")    
            results[gene_list] = [gene_list] + run_fisher_on_gene_lists(m, counts)

        # Only write if something to be written at this point
        if results:
            df = pd.DataFrame.from_dict(results, orient = "index",
                                        columns = ["gene_set", "p_value", "odds_ratio", "ci_95_lower", "ci_95_upper",
                                                   "case_carriers", "case_non_carriers", "control_carriers", "control_non_carriers"])
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
        help = "Path to tab-delimited manifest keyed by samples ('s') with relevant annotations: 'CASECON'",
        type = str,
        required = False
    )
    parser.add_argument(
        "--annotate_casecon",
        help = "Flag to annotate in case control status for each sample",
        action = 'store_true'
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