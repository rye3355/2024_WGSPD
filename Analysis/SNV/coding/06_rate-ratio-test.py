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

def compute_rate_ratio(m: hl.MatrixTable, counts: hl.DictExpression):
    m = m.annotate_cols(total_obs = hl.agg.sum(m.agg)) # Label samples as carriers (any count > 0)
    res = m.cols()
    a = res.aggregate(hl.agg.filter(res.case_con == "CASE", hl.agg.sum(res.total_obs))) # Count carriers
    b = res.aggregate(hl.agg.filter(res.case_con == "CTRL", hl.agg.sum(res.total_obs))) # Count non-carriers
    A = hl.eval(counts.get("CASE"))
    B = hl.eval(counts.get("CTRL"))
    # {fmsb} rateratio
    rate_ratio = (a/A) / (b/B)
    norm_p = hl.qnorm(1 - (1-0.95)/2)
    c = (a - (A/(A+B)) * (a+b))/hl.sqrt((a+b) * (A/(A+B)) * (B/(A+B)))
    p_value = 2 * (1 - hl.pnorm(hl.abs(c)))
    se = hl.sqrt((1/a) + (1/b))
    ci_95_lower = rate_ratio * hl.exp(-norm_p*se)
    ci_95_upper = rate_ratio * hl.exp(norm_p*se)
    p_value = 2*(1-hl.pnorm(hl.abs(hl.log(rate_ratio) / se)))
    return [rate_ratio, ci_95_lower, ci_95_upper, p_value, se, 
            a, b, A, B]




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

        results = {}
        for gene_list in gene_lists: # Iterate through every desired gene list analysis
            print(f"\nAnalyzing {gene_list}: {gene_lists_info[gene_list]['Description']}\n")
            m = mt
            if gene_list != "all":
                g = hl.import_table(gene_lists_info[gene_list]["List path"], 
                                    delimiter = "\t", key = "gene_symbol", impute = True) # Read in list of gene symbols
                m = m.filter_rows(hl.is_defined(g[m.gene_symbol]), keep = True) # Filter down to genes of interest
                
            print(f"\nTotal number of genes: {m.rows().count()}\n")    
            results[gene_list] = [gene_list] + compute_rate_ratio(m, counts)

        # Only write if something to be written at this point
        if results:
            df = pd.DataFrame.from_dict(results, orient = "index",
                                        columns = ["gene_set", "rate_ratio", "ci_95_lower", "ci_95_upper", "p_value", "se", 
                                                   "case_counts", "ctrl_counts", "n_case", "n_ctrl"]
                                        )
            df_ht = hl.Table.from_pandas(df, key = "gene_set")
            print(f"\nWriting result to: {args.out + args.file_prefix + '_rate-ratio_' + '-'.join(gene_lists) + '.tsv'}\n") 
            df_ht.export(args.out + args.file_prefix + '_rate-ratio_' + '-'.join(gene_lists) + '.tsv', delimiter='\t')

            


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