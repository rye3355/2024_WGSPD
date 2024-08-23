"""
Conduct Firth regression test on VEP-counts table (generated in 01_VEP-counts-export.py) 
Currently only supports stratification by ancestry x chip
"""

import hail as hl
import pandas as pd
import argparse
from datetime import date


def main(args):
    # Initialize
    hl.init(default_reference = 'GRCh38',
            tmp_dir = args.tmp)
        
    # Read in MT
    print(f"\nReading in MT ({args.mt})...\n")
    mt = hl.read_matrix_table(args.mt) 
    print(f"\nStarting (genes, samples): {mt.count()}\n")


    # Annotate and process with metadata from manifest
    # Case/control, Population, Chip
    if (args.annotate_casecon | args.annotate_pop | args.annotate_chip | args.filter_pass) and not args.manifest: # Check annotation requirements
        print(f"No manifest provided to pull annotations from...")
        return
    
    # Annotate in fields from manifest
    if args.manifest:
        manifest = hl.import_table(args.manifest,
                                    delimiter = "\t", impute = True, key = "s") # Requires key by "s" and specific fields "CASECON", "POP", and "CHIP" as needed
        if args.filter_pass:
            mt = mt.filter_cols(hl.if_else(hl.is_defined(manifest[mt.s]), # If it's even in manifest
                                                manifest[mt.s][args.filter_pass] == args.pass_value, # Only True if PASS
                                                False), # Otherwise, don't keep
                                                keep = True)
        if args.annotate_casecon: 
            mt = mt.annotate_cols(case_con = manifest[mt.s][args.annotate_casecon]) # Annotate all with "CASE" or "CTRL"
            mt = mt.annotate_cols(is_CASE = mt.case_con == "CASE")
        if args.annotate_pop:
            mt = mt.annotate_cols(pop = manifest[mt.s][args.annotate_pop]) # Annotate all with population
            mt = mt.annotate_cols(is_AFR = hl.if_else(mt.pop == "AFR", 1, 0),
                                  is_AMI = hl.if_else(mt.pop == "AMI", 1, 0),
                                  is_AMR = hl.if_else(mt.pop == "AMR", 1, 0),
                                  is_ASJ = hl.if_else(mt.pop == "ASJ", 1, 0),
                                  is_CSA = hl.if_else(mt.pop == "CSA", 1, 0),
                                  is_EAS = hl.if_else(mt.pop == "EAS", 1, 0),
                                  is_FIN = hl.if_else(mt.pop == "FIN", 1, 0),
                                  is_MID = hl.if_else(mt.pop == "MID", 1, 0),
                                  is_NFE = hl.if_else(mt.pop == "NFE", 1, 0),
                                  is_OTH = hl.if_else(mt.pop == "OTH", 1, 0),
                                  is_SAS = hl.if_else(mt.pop == "SAS", 1, 0)
                                  )
        if args.annotate_chip:
            mt = mt.annotate_cols(chip = manifest[mt.s][args.annotate_chip]) # Annotate all with chip
            mt = mt.annotate_cols(is_Agilent = hl.if_else(mt.chip == "Agilent", 1, 0),
                                  is_Nextera = hl.if_else(mt.chip == "Nextera", 1, 0),
                                  is_Twist = hl.if_else(mt.chip == "Twist", 1, 0),
                                  is_WGS = hl.if_else(mt.chip == "WGS", 1, 0),
                                  is_broad_custom_exome_v1 = hl.if_else(mt.chip == "broad_custom_exome_v1", 1, 0),
                                  )
        
    if args.syn_count:
        syn_counts = hl.import_table(args.syn_count,
                                    delimiter = "\t", impute = True, key = "s")
        mt = mt.annotate_cols(syn_count = syn_counts[mt.s]["syn_count"])
        
    if args.population_only_strat:
        mt = mt.annotate_cols(group = mt.pop,
                              group2 = mt.pop + "_" + mt.case_con)
    else:       
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





    pcs = hl.import_table("gs://bipex2/202407_ancestry-relatedness/outputs/01_bipex2_pca-with-ref_scores_0.90.tsv", 
                          delimiter = "\t", impute = True, key = "s")
    
    mt = mt.annotate_cols(PC1 = pcs[mt.s].PC1,
                          PC2 = pcs[mt.s].PC2,
                          PC3 = pcs[mt.s].PC3,
                          PC4 = pcs[mt.s].PC4,
                          PC5 = pcs[mt.s].PC5,
                          PC6 = pcs[mt.s].PC6,
                          PC7 = pcs[mt.s].PC7,
                          PC8 = pcs[mt.s].PC8,
                          PC9 = pcs[mt.s].PC9,
                          PC10 = pcs[mt.s].PC10)
    ## Next, conduct individual gene-set burden analyses
    if args.gene_lists:
        gene_lists = args.gene_lists.replace(" ", "").split(",") # Parse command-line input


        # Running Firth on each gene individually
        # TODO: Make this more efficient with row aggregations?
        if "individual" in gene_lists:
            genes = mt.gene_symbol.collect() 
            print(f"\nAnalyzing {len(genes)} genes individually\n")
            m = mt
            #m = mt.head(n_rows = 100, n_cols = 10000)
            result_ht = hl.logistic_regression_rows(test = "firth", 
                                                    y = m.is_CASE,
                                                    x = m.agg > 0,
                                                    covariates = [1, 
                                                                  m.is_AFR, m.is_AMI, m.is_AMR, m.is_ASJ, m.is_CSA, m.is_EAS,
                                                                  m.is_FIN, m.is_MID, m.is_NFE, m.is_OTH, m.is_SAS,
                                                                  m.is_Agilent, m.is_Nextera, m.is_Twist, m.is_WGS, m.is_broad_custom_exome_v1,
                                                                  m.syn_count])
            result_ht = hl.logistic_regression_rows(test = "firth", 
                                                    y = m.is_CASE,
                                                    x = m.agg > 0,
                                                    covariates = [1, 
                                                                  m.PC1, m.PC2, m.PC3, m.PC4, m.PC5, 
                                                                  m.PC6, m.PC7, m.PC8, m.PC9, m.PC10,
                                                                  m.syn_count])


            result_ht = result_ht.checkpoint("gs://bipex2/202407_SNV/20240806/20240806_BipEx2_PTV-miss_firth_output.ht", overwrite = True) # PTV+Missense, PCs + syn count
            result_ht = result_ht.flatten().drop("firth_null_fit.b", "firth_null_fit.mu", "fit.b", "fit.mu")
            result_ht.export("gs://bipex2/202407_SNV/20240805_test-firth/20240806_BipEx2_PTV-miss_firth_pcs_output.tsv")

            
            result_ht = result_ht.checkpoint("gs://wes-bipolar-tmp-4day/20240805/firth_output.ht", overwrite = True) # Syn AC 5, PCs + syn count
            result_ht = result_ht.checkpoint("gs://wes-bipolar-tmp-4day/20240805/firth_output_pcs.ht", overwrite = True) # Syn AC 5, PCs
            result_ht = result_ht.flatten().drop("firth_null_fit.b", "firth_null_fit.mu", "fit.b", "fit.mu")
            result_ht.export("gs://bipex2/202407_SNV/20240805_test-firth/firth_output.tsv") #  Syn AC 5, PCs + syn count
            result_ht.export("gs://bipex2/202407_SNV/20240805_test-firth/synAC5_firth_output_pcs.tsv") #  Syn AC 5, PCs
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
        help = "Flag to filter to samples with 'PASS' in manifest; flag value is name of field in manifest",
        type = str,
        required = False
    )
    parser.add_argument(
        "--pass_value",
        help = "Value to keep in filter_pass field",
        type = str,
        required = False
    )
    parser.add_argument(
        "--annotate_casecon",
        help = "Flag to annotate in 'CASE' or 'CTRL' for each sample; flag value is name of field in manifest",
        type = str,
        required = False
    )
    parser.add_argument(
        "--annotate_pop",
        help = "Flag to annotate in population for each sample; flag value is name of field in manifest",
        type = str,
        required = False
    )
    parser.add_argument(
        "--annotate_chip",
        help = "Flag to annotate in chip for each sample; flag value is name of field in manifest",
        type = str,
        required = False
    )
    parser.add_argument(
        "--syn_count",
        help = "Path to tab-delimited manifest keyed by samples ('s') with relevant synonymous counts: 'syn_count'",
        type = str,
        required = False
    )
    parser.add_argument(
        "--population_only_strat",
        help = "Flag to stratify by population only (not chip)",
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