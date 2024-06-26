"""
Get counts for specific VEP annotations and write 


gs://2024-wgspd/snv/coding/202240618_subset_post-qc_protein-coding.mt

"""

import hail as hl
import logging
import argparse
from datetime import date

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s: %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def main(args):
    # Initialize
    hl.init(default_reference = 'GRCh38',
            tmp_dir = args.tmp)
        
    # Read in MT
    logger.info(f"Reading in MT ({args.mt})...")
    mt = hl.read_matrix_table(args.mt) 
    logger.info(f"Starting (variants, samples): {mt.count()}")

    # Variant QC (for cohort-specific AC filtering)
    mt = hl.variant_qc(mt)

    # Read in VEP HT
    ht = hl.read_table(args.vep_ht)

    # Filter to variants in VEP HT (usually already filtered to consequence category in pLoF, other_missense, damaging_missense, synonymous)
    mt = mt.filter_rows(hl.is_defined(ht[mt.locus, mt.alleles]), keep = True)

    # Annotate in relevant fields from VEP HT
    mt = mt.annotate_rows(gene_symbol = ht[mt.locus, mt.alleles].vep.worst_csq_for_variant_canonical.gene_symbol,
                          consequence_category = ht[mt.locus, mt.alleles].consequence_category,
                          type = ht[mt.locus, mt.alleles].type,
                          infrIndel = ht[mt.locus, mt.alleles].infrIndel)



    # Filter rows based on specified input flags
    f = [args.file_prefix] # Running identifier/run details


    # First, filter to variants below gnomAD and/or RGC AC thresholds
    if args.gnomAD_AC_thresh:
        logger.info(f"Filtering to gnomAD AC <= {args.gnomAD_AC_thresh}...")
        GNOMAD_AC = "gs://gcp-public-data--gnomad/release/4.1/ht/genomes/gnomad.genomes.v4.1.sites.ht/"
        gnomad_AC = hl.read_table(GNOMAD_AC)
        gnomad_AC = gnomad_AC.annotate(raw_freq_AC = gnomad_AC.freq[gnomad_AC.freq_index_dict['raw']].AC)
        mt = mt.annotate_rows(gnomad_ac = gnomad_AC[mt.locus, mt.alleles].raw_freq_AC)
        mt = mt.filter_rows(mt.gnomad_ac <= args.gnomAD_AC_thresh, keep = True)
        f.append(f"gnomadAC{args.gnomAD_AC_thresh}")
    if args.RGC_AC_thresh:
        logger.info(f"Filtering to RGC AC <= {args.RGC_AC_thresh}...")
        RGC_AC = "gs://bipex2/annotations/RGC/rgc_me_variant_frequencies.ht"
        rgc_AC = hl.read_table(RGC_AC)
        mt = mt.annotate_rows(rgc_ac = rgc_AC[mt.locus, mt.alleles].ALL_AC)
        mt = mt.filter_rows(mt.rgc_ac <= args.RGC_AC_thresh, keep = True)



    ## Next, consequence_category
    if args.cons_cat:
        cons_cat = args.cons_cat.replace(" ", "").split(",") # Parse command-line input
        logger.info(f"Filtering to {cons_cat}...")
        mt = mt.filter_rows(hl.set(cons_cat).contains(mt.consequence_category), keep = True)
        f.append(f"{'-'.join(cons_cat)}")

    ## Next, AC threshold
    if args.ac == 1:
        logger.info(f"Filtering to singletons...")
        mt = mt.filter_rows(mt.variant_qc.n_non_ref==1, keep = True)
        f.append(f"singletons")
    elif args.ac:
        logger.info(f"Filtering to AC <= {args.ac}...")
        mt = mt.filter_rows(mt.variant_qc.AC[1] <= args.ac, keep = True)
        f.append(f"AC{args.ac}")
    

    ## Next, missense-specific thresholds
    if args.mpc:
        logger.info(f"Filtering to MPC >= {args.mpc}...")
        mt = mt.filter_rows(hl.set(["damaging_missense", "other_missense"]).contains(mt.consequence_category) & # Missense 
                            mt.MPC < args.mpc, # And MPC under threshold
                            keep = False) # Filter out
        f.append(f"MPC{args.mpc}")
    if args.am:
        logger.info(f"Filtering to AM >= {args.am}...")
        mt = mt.filter_rows(hl.set(["damaging_missense", "other_missense"]).contains(mt.consequence_category) & # Missense 
                            mt.AM < args.am, # And AM under threshold
                            keep = False) # Filter out
        f.append(f"AM{args.am}")
    
    
    # Finally, aggregate/count and write
    if args.non_gnomAD_psych:
        logger.info(f"Filtering to non-GnomAD...")
        mt_agg = mt.group_rows_by(mt.gene_symbol).aggregate(
            agg = hl.agg.count_where(mt.GT.is_non_ref() & ~mt.inGnomAD_nonpsych)
        )

        logger.info(f"Writing to {args.out + '_'.join(f) + '_non-gnomAD-psych_counts.mt'}...")
        mt_agg.write(args.out + '_'.join(f) + '_non-gnomAD-psych_counts.mt', overwrite = True)


    logger.info(f"Generating counts...")
    mt_agg = mt.group_rows_by(mt.gene_symbol).aggregate(
        agg = hl.agg.count_where(mt.GT.is_non_ref())
    )

    logger.info(f"Writing to {args.out + '_'.join(f) + '_counts.mt'}...")
    mt_agg.write(args.out + '_'.join(f) + '_counts.mt', overwrite = True)

    logger.info(f"Copying log to {args.out + 'logs/' + str(date.today()) + '_' + '_'.join(f) + '.log'}...")
    hl.copy_log(f"{args.out + 'logs/' + str(date.today()) + '_' + '_'.join(f) + '.log'}")




if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mt",
        help = "Path to input MT (computed in 00_filter-to-coding.py)",
        type = str,
        required = True
    )
    parser.add_argument(
        "--vep_ht",
        help = "Path to VEP-annotation HT for input MT (computed in 01_vep-annotate_[run].py)",
        type = str,
        required = True
    )
    parser.add_argument(
        "--gnomAD_AC_thresh",
        help = "gnomAD AC threshold to filter to (keep entries <= thresh)",
        type = int,
        required = False
    )   
    parser.add_argument(
        "--RGC_AC_thresh",
        help = "RGC AC threshold to filter to (keep entries <= thresh)",
        type = int,
        required = False
    )  
    parser.add_argument(
        "--cons_cat",
        help = "Comma-separated consequence categories to filter to: pLoF, other_missense, damaging_missense, synonymous",
        type = str,
        required = False
    )
    parser.add_argument(
        "--ac",
        help = "AC threshold to filter to (remove entries < thresh), usually 5 or 10",
        type = int,
        required = False
    )
    parser.add_argument(
        "--mpc",
        help = "MPC threshold to filter to (remove entries < thresh), usually 2 (with AM 0.98)",
        type = int,
        required = False
    )
    parser.add_argument(
        "--am",
        help = "AlphaMissense threshold to filter to (>= thresh), usually 0.98 (with MPC 2)",
        type = float,
        required = False
    )    
    parser.add_argument(
        "--non_gnomAD_psych",
        help = "Flag for filtering to non-gnomAD-psych variants",
        action = 'store_true'
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