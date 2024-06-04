"""
Get VEP annotations for all variants in a provided VDS and write to specified HT for later

Make sure to start hail cluster with VEP (and probably gnomad). Something like:
hailctl dataproc start [cluster-name] --vep GRCh38 --packages gnomad
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
    # Initialize=
    hl.init(default_reference = 'GRCh38',
            tmp_dir = args.tmp)
        
    # Read in VDS 
    logger.info(f"Reading in VDS ({args.vds})...")
    ht = hl.vds.read_vds(args.vds).variant_data.rows() # Only interested in variants from the variant_data
    logger.info(f"n variants: {ht.count()}")

    # VEP annotation
    logger.info(f"Generating VEP annotations...")
    ht_vep = hl.vep(ht, args.vep)

    # Write VEP annotations to HT output
    logger.info(f"Writing VEP annotations to {args.out}...")
    ht_vep.write(args.out, overwrite = True)

    logger.info("Copying log to logging bucket...")
    hl.copy_log(f"{args.tmp}logs/{date.today()}_00_VEP-HT-export.log")




if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--vds",
        help = "Path to input VDS for which we want to generate VEP annotations",
        type = str,
        required = True
    )
    parser.add_argument(
        "--vep",
        help = "Path to the VEP annotation .json",
        type = str,
        required = True
    )
    parser.add_argument(
        "--out",
        help = "Path to the VEP HT",
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