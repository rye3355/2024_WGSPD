"""
Just filter to subset and write
"""
import argparse
import sys
import pandas as pd
import hail as hl

# https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v3/resources
from gnomad_qc.v3.resources.annotations import get_info
from gnomad_qc.v3.resources.basics import get_gnomad_v3_vds

# https://github.com/broadinstitute/gnomad_methods/tree/main/gnomad
from gnomad.resources.grch38.gnomad import GENOME_POPS
from gnomad.utils.annotations import annotate_adj

hl.init(gcs_requester_pays_configuration = 'wes-bipolar',
        tmp_dir = 'gs://wes-bipolar-tmp-4day',
        default_reference = 'GRCh38')


vds = get_gnomad_v3_vds(remove_hard_filtered_samples=True)

sample_ht = hl.import_table("gs://gnomad-subsets-2024/gnomad-v3/202402/scz_bp_samplelist_gnomad3_FINAL.tsv", key="s")
vds = hl.vds.filter_samples(vds, sample_ht, keep = True, remove_dead_alleles = True)


vds.write("gs://gnomad-subsets-2024/gnomad-v3/v3_filtered-samples.vds")