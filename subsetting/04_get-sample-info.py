"""
Use gnomad sample HT to get info for our subset
"""
import argparse
import sys
import pandas as pd
import hail as hl

from gnomad_qc.v3.resources.meta import meta as metadata

hl.init(gcs_requester_pays_configuration = 'wes-bipolar',
        tmp_dir = 'gs://wes-bipolar-tmp-4day',
        default_reference = 'GRCh38')

sample_ht = hl.import_table("gs://gnomad-subsets-2024/gnomad-v3/202402/scz_bp_samplelist_gnomad3_FINAL.tsv", key="s")

meta = metadata.ht()
meta = meta.filter(hl.is_defined(sample_ht[meta.s]))

meta.export("gs://2024-wgspd/gnomad_v3.1_subset-metadata.tsv", delimiter = "\t")