"""
Grab hard filtered samples from gnomadv3 VDS
"""
import argparse
import sys
import pandas as pd
import hail as hl

# https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v3/resources
from gnomad_qc.v3.resources.basics import get_gnomad_v3_vds

hl.init(gcs_requester_pays_configuration = 'wes-bipolar',
        tmp_dir = 'gs://wes-bipolar-tmp-4day',
        default_reference = 'GRCh38')


vds = get_gnomad_v3_vds(remove_hard_filtered_samples=True)
hf_s = vds.variant_data.s.collect()

vds = get_gnomad_v3_vds(remove_hard_filtered_samples=False)
all_s = vds.variant_data.s.collect()

hard_filtered = set(all_s).difference(set(hf_s))

with open("../files/hard-filtered-samples", "w") as f:
    for a in hard_filtered:
        f.write(f"{a}\n")