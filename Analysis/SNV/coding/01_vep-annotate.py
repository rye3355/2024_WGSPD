"""
Annotate with VEP consequences
"""
import hail as hl
from gnomad.utils.vep import vep_or_lookup_vep
from gnomad.utils.vep import process_consequences

hl.init(default_reference = 'GRCh38',
                tmp_dir = "gs://wes-bipolar-tmp-4day/")


# Read MT
MT = 'gs://2024-wgspd/snv/coding/202240604_subset_post-qc_protein-coding.mt'
mt = hl.read_matrix_table(MT)
print(f"Original data count: {mt.count()}")
#Original data count: (144360038, 28554)


# Annotate 
ht = mt.rows()
ht_vep = vep_or_lookup_vep(ht, reference="GRCh38")


ht_vep.write("gs://2024-wgspd/snv/coding/202240604_subset_post-qc_protein-coding_VEP-annotated.ht", overwrite = True)



