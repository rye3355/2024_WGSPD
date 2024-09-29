

import hail as hl
hl.init(tmp_dir='gs://wes-bipolar-tmp-4day/20240928/generate-counts/',
        default_reference = 'GRCh38')

# Input
MT= "gs://2024-wgspd/snv/coding/202240928_subset_post-qc_protein-coding.mt"
SITES_TABLE = 'gs://2024-wgspd/snv/coding/20240928_subset_post-qc_protein-coding_VEP-annotated.ht'


############
print(f"Reading {MT}")
mt = hl.read_matrix_table(MT)
print(f"Initial: {mt.count()}")
# Initial: (143350801, 29124)


ht_site = hl.read_table(SITES_TABLE)
mt = mt.annotate_rows(site = ht_site[mt.row_key])


mt = hl.variant_qc(mt, name = "generate_counts_variant_qc")
mt = mt.filter_rows(hl.min(mt["generate_counts_variant_qc"].AC) != 0, keep = True)
#print(f"After filtering invariant: {mt.count()}")
# After filtering invariant: (133573182, 29124)

# gnomad AC
gnomad_AC = hl.read_table("gs://gcp-public-data--gnomad/release/4.1/ht/joint/gnomad.joint.v4.1.sites.ht/")
gnomad_AC = gnomad_AC.annotate(raw_freq_AC = gnomad_AC.joint.freq[gnomad_AC.joint_globals.freq_index_dict['raw']].AC)
mt = mt.annotate_rows(raw_freq_AC = hl.if_else(hl.is_defined(gnomad_AC[mt.locus, mt.alleles]), # Check if locus/allele is even in gnomad table
                                             gnomad_AC[mt.locus, mt.alleles].raw_freq_AC,  # If so, use the value
                                             0), # Otherwise, set to AC = 0
                        )

mt = mt.annotate_rows(
    # Consequence category
    ptv = mt.site.consequence_category == "pLoF",
    mis1 = mt.site.consequence_category == "other_missense",
    mis3 = mt.site.consequence_category == "damaging_missense",
    syn = mt.site.consequence_category == "synonymous",
    noncoding = mt.site.consequence_category == "non_coding",
    # gnomad AC
    gnomad_AC = mt.raw_freq_AC <= 10,
    # AC
    isSingleton = (mt["generate_counts_variant_qc"].n_non_ref==1),
    AC= hl.min(mt["generate_counts_variant_qc"].AC) <= 10,
)

mt = mt.filter_rows(mt.AC & mt.ptv & mt.gnomad_AC, keep = True)

mt2 = mt.group_rows_by(mt.site.vep.worst_csq_for_variant_canonical.gene_symbol).aggregate(
    agg = hl.agg.count_where(mt.GT.is_non_ref()),
    )

print("Export counts table")
mt2.repartition(2000, shuffle = False).write("gs://2024-wgspd/snv/coding/outputs/202409028_WCPG/20240929_WGSPD_pLoF-MAC10-gnomadAC10", overwrite = True)
