import hail as hl
ht = hl.read_table("gs://2024-wgspd/sv/20241104_WGSPD_Xuefang/gnomAD_SV_v3.ht") # (2063885)
mt = hl.read_matrix_table("gs://2024-wgspd/sv/20241104_WGSPD_sample-GT/gnomAD_SV_v3.release_4_1.1KGP_merged.mt") # (2154486, 63046)
manifest = hl.import_table("gs://2024-wgspd/files/20240905_WGSPD_final-qcd-manifest.tsv", key = "s") # (29124)
ht = ht.key_by(ht.sv_id)


#############################################
# Process aggregate SV HT
#############################################
"""
{'CNV': 647,
 'CPX': 12541,
 'CTX': 91,
 'DEL': 598193,
 'DEL:ME:HERVK': 107,
 'DEL:ME:LINE1': 482,
 'DUP': 243317,
 'INS': 72821,
 'INS:ME:ALU': 167875,
 'INS:ME:LINE1': 28763,
 'INS:ME:SVA': 17084,
 'INV': 2119}

"""
# Filter to passing
ht = ht.filter(ht.filters == ["PASS"], keep = True) # (1144040)

# Filter to deletions
ht = ht.filter(ht.sv_type == "DEL", keep = True) # (598193)

# Impose length filter
ht = ht.annotate(SVLEN = hl.int(ht.info["SVLEN"]))
ht = ht.filter(ht.SVLEN > 5000, keep = True) # (129491)

# Filter to AC 5
ht = ht.annotate(AC = hl.int(ht.info["AC"]))
ht = ht.filter(ht.AC <= 2, keep = True) # (68187), ~0.01% MAF




# Filter MT to releasable samples
# NOTE: not all samples from final-qcd manifest are present in MT, some are not releasable
# CASECON releasable unreleasable
#   CASE       6866         1061
#   CTRL      15969         5228
# Filter to samples in manifest
mt = mt.filter_cols(hl.is_defined(manifest[mt.s]), keep = True) # (2154486, 22835)


# Filter MT to wanted SVs
mt = mt.filter_rows(hl.is_defined(ht[mt.rsid]), keep = True) # (89262, 22835)