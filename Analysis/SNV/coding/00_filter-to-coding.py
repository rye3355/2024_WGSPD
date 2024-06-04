"""
Filter to protein coding genes as defined by Gencode v46
Also min rep

https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz
"""
import hail as hl

hl.init(default_reference = 'GRCh38',
                tmp_dir = "gs://wes-bipolar-tmp-4day/")

"""

# Read original data MT
MT = 'gs://2024-wgspd/qc/20240426_subset_final-qcd.mt'
mt = hl.read_matrix_table(MT)
print(f"Original data count: {mt.count()}")
#Original data count: (332778683, 28554)


# Annotate case/control info in
MANIFEST = 'gs://2024-wgspd/20240523_WGSPD_final-qcd-manifest.tsv'
manifest = hl.import_table(MANIFEST, delimiter='\t',
                          key = "s", impute = True)
mt = mt.annotate_cols(is_case = manifest[mt.s].CASECON == "CASE")



# Protein coding
## Read Gencode annotations (https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz)
GTF = 'gs://2024-wgspd/files/gencode.v46.annotation.gtf'
gtf = hl.experimental.import_gtf(GTF,
                                reference_genome='GRCh38',
                                skip_invalid_contigs=True)
print(f"Number of genes from Gencode: {gtf.count()}")
#Number of genes from Gencode: 3467156

## Filter Gencode table down to protein coding (gene/transcript type)
protein_coding = gtf.filter(gtf.gene_type == "protein_coding")
print(f"Number of protein-coding (gene_type) genes from Gencode: {protein_coding.count()}")
#Number of protein-coding (gene_type) genes from Gencode: 3078241
protein_coding = protein_coding.filter(protein_coding.transcript_type == "protein_coding")
print(f"Number of protein-coding (gene/transcript type) genes from Gencode: {protein_coding.count()}")
#Number of protein-coding (gene/transcript type) genes from Gencode: 2113938

# Filter variants in MT 
mt = mt.filter_rows(hl.is_defined(protein_coding[mt.locus]))
print(f"Final data (coding) count: {mt.count()}")
#Final data (coding) count: (144360038, 28554)


# Min rep
mt = mt.annotate_rows(min_rep = hl.min_rep(mt.locus, mt.alleles))
mt = mt.key_rows_by('min_rep')
mt = mt.drop('locus', 'alleles')
mt = mt.annotate_rows(locus = mt.min_rep.locus,
                      alleles = mt.min_rep.alleles)
mt = mt.key_rows_by('locus', 'alleles')
mt = mt.drop('min_rep')
print(f"min rep count: {mt.count()}")
#min rep count: (144360038, 28554)


# Write
mt.write("gs://2024-wgspd/snv/coding/202240604_subset_post-qc_protein-coding.mt", overwrite = True)


"""

mt = hl.read_matrix_table("gs://2024-wgspd/snv/coding/202240604_subset_post-qc_protein-coding.mt")
manifest = hl.import_table('gs://2024-wgspd/files/20240523_WGSPD_final-qcd-manifest.tsv', delimiter='\t',
                          key = "s", impute = True)
mt = mt.annotate_cols(is_case = manifest[mt.s].CASECON == "CASE")
mt.write("gs://2024-wgspd/snv/coding/202240603_subset_post-qc_protein-coding.mt", overwrite = True)