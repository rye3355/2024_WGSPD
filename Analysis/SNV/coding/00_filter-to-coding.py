"""
Filter to protein coding genes as defined by Gencode v46
Also min rep

Also annotate samples with case/control assignments,
variants with MPC, AM, in gnomAD nonpsych, and in discovEHR

https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz
"""
import hail as hl

hl.init(default_reference = 'GRCh38',
                tmp_dir = "gs://wes-bipolar-tmp-4day/")

# Read original data MT
MT = 'gs://2024-wgspd/qc/20240426_subset_final-qcd.mt'
mt = hl.read_matrix_table(MT)
print(f"Original data count: {mt.count()}")
#Original data count: (332778683, 28554)


# Annotate case/control info in
MANIFEST = 'gs://2024-wgspd/files/20240523_WGSPD_final-qcd-manifest.tsv'
manifest = hl.import_table(MANIFEST, delimiter='\t',
                          key = "s", impute = True)
mt = mt.annotate_cols(case_con = manifest[mt.s].CASECON)



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


# Updated MPC annotation
MPC = 'gs://bipex2/annotations/mpc_grch38_deduped_with_outliers_2024-04-30.ht'
mpc = hl.read_table(MPC)
mpc = mpc.key_by('locus', 'alleles')
mt = mt.annotate_rows(MPC = mpc[mt.locus, mt.alleles].mpc)

# gnomAD nonpsych
GNOMAD_NONPSYCH='gs://raw_data_bipolar_dalio_w1_w2/inputs/gnomad.exomes.r2.1.1.non_psych_sites_GRCh38.ht'
gnomAD_nonpsych = hl.read_table(GNOMAD_NONPSYCH)
mt = mt.annotate_rows(inGnomAD_nonpsych = hl.is_defined(gnomAD_nonpsych[mt.locus, mt.alleles]))

# discovEHR
DISCOVEHR = 'gs://bd_scz/BGE_Callset_PAISA_QIMR_NeuroMex_KenyaPsych/annotations/DiscovEHR_GHS_Freeze_50.L3DP10.pVCF.frq_sites_grch38.ht'
discovehr = hl.read_table(DISCOVEHR)
mt = mt.annotate_rows(inDiscovEHR = hl.is_defined(discovehr[mt.locus, mt.alleles].info))

# OS
OS = 'gs://bipex2/annotations/20240618_GRCh38_OS.ht'
os = hl.read_table(OS)
mt = mt.annotate_rows(inOS = hl.is_defined(os[mt.locus, mt.alleles]))

# Other miss annotations
OTH = 'gs://bipex2/annotations/AlphaMissense_deduplicated_hg38_with_ps_mf_2024-05-21.ht'
oth = hl.read_table(OTH)
oth = oth.key_by(oth.locus, oth.alleles)
mt = mt.annotate_rows(AM = oth[mt.locus, mt.alleles].am_pathogenicity,
                      score_ml = oth[mt.locus, mt.alleles].score_ml,
                      MisFit_D = oth[mt.locus, mt.alleles].MisFit_D,
                      MisFit_S = oth[mt.locus, mt.alleles].MisFit_S,
                      protein_variant = oth[mt.locus, mt.alleles].protein_variant)


# Write
mt.write("gs://2024-wgspd/snv/coding/202240618_subset_post-qc_protein-coding.mt", overwrite = True)

"""
print(mt.aggregate_rows(hl.agg.counter(mt.inOS)))

#{'NAGNAG_SITE': 754, 'NON_CAN_SPLICE': 339, 'NON_CAN_SPLICE,NAGNAG_SITE': 4, 'PHYLOCSF_UNLIKELY_ORF': 125, 'PHYLOCSF_WEAK': 23076, 'SINGLE_EXON': 3523, 'SINGLE_EXON,PHYLOCSF_WEAK': 928, None: 144331289}
"""

