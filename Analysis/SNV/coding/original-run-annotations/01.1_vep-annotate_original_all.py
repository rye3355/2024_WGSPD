"""
Annotate with VEP consequences

Input (from 00_filter-to-coding.py):
gs://2024-wgspd/snv/coding/202240618_subset_post-qc_protein-coding.mt

Output (feed into 02_VEP-counts-export.py):
gs://2024-wgspd/snv/coding/20240625_subset_post-qc_protein-coding_VEP-annotated_original.ht
"""
import hail as hl
from gnomad.utils.vep import vep_or_lookup_vep
from gnomad.utils.vep import process_consequences

hl.init(default_reference = 'GRCh38',
                tmp_dir = "gs://wes-bipolar-tmp-4day/")

# Some hardcoded cases
PLOF_CSQS = ["transcript_ablation", "splice_acceptor_variant",
            "splice_donor_variant", "stop_gained", "frameshift_variant"]
MISSENSE_CSQS = ["stop_lost", "start_lost", "transcript_amplification",
                "inframe_insertion", "inframe_deletion", "missense_variant",
                "splice_region_variant"] # Maybe remove splice_region_variant
SYNONYMOUS_CSQS = ["stop_retained_variant", "synonymous_variant"]
OTHER_CSQS = ["coding_sequence_variant", "mature_miRNA_variant", "5_prime_UTR_variant",
            "3_prime_UTR_variant", "non_coding_transcript_exon_variant", "intron_variant",
            "NMD_transcript_variant", "non_coding_transcript_variant", "upstream_gene_variant",
            "downstream_gene_variant", "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant",
            "regulatory_region_ablation", "regulatory_region_amplification", "feature_elongation",
            "regulatory_region_variant", "feature_truncation", "intergenic_variant"]
POSSIBLE_LOFTEE_FLAGS = ["SINGLE_EXON",
                         "NAGNAG_SITE",
                         "PHYLOCSF_WEAK", "PHYLOCSF_UNLIKELY_ORF",
                         "NON_CAN_SPLICE"]
# Updated case builder 
# included_loftee_flags should be a list of loftee flags that we will tolerate (not automatically counted as LC)
def annotation_case_builder_updated(worst_csq_for_variant_canonical_expr, 
                                    lof_use_loftee: bool = True, 
                                    included_loftee_flags: list[str] = None,
                                    mis_use_polyphen_and_sift: bool = False, 
                                    mis_use_strict_def: bool = False, 
                                    syn_use_strict_def: bool = False):
    case = hl.case(missing_false=True)
    if lof_use_loftee:
        if included_loftee_flags is None: # If included_loftee_flags is None (WE DON'T WANT TO CONSIDER ANY FLAGS)
            case = (case
                    .when(worst_csq_for_variant_canonical_expr.lof == 'HC', 'pLoF') 
                    .when(worst_csq_for_variant_canonical_expr.lof == 'LC', 'LC'))
        elif len(included_loftee_flags) == 0: # If included_loftee_flags is empty (WE WANT TO FILTER OUT ALL FLAGS)
            case = (case
                    .when(worst_csq_for_variant_canonical_expr.lof == 'HC', 'pLoF') #predicted loss-of-function
                    .when((worst_csq_for_variant_canonical_expr.lof == 'LC') | 
                          (hl.is_defined(worst_csq_for_variant_canonical_expr.lof_flags)), 'LC')) # Low confidence if labelled so, or if there is an additional flag
        elif len(included_loftee_flags) > 0: # If included_loftee_flags has values (WE WANT TO TOLERATE SOME FLAGS)
             case = (case
                    .when(worst_csq_for_variant_canonical_expr.lof == 'HC', 'pLoF') 
                    .when((worst_csq_for_variant_canonical_expr.lof == 'LC') | 
                          (hl.is_defined(worst_csq_for_variant_canonical_expr.lof_flags) & ~hl.set(included_loftee_flags).contains(worst_csq_for_variant_canonical_expr.lof_flags)), 'LC')) # Low confidence if labelled so, or if there is an additional BAD flag
    else:
        case = case.when(hl.set(PLOF_CSQS).contains(worst_csq_for_variant_canonical_expr.most_severe_consequence), 'pLoF')
    if mis_use_polyphen_and_sift:
        case = (case
                .when(hl.set(MISSENSE_CSQS).contains(worst_csq_for_variant_canonical_expr.most_severe_consequence) &
                      (worst_csq_for_variant_canonical_expr.polyphen_prediction == "probably_damaging") &
                      (worst_csq_for_variant_canonical_expr.sift_prediction == "deleterious"), "damaging_missense")
                .when(hl.set(MISSENSE_CSQS).contains(worst_csq_for_variant_canonical_expr.most_severe_consequence), "other_missense"))
    else:
        if mis_use_strict_def:
            case = case.when(worst_csq_for_variant_canonical_expr.most_severe_consequence == 'missense_variant', 'missense')
        else:
            case = case.when(hl.set(MISSENSE_CSQS).contains(worst_csq_for_variant_canonical_expr.most_severe_consequence), 'missense')
    if syn_use_strict_def:
        case = case.when(worst_csq_for_variant_canonical_expr.most_severe_consequence == 'synonymous_variant', 'synonymous')
    else:
        case = case.when(hl.set(SYNONYMOUS_CSQS).contains(worst_csq_for_variant_canonical_expr.most_severe_consequence), 'synonymous')
    case = case.when(hl.set(OTHER_CSQS).contains(worst_csq_for_variant_canonical_expr.most_severe_consequence), 'non_coding')
    return case.or_missing()


# Read MT
MT = 'gs://2024-wgspd/snv/coding/202240618_subset_post-qc_protein-coding.mt'
mt = hl.read_matrix_table(MT)
print(f"Original data count: {mt.count()}")
#Original data count: (144360038, 28554)

# VEP
ht = mt.rows()
ht_vep = vep_or_lookup_vep(ht, reference="GRCh38")
ht_vep = ht_vep.select('vep')

# Process consequences
ht_vep = process_consequences(ht_vep)
ht_vep = ht_vep.annotate(consequence_category = annotation_case_builder_updated(ht_vep.vep.worst_csq_for_variant_canonical, 
                                                                                False, None, 
                                                                                True, False, False)) # Matching original analysis for now


# Add SNP, indel, none; inframe indel
ht_vep = ht_vep.annotate(
    type = (hl.case()
            .when((hl.len(ht_vep.alleles[0]) == 1) & (hl.len(ht_vep.alleles[1]) == 1), "SNP")
            .when(hl.len(ht_vep.alleles[0]) < hl.len(ht_vep.alleles[1]), "Insertion")
            .when(hl.len(ht_vep.alleles[0]) > hl.len(ht_vep.alleles[1]), "Deletion")
            .or_missing()),
    infrIndel = ((ht_vep.vep.worst_csq_for_variant_canonical.most_severe_consequence == "inframe_insertion") | 
                 (ht_vep.vep.worst_csq_for_variant_canonical.most_severe_consequence == "inframe_deletion"))
)

ht_vep = ht_vep.repartition(500, shuffle = False)

# Write
ht_vep.write("gs://2024-wgspd/snv/coding/20240716_subset_post-qc_protein-coding_VEP-annotated_original_all-types.ht", overwrite = True)


