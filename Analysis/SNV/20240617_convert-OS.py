import hail as hl

hl.init(gcs_requester_pays_configuration = 'wes-bipolar',
        tmp_dir = 'gs://wes-bipolar-tmp-4day')

"""
ht = hl.read_table("gs://gcp-public-data--gnomad/resources/context/grch37_context_vep_annotated.ht/")

ht.describe()
print(ht.count()) # 8575974294
ht = ht.annotate(OS = ht.vep.transcript_consequences.lof.contains("OS"))


rg37 = hl.get_reference('GRCh37')  
rg38 = hl.get_reference('GRCh38')  
rg37.add_liftover('gs://hail-common/references/grch37_to_grch38.over.chain.gz', rg38)  

ht = ht.annotate(new_locus=hl.liftover(ht.locus, 'GRCh38'))
h = ht.filter(ht.OS == True).persist() # Filter to only those with OS annotation
# gs://wes-bipolar-tmp-4day/persist_TableGL21x2S2M8

h = h.filter(hl.is_defined(h.new_locus)) # Filter to only those with successful liftover
print(h.count()) # 931684

h = h.annotate(GRCh37_locus = h.locus, GRCh38_locus = h.new_locus)

h = h.annotate(transcript_id = h.vep.transcript_consequences.transcript_id,
               gene_id = h.vep.transcript_consequences.gene_id,
               gene_symbol = h.vep.transcript_consequences.gene_symbol)
h = h.drop(h.context, h.vep, 
           h.a_index, h.was_split, 
           h.methylation, h.coverage,
           h.gerp)
           
h = h.checkpoint("gs://bipex2/annotations/20240617_GRCh38_OS.ht", overwrite = True)
"""
ht = hl.read_table("gs://bipex2/annotations/20240617_GRCh38_OS.ht")

h = ht.repartition(500, shuffle = False)
h = h.key_by(locus = h.GRCh38_locus, alleles = h.alleles)

#h = h.checkpoint("gs://wes-bipolar-tmp-4day/20240618_GRCh38_OS_reparititioned.ht", overwrite = True)
h = h.checkpoint("gs://gs://bipex2/annotations/20240618_GRCh38_OS.ht")
#ht = ht.key_by(locus=ht.new_locus, alleles = ht.alleles) 
#ht = ht.select('OS')


