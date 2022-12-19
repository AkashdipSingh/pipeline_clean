
## check how many peptides correspond to how many genes
#for each individual step

annot_allpeps <- getBM(attributes = c("ensembl_peptide_id","hgnc_symbol"),
                       filters = "ensembl_peptide_id",
                       values = peptides$ensembl_peptide_id,
                       mart = ensembl)

pipeline_steps <- annot_allpeps %>%
  mutate(ITIMcontaining = ifelse(ensembl_peptide_id %in% peptides_filt$peptide_id, TRUE, FALSE),
         membraneprot = ifelse(ensembl_peptide_id %in% topologies$peptide_id, TRUE, FALSE),
         innerITIM = ifelse(ensembl_peptide_id %in% topologies_filt$peptide_id, TRUE, FALSE),
         uniquegene = !duplicated(hgnc_symbol))

pipeline_steps %>% 
  filter(uniquegene,
         any(ITIMcontaining,membraneprot,innerITIM) == TRUE) %>% 
  distinct(hgnc_symbol, ITIMcontaining, membraneprot, innerITIM) %>% 
  group_by(ITIMcontaining, membraneprot, innerITIM) %>% 
  summarise(n(unique))


annot_allpeps %>%
  filter(ensembl_peptide_id %in% peptides_filt$peptide_id) %>%
  pull(hgnc_symbol) %>%
  unique() %>% 
  length()

annot_allpeps %>%
  filter(ensembl_peptide_id %in% topologies$peptide_id) %>%
  pull(hgnc_symbol) %>%
  unique() %>%
  length()

annot_allpeps %>%
  filter(ensembl_peptide_id %in% topologies_filt$peptide_id) %>%
  pull(hgnc_symbol) %>%
  unique() %>% 
  length()

check_genes <- c("TNF","SYNCRIP", "PCDH11Y") 

genes_noAF <- read_lines("../pipeline_paper_analyses/output/220829/genes_noITIM_AlphaFold.txt")

peptides_noAF <- annot_allpeps %>%
  filter(hgnc_symbol %in% genes_noAF)

isoform_sequences_noAF <- topologies_filt %>%
  filter(peptide_id %in% peptides_noAF$ensembl_peptide_id) %>%
  left_join(peptides_noAF, by = c("peptide_id" = "ensembl_peptide_id"))

write_fasta <- function(x, y){
  fasta <- c()
  for (i in 1:nrow(x)){
    print(i)
    fasta = c(fasta, as.character(paste(">", x[i,"peptide_id"], sep = "")))
    fasta = c(fasta,as.character(x[i,"psequence"]))
  }
  write_lines(fasta, y)
}

isoform_sequences_noAF %>%
  mutate(peptide_id = paste(peptide_id, hgnc_symbol, sep = "_")) %>%
  dplyr::select(peptide_id, psequence) %>% 
  write_fasta("test.fasta")
