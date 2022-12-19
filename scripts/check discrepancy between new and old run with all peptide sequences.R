df_odds <- readRDS(paste0(processedF, "df_odds",n_perm,".rds"))

df_annot <- readRDS(paste0("../pipeline_paper_analyses/processed/220829/220829_10000_ITIModds.rds")) %>% 
  dplyr::select(gene_symbol, totalcounts2,odds_filt,TMtype)

df_inners <-  readRDS(paste0("../pipeline_paper_analyses/processed/220829/220829_localization_annot_all.rds")) %>%
  mutate(check_TM = ifelse(is.na(TM_localisation),FALSE,TRUE)) %>% 
  dplyr::select(gene_symbol, check_TM)

df_old <- df_inners %>% 
  left_join(df_annot) %>% 
  mutate(noITIM = ifelse(totalcounts2 == 1, TRUE, FALSE))



df_new <- df_odds %>% 
  dplyr::select(peptide_id, gene_symbol,total_odds, odds_filt)

all(df_old$gene_symbol %in% df_new$gene_symbol)
all(df_new$gene_symbol %in% df_old$gene_symbol)

setdiff(df_old$gene_symbol, df_new$gene_symbol)
setdiff(df_new$gene_symbol, df_old$gene_symbol) %in% unlist(dropouts)

dropouts <- readRDS(paste0("../pipeline_paper_analyses/processed/220829/220829_dropouts.rds"))


no_odds_old <- df_old %>%
  filter(noITIM | !check_TM) %>% 
  pull(gene_symbol) %>% 
  unique() %>% 
  c(unlist(dropouts))

odds_old <- df_old %>% 
  filter(odds_filt == "incl") %>% 
  pull(gene_symbol) %>% 
  unique()
  
new_odds_genes <- df_new %>% 
  filter(odds_filt == "incl") %>% 
  pull(gene_symbol) %>% 
  unique()

new_genes <- setdiff(new_odds_genes,c(odds_old,no_odds_old))

df_old %>% 
  filter(gene_symbol %in% new_genes)

df_new %>%
  filter(gene_symbol %in% new_genes) %>% 
  arrange(gene_symbol) %>% 
  group_by(gene_symbol, odds_filt) %>% 
  summarise(n()) %>% 
  print(n = 200)

new_genes_oneseq <- df_new %>%
  filter(gene_symbol %in% new_genes) %>% 
  arrange(gene_symbol) %>% 
  group_by(gene_symbol) %>%
  mutate(check = all(odds_filt == "incl")) %>% 
  filter(check) %>% 
  pull(gene_symbol) %>% 
  unique()

df_old %>% 
  filter(gene_symbol %in% new_genes_oneseq)

df_new %>%
  filter(gene_symbol %in% new_genes) %>% 
  arrange(gene_symbol) %>% 
  group_by(gene_symbol) %>%
  mutate(check = all(odds_filt == "incl")) %>% 
  filter(check)

df_odds_clean <- readRDS(paste0(processedF, "df_odds",n_perm,".rds")) %>%
  dplyr::select(peptide_id, gene_symbol, psequence_new = psequence, topology_new = topology)
  

readRDS(paste0("../pipeline_paper_analyses/processed/220829/220829_10000_ITIModds.rds")) %>% 
  filter(gene_symbol %in% new_genes_oneseq) %>% 
  dplyr::select(gene_symbol,psequence_old = psequence, topology_old = topology) %>% 
  left_join(df_odds_clean) %>% 
  mutate(clean = psequence_old == psequence_new,
         clean2 = topology_old == topology_new)

