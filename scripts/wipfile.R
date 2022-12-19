

for(i in 1:nrow(df_odds_filt[1:2,])){
  gene_symbol <- df_odds_filt$gene_symbol[i]
  peptide_id <- df_odds_filt$peptide_id[i]
  psequence <- df_odds_filt$psequence[i]
  topology <- df_odds_filt$topology[i]
  print(paste("running:", gene_symbol))
  
  ITIM_locs <- str_locate_all(psequence, ITIMseq)[[1]]
  inner_ITIMs <- ITIM_locs[str_sub(topology, ITIM_locs[,"start"],ITIM_locs[,"end"]) == "iiiiii",]
  
  find_pdb <- str_detect(pdb_files, peptide_id)
  if(any(find_pdb)){
    pdb <- read.pdb(pdb_files[find_pdb])
    
    conf_scores <- pdb$atom %>%
      distinct(resno, b) %>%
      pull(b)
    
    if(length(inner_ITIMs) == 0){
      print("no ITIMs?")
      break()
    }
    
    if (!length(inner_ITIMs) == 2){
      print(paste("nr of ITIMs:", length(inner_ITIMs)/2))
      means <- c()
      for (j in 1:nrow(inner_ITIMs)){
        print(paste("j:", j))
        means <- c(means, mean(conf_scores[inner_ITIMs[j, "start"]:inner_ITIMs[j, "end"]]))
      }
      confidence_scores[[i]] <- means
    } else if(length(inner_ITIMs) == 2) {
      print(paste("nr of ITIMs:", length(inner_ITIMs)/2))
      confidence_scores[[i]] <- mean(conf_scores[inner_ITIMs["start"]:inner_ITIMs["end"]])
    }
    print(paste("confidence score(s):", confidence_scores[[i]]))
    df_AF[[i]] <- confidence_scores[[i]]
    names(df_AF)[i] <- df_odds_filt$peptide_gene[i]
    
    print(paste(i,"/",nrow(df_odds_filt),"sequences done - ", gene_symbol))
    
  }else{
    print(".pdb file not found")
    df_AF[[i]] <- NA
    names(df_AF)[i] <- df_odds_filt$peptide_gene[i]
    #break()
  }
}

df_AF_all <- df_AF %>%
  append(list(1:5)) %>%  
  as_tibble_col("ITIM") %>% 
  unnest_wider(ITIM,names_sep = "") %>% 
  mutate(peptide_gene = c(names(df_AF),1))
