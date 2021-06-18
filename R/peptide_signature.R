generate_substrates <- function(fisher_tables, uniprot_dt, screener_dt, target_kinase, 
                               method = c("prod", "log2_sum", "w_prod"), 
                               screening_kinase = "ALL", n_hits = NULL){
  if (is.null(n_hits)){
    n_hits = 10
  }
  cat("Identifying candidate amino acids...")
  fisher_candidates <- get_candidateaa(fisher_tables, screener_dt, target_kinase)
  candidate_matrix <- score_candidates(fisher_candidates, screener_dt, method, 
                                       screening_kinase)
  spec_barcodes <- candidate_matrix[[1]][kinase == target_kinase & n_active == 
                                        min(n_active) & active == "TRUE"][, 
                                        substrate_barcode]
  
  top_barcodes <- candidate_matrix[[1]][substrate_barcode %in% spec_barcodes & 
                                       kinase != target_kinase][, sum(score), 
                                       by = substrate_barcode][order(V1)][1:n_hits, 
                                       substrate_barcode]
  
  output <- list(candidate_aa = fisher_candidates, candidate_matrix[[1]], 
                 candidate_matrix[[2]], top_hits = candidate_matrix[[2]][
                 substrate_barcode %in% top_barcodes])
  return(output)
}

get_candidateaa <- function(fisher_tables, screener_dt, target_kinase){
  # aa_cols <- c(paste0("-",rev(seq(1:7))),"0", paste0(seq(1:7)))
  ptm_pos <- paste0(names(which.max(fisher_tables[[2]][,8])),":0")
  
  valid_cols <- which(apply(fisher_tables[[2]], 2, function(x) 
                     length(unique(x)) != 1)) #remove columns where there is no data
  kinase_fisher <- lapply(fisher_tables, function(x) x[, valid_cols])
  
  fisher_combn <- fisher_long(kinase_fisher)
  fisher_combn[, sig_pos:= sum(fisher_pval <= 0.05), by=flank_pos]
  
  sig_fav_pssm <- fisher_combn[fisher_pval <= 0.05 & fisher_odds > 1 & 
                              sig_pos >0][, .(barcode, flank_pos)]
  
  sig_disfav_pssm <- fisher_combn[fisher_pval <= 0.05 & fisher_odds < 1 & 
                                 sig_pos >0][, .(barcode, flank_pos)]
  
  if (target_kinase %in% screener_dt[[2]][, kinase]){
    sig_fav_screener <- screener_dt[[2]][kinase == target_kinase & 
                                        fisher_pval <= 0.05 & 
                                        fisher_odds > 1][, .(
                                        barcode, flank_pos)]
    sig_disfav_screener <- screener_dt[[2]][kinase == target_kinase & 
                                           fisher_pval <= 0.05 & 
                                           fisher_odds < 1][, .(
                                           barcode, flank_pos)]
    
    sig_fav_all <- unique(rbind(sig_fav_pssm, sig_fav_screener))
    sig_disfav_all <- unique(rbind(sig_disfav_pssm, sig_disfav_screener))
  }else{
    sig_fav_all <- sig_fav_pssm
    sig_disfav_all <- sig_disfav_pssm
  }
  
  aa_spec <- setdiff(aa_cols, c(sig_fav_all[, flank_pos], "0"))
  
  disfav_screener <- screener_dt[[2]][kinase != target_kinase & flank_pos %in% 
                                     aa_spec & fisher_odds < 1][, .(
                                     length(fisher_odds), mean(fisher_odds, 
                                     na.rm = T)), by = .(barcode, flank_pos)]
  
  if (target_kinase %in% screener_dt[[2]][, kinase]){
    kinase_crossref <- screener_dt[[2]][kinase == target_kinase & barcode %in% 
                                       disfav_screener[, barcode] & 
                                       fisher_odds < 1][, barcode]
    disfav_screener <- disfav_screener[!barcode %in% kinase_crossref]
  }
  
  disfav_screener <- disfav_screener[, .SD[V1 %in% rev(sort(V1))][1:5], 
                                    by = flank_pos]
  
  fav_pssm <- fisher_combn[barcode %in% disfav_screener[, barcode] & 
                             fisher_odds > 1 & !barcode %in% sig_disfav_all | 
                             barcode %in% sig_fav_all[, barcode]]
  
  ptm_pssm <- fisher_combn[barcode == ptm_pos]
  
  aa_candidates <- rbind(fav_pssm, ptm_pssm)
  
  aa_nomatch <- setdiff(aa_cols, aa_candidates[, flank_pos])
  
  if (length(aa_nomatch) > 0){
    pssm_nomatch <- fisher_combn[flank_pos %in% aa_nomatch][, 
                                fisher_odds*1/fisher_pval, by = .(
                                amino_acid, flank_pos, barcode)][, .SD[
                                V1 %in% rev(sort(V1))[1:3]], by = .(
                                flank_pos)][, barcode]
    aa_candidates <- rbind(aa_candidates, fisher_combn[barcode %in% pssm_nomatch])
  }
  return(aa_candidates)
}

score_candidates <- function(fisher_candidates, screener_dt, method = 
                            c("prod", "log2_sum", "w_prod"), kinase = "ALL"){
  candidate_aa <- split(fisher_candidates[, amino_acid], fisher_candidates[, 
                     flank_pos])
  
  cat("\nCreating substrate permutations...")
  new_substrates <- do.call(expand.grid,candidate_aa)
  
  new_substrates <- reshape2::melt(t(new_substrates), varnames = c("flank_pos", 
                                  "substrate_barcode"), value.name = 
                                  "amino_acid")
  new_substrates <- data.table(new_substrates)
  new_substrates[, barcode:= paste0(amino_acid, ":", flank_pos)]
  new_substrates<- merge(fisher_candidates[, .(fisher_odds), by = barcode], 
                         new_substrates, by = "barcode")
  
  cat("\nScoring candidate sequences...")
  
  if (method == "w_prod"){
    new_substrates[, raw_score:= get_score(fisher_odds, method, fisher_pval), 
                   by = substrate_barcode]
  }else{
    new_substrates[, raw_score:= get_score(fisher_odds, method), 
                   by = substrate_barcode]
  }
  
  substrate_score <- multi_candidate_screener(screener_dt, new_substrates, 
                                              kinase, method, FALSE)
  
  n_hits <- unique(substrate_score[, n_active, by = substrate_barcode])
  n_hits[, substrate_barcode:= as.numeric(substrate_barcode)]
  
  candidate_matrix<- data.table(reshape2::dcast(new_substrates, 
                               substrate_barcode + raw_score ~ flank_pos, 
                               value.var = "amino_acid"))
  candidate_matrix <- candidate_matrix[rev(order(raw_score))]
  
  candidate_matrix <- merge(candidate_matrix, n_hits, by = "substrate_barcode")
  output <- list(substrate_score, candidate_matrix)
}

get_similarity <- function(aa_seq1, aa_seq2){
  common_pos <- intersect(colnames(aa_seq1), colnames(aa_seq2))
  aa_seq1 <- aa_seq1[, ..common_pos]
  aa_seq2 <- aa_seq2[, ..common_pos]
  
  y<- sapply(aa_seq1, function(x) which(colnames(blosum62) == x))
  x<- sapply(aa_seq2, function(x) which(colnames(blosum62) == x))
  idx <- cbind(y,x)
  similarity_score <- sapply(1:nrow(idx), function(i) 
                            blosum62[idx[i,1], idx[i, 2]])
  
  cat(paste0("The similarity scores between sequences:\n", 
      apply(aa_seq1, 1, paste0, collapse=""), " and ", apply(aa_seq2, 1, 
      paste0, collapse=""), "\nare "), paste(similarity_score))
}

get_normscores <- function(candidates, screener_dt, n_hits = NULL){
  if (is.null(n_hits)){
    n_hits = 10
  }
  top_hits <- candidates[[3]][1:n_hits,]
  sub_score <- candidates[[2]][substrate_barcode %in% 
                              top_hits[, substrate_barcode]]
  normscores <- norm_scores(sub_score, screener_dt)
}

norm_scores <- function(scoretable, screener_dt){
  score_quantile <- screener[[3]]
  
  map_table <- data.table(merge(scoretable, score_quantile, 
                         by = c("kinase", "cutpoint"), all = T))
  map_table[, converted_score:= round((log(score) - log(min)) / (log(max) - 
                                     log(min))*(100), 0), by = kinase]
  scaled_hits <- reshape2::dcast(map_table, substrate_barcode ~ kinase, 
                       value.var = "converted_score")
  scaled_threshold <- data.table(screener[[3]])[,round((log(cutpoint) - 
                                log(min))/(log(max) - log(min))*(100),0), 
                                by = kinase]
  setnames(scaled_threshold, "V1", "norm_thresholds")
  output <- list(scaled_hits, scaled_threshold)
}
