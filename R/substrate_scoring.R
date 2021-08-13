multi_candidate_screener <- function(screener_dt, 
                                     candidates_dt, 
                                     kinase, 
                                     family = FALSE){
  
  method <- screener_dt[[4]]$method
  pval_corr <- screener_dt[[4]]$pval_corr
  norm_method <- screener_dt[[4]]$norm_method
  
  n_candidates <- unique(candidates_dt[, substrate_barcode])
  output <- list()
  pb <- txtProgressBar(min = 0, 
                       max = length(n_candidates), 
                       initial = 0, 
                       style = 3)
  
  for (i in 1:length(n_candidates)){
    output[[i]] <- peptide_screener(screener_dt, 
                                    candidates_dt[substrate_barcode == 
                                                  n_candidates[i]][, barcode],
                                    kinase,
                                    family)
    setTxtProgressBar(pb, i)
  }
  names(output) <- paste("substrate_barcode", n_candidates) 
  candidate_hits <- Filter(Negate(is.null), output)
  names <- gsub("substrate_barcode ", "", names(candidate_hits))
  candidate_dt <- data.table(reshape2::melt(setNames(candidate_hits, names), 
                            id.vars = c("kinase", 
                                        "active", 
                                        "perf", 
                                        "score", 
                                        "cutpoint")))
  colnames(candidate_dt)[length(colnames(candidate_dt))] <- "substrate_barcode"
  candidate_dt[, n_active:= sum(active), by = substrate_barcode]
  return(candidate_dt)
}

peptide_screener <- function(screener_dt, 
                             candidates_dt, 
                             target_kinase, 
                             family = FALSE){
  
  method <- screener_dt[[4]]$method
  pval_corr <- screener_dt[[4]]$pval_corr
  norm_method <- screener_dt[[4]]$norm_method
  
  scores_dt <- screener_dt[[2]]
  if (isTRUE(pval_corr)){
    scores_dt[, fisher_odds:= ifelse(fisher_pval > 0.05, 1, fisher_odds)]
  }
  sub_score <- scores_dt[barcode %in% candidates_dt]
  if (method == "w_prod"){
    sub_score[, score:= get_score(fisher_odds, method, fisher_pval), by = kinase]
  }else{
    sub_score[, score:= get_score(fisher_odds, method), by = kinase]
  }
  
  scores_quantiles <- screener_dt[[3]]
  scores_merge <- data.table(merge(scores_quantiles, 
                                   sub_score[, score, by = kinase],
                                   by = "kinase"))
  if (norm_method == "bkgrnd"){
    scores_merge[, score:= (score - bkgrnd_mean) / bkgrnd_sd]
  }  
  
  scores_merge[, active:= ifelse(score >= cutpoint, TRUE, FALSE)]
  scores_merge[, perf:= ifelse(score > Q90, "high", 
                       ifelse(score < cutpoint, "inactive", 
                       ifelse(score < Q10, "low", "medium")))]
  
  candidate_results <- unique(scores_merge[, .(score, cutpoint, active, perf), 
                                          by = kinase])  
  if (!"ALL" %in% target_kinase & !isTRUE(family)){
    candidate_results <- candidate_results[kinase %in% target_kinase]
  } else if (isTRUE(family)){
    kinase_family <- kinase_anno[Name %in% target_kinase | 
                                GENENAME %in% target_kinase][, 
                                Family]
    family_set <- unique(unlist(kinase_anno[Family %in% kinase_family][, .(GENENAME, Name)]))
    candidate_results <- candidate_results[kinase %in% family_set]
  }
  return(candidate_results)
}


get_score <- function(odds, method = c("prod", "log2_sum", "w_prod"), pval = NULL){
  method <- match.arg(method)
  if (method == "prod"){
    score <- prod(odds, na.rm = T)
  }else if (method == "log2_sum"){
    pseudo_count = 1
    score <- sum(log(odds), na.rm = T)
  }else{
    # pval_inv <- 1/pval
    # log_pval <- log2(pval_inv)/sum(log2(pval_inv))
    # w_odds <- mapply(function(x) odds[x]^log_pval[x], seq_along(log_pval))
    w_odds <- mapply(function(x) odds[x]^exp(-pval[x]), seq_along(pval))
    
    score <- prod(w_odds)
    # matrix <- data.table(odds, pval)
    # sig <- matrix[pval <= 0.05][, prod(odds)]
    # nonsig <- matrix[pval > 0.05][, prod(odds)]
    # score <- (sig^0.8)*(nonsig^0.2)
  }
  return(score)
}

geo_mean <- function(data){
  log_data <- log(data)
  gm <- exp(mean(log_data[is.finite(log_data)]))
  return(gm)
}

fisher_long <- function(fisher_dt, n = NULL){
  if (is.null(n)){
    pval_long <- data.table(reshape2::melt(fisher_dt[[1]], 
                                           varnames=c("amino_acid", "flank_pos"),
                                           value.name="fisher_pval")
    )
    odds_long <- data.table(reshape2::melt(fisher_dt[[2]], 
                                           varnames=c("amino_acid", "flank_pos"),
                                           value.name="fisher_odds")
    )
    merge_dt <- merge(odds_long, pval_long, by = c("amino_acid", "flank_pos")) 
  }else{
    merge_dt <- data.table(reshape2::melt(fisher_dt, 
                                          varnames=c("amino_acid", "flank_pos"), 
                                          value.name="fisher_odds"))
  }
  
  merge_dt[,barcode:= paste0(amino_acid, ":", flank_pos)]
  merge_dt[, fisher_odds:= ifelse(fisher_odds == 0 & 
                                    fisher_pval > 0.05, NA,
                                  fisher_odds)]
  
  return(merge_dt)
}
