substrate_fisher_test <- function(substrates_dt, uniprot_dt, type = 
                                 c("aa", "aa_property"), property = 
                                 c("property_chemical", "property_hydropathy", 
                                 "property_volume", "property_polar")){
  type <- match.arg(type)
  property <- match.arg(property)
  
  if (type == "aa"){
    aa_rows <- aa_freq(substrates_dt)[[1]]
    # aa_cols <- c(paste0("-", rev(seq(1:7))), "0", paste0(seq(1:7)))
    sub_freq <- aa_freq(substrates_dt)
    uniprot_freq <- aa_freq(uniprot_dt)
    uniprot_freq <- uniprot_freq[aa %in% sub_freq[, aa]]
  }else{
    # aa_cols <- c(paste0("-", rev(seq(1:7))), "0", paste0(seq(1:7)))
    sub_freq <- aa_property_freq(substrates_dt, property)
    aa_rows <- sub_freq[[1]]
    uniprot_freq <- aa_property_freq(uniprot_dt, property)
    uniprot_freq <- uniprot_freq[eval(as.name(property)) %in% 
                                 sub_freq[, ..property][[1]]]
  }
  
  sub_freq <- sub_freq[, -1]
  uniprot_freq <- uniprot_freq[, -1]
  
  sub_row_sums <- colSums(sub_freq)
  uniprot_row_sums <- colSums(uniprot_freq)
  
  fisher_pval <- matrix(nrow=nrow(sub_freq), ncol = ncol(sub_freq))
  colnames(fisher_pval) <- aa_cols
  rownames(fisher_pval) <- aa_rows
  
  fisher_odds <- matrix(nrow = nrow(sub_freq), ncol = ncol(sub_freq))
  colnames(fisher_odds) <- aa_cols
  rownames(fisher_odds) <- aa_rows 
  
  for (i in 1:ncol(sub_freq)){
    for (j in 1:nrow(sub_freq)){
      fisher_matrix <- matrix(unlist(c(sub_freq[j, ..i], 
                                       sub_row_sums[i], 
                                       uniprot_freq[j, ..i], 
                                       uniprot_row_sums[i])), ncol = 2)
      fisher_results <- fisher.test(fisher_matrix)
      fisher_pval[j, i] <- fisher_results$p.value
      fisher_odds[j, i] <- fisher_results$estimate
    }
  }
  fisher_odds[fisher_odds == 0] <- 0.05
  fisher_list <- list(fisher_pval,fisher_odds)
  names(fisher_list) <- c("fisher_p_values", "fisher_odds")
  return(fisher_list)
}

aa_freq <- function(substrates_dt){
  # aa_cols <- c(paste0("-",rev(seq(1:7))),"0", paste0(seq(1:7)))
  pseudo_count <- 1
  substrate_melt <- unique(data.table(melt(substrates_dt,id.vars = 
                          c("substrate_barcode"), measure.vars = aa_cols, 
                          value.name = "aa")))
  substrate_cast <- data.table::dcast(substrate_melt,aa~variable, value.var = "aa", 
                         fun.aggregate = length)
  substrate_cast <- substrate_cast[!aa %in% c("NA", " ", "-")]
  substrate_cast <- substrate_cast[!is.na(aa)]
  
  return(substrate_cast)
}

aa_property_freq <- function(substrates_dt, property = c("property_chemical", 
                            "property_hydropathy", "property_volume", 
                            "property_polar")){
  property <- match.arg(property)
  
  # aa_cols <- c(paste0("-",rev(seq(1:7))),"0", paste0(seq(1:7)))
  substrate_melt<- unique(data.table(melt(substrates_dt,id.vars = 
                         c("substrate_barcode"), measure.vars = aa_cols, 
                         value.name = "aa")))
  substrate_melt<- merge(substrate_melt,aa_classification, by="aa", all.x=T)
  substrate_cast<- data.table::dcast(substrate_melt, 
                                   eval(as.name(property)) ~ variable,
                                   fun.aggregate = length,
                                   value.var = property)
  colnames(substrate_cast)[[1]] <- property
  substrate_cast<- substrate_cast[!is.na(eval(as.name(property)))]
  return(substrate_cast)
}
