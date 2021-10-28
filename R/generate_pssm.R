#' Construct PSSM and perform Fisher Exact Test
#'
#' @description This function uses the Fisher Exact Test to calculate the 
#' odds ratio and significance of each residue across the flanking sequence.
#'  
#' @param substrates_dt An enzymatic preference data.table. 
#' 
#' @param uniprot_dt A data.table containing in silico negative peptide 
#' sequences for each `uniprot_id` found in the `substrates_dt` file.
#' 
#' @param type Character of either `aa` or `aa_property` to indicate whether 
#' amino acid residues or amino acid properties should be used.
#' 
#' @param property Character indicating which amino acid property should be 
#' analyzed. Only used if type = `aa`. See available property options 
#' by typing `aa_classificiation`.
#'
#' @return A list containing two data tables: 1) Exact Test P-values and 
#' 2) Fisher Exact Odds
#' @export
#'
#' @examples
#' pssm <- substrate_fisher_test(substrates_dt = uniq_subs, uniprot_dt = uniprot,
#' type = "aa")
substrate_fisher_test <- function(substrates_dt, uniprot_dt, 
                                  type = c("aa", "aa_property"), 
                                  property = c("property_chemical", 
                                               "property_hydropathy", 
                                               "property_volume", 
                                               "property_polar", 
                                               "property_charge")){
  type <- match.arg(type)
  
  if (type == "aa_property"){
    property <- match.arg(property)
  }
  
  if (type == "aa"){
    aa_rows <- aa_freq(substrates_dt)[[1]]
    sub_freq <- aa_freq(substrates_dt)
    uniprot_freq <- aa_freq(uniprot_dt)
    
    #Match background AAs with those in data file 
    uniprot_freq <- uniprot_freq[aa %in% sub_freq[, aa]]
  }else{
    sub_freq <- aa_property_freq(substrates_dt, property)
    aa_rows <- sub_freq[[1]]
    uniprot_freq <- aa_property_freq(uniprot_dt, property)
    uniprot_freq <- uniprot_freq[eval(as.name(property)) %in% 
                                 sub_freq[, ..property][[1]]]
  }
  
  #Subset numeric portion of frequency tables
  sub_freq <- sub_freq[, -1]
  uniprot_freq <- uniprot_freq[, -1]
  
  sub_row_sums <- colSums(sub_freq)
  uniprot_row_sums <- colSums(uniprot_freq)
  
  #Initialize PSSM tables
  fisher_pval <- matrix(nrow = nrow(sub_freq), ncol = ncol(sub_freq))
  colnames(fisher_pval) <- aa_cols
  rownames(fisher_pval) <- aa_rows
  
  fisher_odds <- matrix(nrow = nrow(sub_freq), ncol = ncol(sub_freq))
  colnames(fisher_odds) <- aa_cols
  rownames(fisher_odds) <- aa_rows 
  
  for (i in 1:ncol(sub_freq)){
    for (j in 1:nrow(sub_freq)){
      #Set up 2x2 contingency table
      fisher_matrix <- matrix(unlist(c(sub_freq[j, ..i], 
                                       sub_row_sums[i], 
                                       uniprot_freq[j, ..i], 
                                       uniprot_row_sums[i])), 
                              ncol = 2)
      
      colnames(fisher_matrix)<-c("expt", "bkgrnd")
      rownames(fisher_matrix)<-c("frequency", "total")
      
      fisher_results <- fisher.test(fisher_matrix)
      
      fisher_pval[j, i] <- fisher_results$p.value
      
      if(is.finite(fisher_results$estimate)){
        fisher_odds[j, i] <- fisher_results$estimate
      }else{
        fisher_odds[j, i] <- fisher_results$conf.int[[1]]
      }
    }
  }
  
  #Replace exact 0 odds value with small pseudo-odds 
  fisher_odds[fisher_odds == 0] <- 0.05
  
  fisher_list <- list(fisher_pval,fisher_odds)
  
  names(fisher_list) <- c("fisher_p_values", "fisher_odds")
  
  return(fisher_list)
}


#' Calculate frequency of each amino acid residue within PSSM
#'
#' @description This lower-level function is called within `substrate_fisher_test`. 
#' End users should use the higher-level function instead.
#' 
#' @param substrates_dt An enzymatic preference data.table. 
#'
#' @return A data.table object
#' @export
#'
aa_freq <- function(substrates_dt){
  pseudo_count <- 1
  substrate_melt <- unique(data.table(melt(substrates_dt,id.vars = 
                          c("substrate_barcode"), measure.vars = aa_cols, 
                          value.name = "aa")))
  substrate_cast <- data.table::dcast(substrate_melt, 
                                      aa ~ variable, 
                                      value.var = "aa", 
                                      fun.aggregate = length)
  substrate_cast <- substrate_cast[!aa %in% c("NA", " ", "-")]
  substrate_cast <- substrate_cast[!is.na(aa)]
  
  return(substrate_cast)
}


#' Calculate frequency of each amino acid property within PSSM
#'
#' @description This lower-level function is called within `substrate_fisher_test`. 
#' End users should use the higher-level function instead.
#' @param substrates_dt An enzymatic preference data.table.
#' @param property Character indicating which amino acid property should be 
#' analyzed. Only used if type = `aa`. See available property options 
#' by typing `aa_classificiation`.
#'
#' @return A data.table object
#' @export
#'
aa_property_freq <- function(substrates_dt, 
                             property = c("property_chemical",
                                          "property_hydropathy", 
                                          "property_volume", 
                                          "property_polar",
                                          "property_charge")){
  property <- match.arg(property)
  
  substrate_melt<- unique(data.table(melt(substrates_dt,
                                          id.vars = c("substrate_barcode"), 
                                          measure.vars = aa_cols,
                                          value.name = "aa")))
  substrate_melt<- merge(substrate_melt, aa_classification, by = "aa", all.x = T)
  substrate_cast<- data.table::dcast(substrate_melt, 
                                   eval(as.name(property)) ~ variable,
                                   fun.aggregate = length,
                                   value.var = property)
  colnames(substrate_cast)[[1]] <- property
  substrate_cast<- substrate_cast[!is.na(eval(as.name(property)))]
  return(substrate_cast)
}
