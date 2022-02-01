generate_flank_seq <- function(central_ptm, peptide_split, peptide_length){
  #Calculate PTM position relative to C- and N-terminus flanks (assuming -7 to +7)
  #These are flipped (-7 should be n_term need to fix this)
  c_term_flank_start <- c(central_ptm - 7)
  n_term_flank_end <- c(central_ptm + 7)
  
  c_backfill <- c()
  n_backfill <- c()
  
  #If peptide truncated on the n-terminus fill missing flank positions with '-'
  if (c_term_flank_start <= 0){
    c_backfill <- rep("-",c(abs(c_term_flank_start) + 1))
    
    #Start position of relevant peptide substring
    adj_c_term_flank_start <- 1
  }else{
    adj_c_term_flank_start <- c_term_flank_start
  }
  
  #If peptide truncated on the c-terminus fill missing flank positions with '-'
  if (n_term_flank_end > peptide_length){
    n_backfill <- rep("-", abs(n_term_flank_end - peptide_length))
    
    #End position of relevant peptide substring
    adj_n_term_flank_end <- peptide_length
  }else{
    adj_n_term_flank_end <- n_term_flank_end
  }
  
  #Extract peptide substring that falls within flank positions
  centered_seq <- peptide_split[[1]][adj_c_term_flank_start:adj_n_term_flank_end]
  
  #Construct full sequence with '-' as necessary for truncations
  flank_seq <- c(c_backfill, centered_seq, n_backfill)
  return(flank_seq)
}

peptide_align <- function(file_dt, ptm_aa){
  #Generate empty data table with number of columns matching total length of flank sequence
  aa_dt <- data.table(matrix(ncol=15, nrow=0))
  
  #Column name from get_uniprot()
  peptide_colname = "peptides"
  
  centered_peptides <- list()
  
  uniprot_id <- file_dt[, uniprot_id]
  
  for (i in 1:nrow(file_dt)){
    peptide_split <- strsplit(file_dt[, ..peptide_colname][i][[1]], "")
    peptide_length <- length(peptide_split[[1]]) 
    
    #Which positions contain the modified amino acid of interest
    central_ptm <- which(peptide_split[[1]] == ptm_aa) #THIS NEEDS TO BE FLEXIBLE
    
    if (length(central_ptm) == 0){
      next
    }else{
        for (j in 1:length(central_ptm)){
        flank_seq <- generate_flank_seq(central_ptm[j], 
                                        peptide_split, 
                                        peptide_length)
        
        peptide_num <- paste0(i, "_", j)
        
        seq_dt <- data.table(uniprot_id = rep(uniprot_id[i], length(flank_seq)), 
                             peptide_no = rep(peptide_num,length(flank_seq)), 
                             aa_cols, 
                             flank_seq)
        
        #This could blow up quickly if the file is large, need an alternative way to determine position
        centered_peptides[[length(centered_peptides) + 1]] <- seq_dt
      }
    }
  }
  centered_dt<- do.call(rbind, centered_peptides)
}
