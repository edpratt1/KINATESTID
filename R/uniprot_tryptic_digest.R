import_uniprot <- function(substrates_dt, uniprot_digest, path){
  n_ids <- unique(unlist(substrates_dt[, uniprot_id]))
  n_ids <- n_ids[order(n_ids)]
  if (length(n_ids) <= 200){
    full_uniprot_dt <- get_uniprot(substrates_dt, uniprot_digest)
  }else {
    numrows <- length(n_ids) 
    idx <- seq(from = 0, to = numrows, by = 200)
    if (idx[length(idx)] != numrows){
      idx[length(idx) + 1] <- numrows
    }
    
    if ("file_name" %in% colnames(substrates_dt)){
      kinase = sub("_.*", "", substrates_dt[, file_name][1])
      output_dir = file.path(path, paste0(kinase, "_output")) 
    }else{
      kinase = "na"
      output_dir = file.path(path, paste0("_output")) 
    }
    
    
    if (!dir.exists(output_dir)){
      dir.create(output_dir)
    }
    
    for (i in 1:c(length(idx)-1)){
      start <- c(idx[i]+1)
      end <- idx[i+1]
      uniprot_dt <- get_uniprot(substrates_dt[mapply(function(x) any(
                               str_detect(x, paste(n_ids[start:end], collapse = 
                               "|"))), uniprot_id) == "TRUE"], uniprot_digest)
      if (is.null(uniprot_dt)){
        next
      } else{
        label <- paste(start, "to", end, sep = "_")
        file_name <- paste0(kinase, "_uniprot_", label, ".csv")
        write.csv(uniprot_dt, file.path(output_dir, file_name) , row.names=F)
      }
      
    }
    full_uniprot_dt <- multimerge(output_dir)
    write.csv(full_uniprot_dt, file.path(output_dir, paste0(kinase, 
                                        "_uniprot_merge.csv")), row.names=F)    
  }
  return(full_uniprot_dt)
}


get_uniprot<- function(substrates_dt, uniprot_digest){
  col_order <- c("uniprot_id", "peptide_no", aa_cols)
  ptm_aa <- unique(substrates_dt[, `0`]) 
  
  uniprot_key <- uniprot_digest[[1]]
  substrate_uniprot <- unique(unlist(substrates_dt[, .SD, .SDcols = "uniprot_id"]))
  uniprot_match <- uniprot_key[uniprot_id %in% substrate_uniprot]
  uniprot_nomatch <- setdiff(substrate_uniprot, uniprot_match[, uniprot_id])
  if (nrow(uniprot_match) > 0){
    uniprot_subset <- uniprot_digest[[2]][names(uniprot_digest[[2]]) %in% 
                                         uniprot_match[, "accession_number"][[1]]]
    
    uniprot_long<- data.table(reshape2::melt(uniprot_subset))
    colnames(uniprot_long) <- c("peptides", "accession_number")
    uniprot_long<- merge(uniprot_key, uniprot_long, by = "accession_number")
  }else{
    warning(paste("No match found for ", length(substrate_uniprot), "substrates"))
    uniprot_long <- NULL
  }
  
  isoform_nomatch <- sum(str_count(uniprot_nomatch, "-.*"), na.rm = T)
  protein_nomatch <- length(unique(gsub("-.*","", uniprot_nomatch)))
  percent_nomatch <- round(length(uniprot_nomatch) / (length(uniprot_match[,
                          uniprot_id]) + length(uniprot_nomatch))*100, 2)
  nomatch_msg <- 
    paste0("No match found for ", percent_nomatch, "% of uniprot ids (", 
           length(uniprot_nomatch), "/",  length(substrate_uniprot), "). ", 
           isoform_nomatch, " mismatches are protein isoforms.",  
           " Consider downloading a more complete UniprotKB file.")
  
  message(nomatch_msg)
  message("\nPTM-centering matching peptides...")
  
  if (!is.null(uniprot_long)){
    uniprot_centered <- peptide_align(uniprot_long, ptm_aa)
    uniprot_cast <- data.table(reshape2::dcast(uniprot_centered, uniprot_id + 
                              peptide_no ~ aa_cols, value.var = "flank_seq"))
    uniprot_cast <- uniprot_cast[, ..col_order]
    uniprot_cast <- add_barcode(uniprot_cast)
    uniprot_cast <- uniprot_cast[!substrate_barcode %in% substrates_dt[, 
                                substrate_barcode]]
  }else{
    warning(paste("No match found for substrates, returning NULL value"))
    return(NULL)
  }
  return(uniprot_cast)
}

generate_uniprot_digest <- function(path) {
  uniprot_tryptic_digest <- readRDS(path)
  uniprot_key <- data.table(accession_number = 
                           names(uniprot_tryptic_digest[-1]))
  accession_values <- uniprot_key[[1]]
  uniprot_id <- gsub("[;].*","",accession_values)
  uniprot_id <- as.factor(gsub(".*[|]([^.]+)[|].*","\\1",uniprot_id))
  uniprot_key <- cbind(uniprot_key, uniprot_id)
  uniprot_digest <- list(uniprot_key, uniprot_tryptic_digest)
  return(uniprot_digest)
}

multimerge = function(path){
  filenames = list.files(path = path, full.names=TRUE)
  rbindlist(lapply(filenames, function(x) data.table::fread(x, header = T)))
}