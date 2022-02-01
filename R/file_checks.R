check_pssm <- function(pssm_file){
  if(!inherits(pssm_file, "list")){
  class_error <- paste0("PSSM must be of type list, provided file is of type ",
                        toString(class(pssm_file)))
  stop(class_error, call. = TRUE)
  }
  
  if(length(pssm_file) !=2){
    length_error <- paste0("File must be two element list containing ",
                      "matrix arrays containing 1) p-values and 2) odds")
    stop(length_error, call. = TRUE)
  }
  
  if(!isTRUE("fisher_p_values" %in% names(pssm_file)) | !isTRUE("fisher_odds" %in% names(pssm_file))){
    label_error <- paste0("PSSM must contain names `fisher_p_values` and `fisher_odds`")
    stop(label_error, call.= TRUE)
  }
  return(TRUE)
}

check_screener <- function(screener_file){
  if(!inherits(screener_file, "data.table")){
    class_error <- paste0("Screener must be of type data.table, provided file is of type ",
                          toString(class(screener_file)))
    stop(class_error, call. = TRUE)
  }
  
  valid_cols <- c("substrate_barcode", "file_name", "kinase", 
                           "uniprot_id", "class", "Species", aa_cols)
  
  if (!all(valid_cols %in% colnames(screener_file))){
    missing_cols <- setdiff(valid_cols, colnames(screener_file))
    colname_error <- paste0("This file is missing required columns: ", 
                            toString(missing_cols))
    stop(colname_error, call. = TRUE)
  }
  return(TRUE)
}

check_screenerList <- function(screener_list){
  if(!inherits(screener_list, "list")){
    class_error <- paste0("Screener must be of type data.table, provided file is of type ",
                          toString(class(screener_list)))
    stop(class_error, call. = TRUE)
  }
  
  if(length(screener_list) != 4){
    length_error <- paste0("File must be four element list.")
    stop(length_error, call. = TRUE)
  }
  
  if(!inherits(screener_list[[2]], "data.table")){
    class_error <- paste0("Second element of file must be a data.table",
                          "provided file is of type  ",
                          toString(class(screener_list)))
    stop(class_error, call. = TRUE)
  }
  
  if(!inherits(screener_list[[3]], "data.table")){
    class_error <- paste0("Third element of file must be a data.table",
                          "provided file is of type  ",
                          toString(class(screener_list)))
    stop(class_error, call. = TRUE)
  }
  
  valid_cols <- c("kinase", "amino_acid", "flank_pos", "fisher_odds", 
                  "fisher_pval", "barcode", "cutpoint")
  
  if (!all(valid_cols %in% colnames(screener_list[[2]]))){
    missing_cols <- setdiff(valid_cols, colnames(screener_list[[2]]))
    colname_error <- paste0("This file is missing required columns: ", 
                            toString(missing_cols))
    stop(colname_error, call. = TRUE)
  }
  
  return(TRUE)
}

check_uniprot <- function(uniprot_file){
  if(!inherits(uniprot_file, "data.table")){
    class_error <- paste0("Uniprot File must be of type data.table, provided file is of type ",
                          toString(class(uniprot_file)))
    stop(class_error, call. = TRUE)
  }
  
  valid_cols <- c("uniprot_id", "substrate_barcode", aa_cols)
  
  if (!all(valid_cols %in% colnames(uniprot_file))){
    missing_cols <- setdiff(valid_cols, colnames(uniprot_file))
    colname_error <- paste0("This file is missing required columns: ", 
                            toString(missing_cols))
    stop(colname_error, call. = TRUE)
  }
  return(TRUE)
}

check_candidates <- function(candidates_file){
  if(!inherits(candidates_file, "data.table")){
    class_error <- paste0("Substrate Candidates File must be of type data.table, provided file is of type ",
                          toString(class(candidates_file)))
    stop(class_error, call. = TRUE)
  }
  
  valid_cols <- c("substrate_barcode", "barcode")
  
  if (!all(valid_cols %in% colnames(candidates_file))){
    missing_cols <- setdiff(valid_cols, colnames(candidates_file))
    colname_error <- paste0("This file is missing required columns: ", 
                            toString(missing_cols))
    stop(colname_error, call. = TRUE)
  }
  return(TRUE)
}
