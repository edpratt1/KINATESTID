#' Import and annotate enzyme substrate preference files
#' 
#' @description This function imports enzymatic preference data sheets and 
#' creates an annotated preference data.table through a combination of parsing 
#' file names and column data. Check the vignette for the required formatting of
#' input data files.
#' 
#' @param path The directory containing the analyzed phosphoproteomics data.
#' 
#' @param ptm_type The post-translational modification to be analyzed and its 
#' targeted the amino acid. Currently, only one ptm type can be analyzed at a 
#' time. A list of ptm options can be accessed by typing `ptm_key`.
#' 
#' @param legacy A logical parameter where TRUE is for files generated using the
#' Galaxy-P KinaMine workflow. Is set to FALSE by default.
#' 
#' @param freq A logical parameter where TRUE indicates substrate frequency 
#' across technical replicates should be recorded. Is set to TRUE by default.
#' 
#' @param ref_col Only applicable to legacy files. A string with the name of the
#' column containing the uniprot identifier information. If none is specified, 
#' it will search for a Reference column automatically.
#'
#' @return A data.table object.
#' @export
#'
#' @examples
#'     path <- system.file("extdata", package = "KINATESTID")
#'     substrates <- import_substrates(path, ptm_type = "Phosphorylation (STY)-Y",
#'                                     legacy =  FALSE, freq = TRUE)
#'                                     
import_substrates <- function(path, 
                              ptm_type = 
                                c("Oxidation (M)-M", "Carbamidomethylation-C", 
                                  "Deamidation (NQ)-N", "Deamidation (NQ)-Q", 
                                  "Phosphorylation (STY)-Y", 
                                  "Phosphorylation (STY)-T", 
                                  "Phosphorylation (STY)-S",
                                  "Acetylation (Protein N-term)-E", 
                                  "Acetylation (Protein N-term)-D",
                                  "Acetylation (Protein N-term)-M", 
                                  "Acetylation (Protein N-term)-S",
                                  "Acetylation (Protein N-term)-A", 
                                  "Acetylation (Protein N-term)-T",
                                  "Acetylation (Protein N-term)-G",
                                  "Acetylation (Protein N-term)-C",
                                  "Acetylation (Protein N-term)-V", 
                                  "Pyro-glu from Q-Q"), 
                              legacy = FALSE, 
                              freq = FALSE, 
                              ref_col = NULL){
  
  #Check if directory exists
  if(!dir.exists(path)){
    stop("Directory provided does not exist")
  }
  
  ptm_type <- match.arg(ptm_type)
  
  #Import files and generate UID using flank sequence 
  substrates_dt <- multi_sample(path, ptm_type, legacy)
  substrates_dt <- add_barcode(substrates_dt, freq)
  
  #Parse reference column to extract Uniprot ID if column does not already exist
  if(isTRUE(legacy) && !"uniprot_id" %in% colnames(substrates_dt)){
    if(is.null(ref_col)){
      warning("No reference column provided, using `Reference` automatically...")
      ref_col <- "Reference"
      substrates_dt<- add_uniprotid(substrates_dt, ref_col)
    }else{
      substrates_dt<- add_uniprotid(substrates_dt, ref_col)
    }
  }
  return(substrates_dt)
}

#' Import enzyme substrate preference files and parse sample type
#' 
#' @description  This lower-level function is called within `import_substrates`. 
#' End users should use the higher-level function instead.
#' 
#' @param path The directory containing the analyzed phosphoproteomics data.
#' 
#' @param ptm_type The post-translational modification to be analyzed and its 
#' targeted the amino acid. Currently, only one ptm type can be analyzed at a 
#' time. A list of ptm options can be accessed by typing `ptm_key`.
#' 
#' @param legacy A logical parameter where TRUE is for files generated using the
#' Galaxy-P KinaMine workflow. Is set to FALSE by default.
#'
#' @return A data.table object
#' @export
#'
multi_sample <- function(path, ptm_type = c("Oxidation (M)-M", "Carbamidomethylation-C",
                        "Deamidation (NQ)-N", "Deamidation (NQ)-Q", 
                        "Phosphorylation (STY)-Y", "Phosphorylation (STY)-T", 
                        "Phosphorylation (STY)-S", "Acetylation (Protein N-term)-E", 
                        "Acetylation (Protein N-term)-D", "Acetylation (Protein N-term)-M", 
                        "Acetylation (Protein N-term)-S", "Acetylation (Protein N-term)-A", 
                        "Acetylation (Protein N-term)-T", "Acetylation (Protein N-term)-G", 
                        "Acetylation (Protein N-term)-C", "Acetylation (Protein N-term)-V", 
                        "Pyro-glu from Q-Q"), legacy = FALSE){
  ptm_type <- match.arg(ptm_type)
  
  #File name keywords used to classify file & sample type
  samp_key <- "PLUS"
  cntrl_key <- "MINUS"
  substrate_key <- "SUBSTRATE"
  freq_key <- c("SBF", "FREQ")
  
  #Quality score cutoff for peptide-spectra
  ascore_cutoff <- 30
  
  filenames<- toupper(list.files(path = path, full.names=T))
  
  #Legacy files are those generated using the Galaxy-P/Kinamine workflow
  if (isTRUE(legacy)){
    id_cols <- c("uniprot_id", "file_name", "Substrates", "Species", "Reference",
                 "Phosphite", "Phosphosite")
    
    valid_cols <- paste0(c(id_cols, aa_cols))
    
    substrate_filenames <- filenames[filenames %like% substrate_key]
    
    #Pulls first 19 columns of csv, if flank sequence length changes this also needs to be changed
    substrates_list <- lapply(substrate_filenames, 
                              function(x) data.table::fread(x, 
                                                            na.strings = c(""), 
                                                            select=c(1:19), 
                                                            header = T))
    for (i in 1:length(substrates_list)){
      substrates_list[[i]][, file_name:= parse_sample_name(substrate_filenames[i])]
    }
    substrates_dt <- rbindlist(substrates_list)
    
    #Subset data input to only necessary columns (defined in valid_cols)
    subset_cols <- colnames(substrates_dt)[colnames(substrates_dt) %in% valid_cols]
    subset_dt <- substrates_dt[, ..subset_cols]
    
    subset_dt[, class:= ifelse(file_name %like% samp_key, "sample", 
                        ifelse(file_name %like% cntrl_key, "control", 
                        "unclassified"))]
    
    #Parse technical replicate number from file name
    subset_dt[, replicate:=  sub(".+?_", "", file_name)]
    subset_dt[, replicate:= as.factor(gsub("[^0-9]", "", replicate))]
  }else{
    substrate_filenames <- filenames 
    substrates_list <- data.table::fread(substrate_filenames, na.strings= c(""), 
                                         header = T)
    
    #Subset input based on PTM type defined by user
    subset_dt <- substrates_list[`Modification Name` == 
                                   ptm_key[ptm_id == ptm_type][[1]] &
                                   `Amino Acid Letter` == 
                                   ptm_key[ptm_id == ptm_type][[2]]]
    
    
    subset_dt[, class:= ifelse(`File Name` %like% samp_key, "sample", 
                               ifelse(`File Name` %like% cntrl_key, 
                                      "control", "unclassified"))]
    
    #Parse technical replicate number from file name
    subset_dt[, replicate:= str_extract(`File Name`, "_R[0-9]")]
    subset_dt[, replicate:= as.factor(gsub("[^0-9]", "", replicate))]
    
    #Subset based on quality score
    subset_dt <- subset_dt[AScore > ascore_cutoff]
    
    #Coerce to valid column name
    setnames(subset_dt, "UniProt Identifier", "uniprot_id")
    setnames(subset_dt, "File Name", "file_name")
  }
  
  #Replace any empty flank positions with '-'
  subset_dt[, (aa_cols):= lapply(.SD, function(x) ifelse(is.na(x), "-", 
                                paste(x))), .SDcols= aa_cols]
  return(subset_dt)
}


parse_sample_name <- function(sample_paths){
  sample_subpaths <- read.table(text = sample_paths,
                                sep = "/", 
                                colClasses="character")
  
  sample_names <- sample_subpaths[, ncol(sample_subpaths)]
}


#' Generate unique barcode identifier for peptide substrates
#' 
#' @description This lower-level function is called within `import_substrates`. 
#' End users should use the higher-level function instead.
#' 
#' @param substrates_dt An enzymatic preference data.table. 
#' 
#' @param freq a logical parameter where TRUE indicates substrate frequency 
#' across technical replicates should be recorded.
#' 
#' @param legacy a logical parameter where TRUE is for files generated using 
#' the old Galaxy-P KinaMine workflow. Is set to FALSE by default.
#'
#' @return A data.table object
#' @export
#'
add_barcode <- function(substrates_dt, freq = FALSE, legacy = FALSE){
  #Isolate flank sequence (defined by aa_cols)
  substrate_aa <- substrates_dt[,.SD, .SDcols = aa_cols]
  
  substrate_barcode <- apply(substrate_aa, 1, paste0, collapse = "")
  substrates_dt<- unique(cbind(substrates_dt, substrate_barcode))
  
  #Calculate frequency of each peptide UID
  if (isTRUE(freq)){
    if (isTRUE(legacy)){
      unique_substrate_dt <- unique(substrates_dt[, -c("Reference")])
      barcode_freq <- data.table(reshape2::dcast(unique_substrate_dt, 
                                                 substrate_barcode ~ class,
                                                 fun.aggregate = length, 
                                                 value.var = "substrate_barcode")) 
    }else{
      unique_substrate_dt <- unique(substrates_dt[, substrate_barcode, 
                                                  by = .(replicate, class)])
      barcode_freq <- 
        data.table(reshape2::dcast(unique_substrate_dt, 
                                   substrate_barcode ~ class, 
                                   fun.aggregate = length, 
                                   value.var = "substrate_barcode"))
      
      if (!"control" %in% colnames(barcode_freq)){
        barcode_freq[, control:= 0]
      }
    }
    substrates_dt <- merge(substrates_dt, barcode_freq, by = "substrate_barcode")    
  }
  return(substrates_dt)
}

#' Parse UniprotKB ID from phosphoproteomics data files
#' 
#' @description This lower-level function is called within `import_substrates`. 
#' End users should use the higher-level function instead.
#' 
#' @param file_dt A data.table.
#' 
#' @param col_name The column name containing UniprotKB IDs.
#'
#' @return A data.table object
#' @export
#'
add_uniprotid <- function(file_dt, col_name){
  file_expand <- data.table(do.call(rbind, 
                                    apply(file_dt, 1, function(x) {
                                      do.call(expand.grid, strsplit(x, ";"))
                                      })))
  
  if (all(grepl("|", file_expand[, ..col_name]))){
    #Expand data.table to create entry for every UniprotKB ID in reference column
    file_expand[, uniprot_id:= str_extract_all(unlist(.SD), "\\|(.*?)\\|"), .SDcols = col_name]
    
    dropped <- nrow(file_expand[uniprot_id == "character(0)"])
    warning(paste0("No Uniprot ID found for ", dropped, " substrates."))
    
    #Extract uniprot_id substring
    file_expand[, uniprot_id:= ifelse(uniprot_id == "character(0)", NA, 
                                      str_sub(unlist(uniprot_id), 
                                              start = 2, 
                                              end = str_length(unlist(uniprot_id)) - 1))]
    
    #Trim white space
    file_expand[, uniprot_id:= str_trim(uniprot_id, "both")]
  }else{
    file_expand[, uniprot_id:= eval(parse(text = paste(col_name)))]
  }
  
  #Coerce replicate, control, and sample columns to numeric
  cols <- c("replicate", "control", "sample")
  if(all(cols %in% colnames(file_dt))){
    file_expand[, (cols):= lapply(.SD, function(x) as.numeric(
                                  as.character(x))), .SDcols = cols]
  }
  return(file_expand)
}

#' Subset common peptide sequences
#'
#' @description This function will automatically select peptide sequences which
#' are present in maximum number of treated replicates and 
#' minimum number of control samples.
#' 
#' @param substrates_dt An enzymatic preference data.table imported 
#' using `import_substrates()`
#'
#' @return A data.table object.
#' @export
#' @examples
#' uniq_subs <- bkgrnd_corr(substrates)
bkgrnd_corr <- function(substrates_dt){
  cntrl_limit <- min(substrates_dt[,"control"])
  sample_limit <- max(substrates_dt[,"sample"])
  corr_dt <- substrates_dt[control == cntrl_limit & sample == sample_limit]
  return(corr_dt)
}