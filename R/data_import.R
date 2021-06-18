import_substrates <- function(path, ptm_type = c("Oxidation (M)-M", 
                              "Carbamidomethylation-C", "Deamidation (NQ)-N", 
                              "Deamidation (NQ)-Q", "Phosphorylation (STY)-Y", 
                              "Phosphorylation (STY)-T", "Phosphorylation (STY)-S",
                              "Acetylation (Protein N-term)-E", "Acetylation (Protein N-term)-D",
                              "Acetylation (Protein N-term)-M", "Acetylation (Protein N-term)-S",
                              "Acetylation (Protein N-term)-A", "Acetylation (Protein N-term)-T",
                              "Acetylation (Protein N-term)-G", "Acetylation (Protein N-term)-C",
                              "Acetylation (Protein N-term)-V", "Pyro-glu from Q-Q"), 
                              legacy = FALSE, freq = FALSE, ref_col = NULL){
  
  ptm_type <- match.arg(ptm_type)
  substrates_dt <- multi_sample(path, ptm_type, legacy)
  substrates_dt <- add_barcode(substrates_dt, freq)
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
  samp_key <- "PLUS"
  cntrl_key <- "MINUS"
  substrate_key <- "SUBSTRATE"
  freq_key <- c("SBF", "FREQ")
  # aa_cols <- c(paste0("-",rev(seq(1:7))),"0", paste0(seq(1:7)))
  ascore_cutoff <- 30
  
  # path <- choose.dir()
  filenames<- toupper(list.files(path = path, full.names=T))
  
  if (isTRUE(legacy)){
    id_cols <- c("uniprot_id", "file_name", "Substrates", "Species", "Reference",
                 "Phosphite", "Phosphosite")
    
    valid_cols <- paste0(c(id_cols, aa_cols))
    
    substrate_filenames <- filenames[filenames %like% substrate_key]
    substrates_list <- lapply(substrate_filenames, function(x) data.table::fread(x, na.strings
                              = c(""), select=c(1:19), header = T))
    for (i in 1:length(substrates_list)){
      substrates_list[[i]][, file_name:= parse_sample_name(substrate_filenames[i])]
    }
    substrates_dt <- rbindlist(substrates_list)
    subset_cols <- colnames(substrates_dt)[colnames(substrates_dt) %in% valid_cols]
    subset_dt <- substrates_dt[, ..subset_cols]
    
    subset_dt[, class:= ifelse(file_name %like% samp_key, "sample", 
                        ifelse(file_name %like% cntrl_key, "control", 
                        "unclassified"))]
    subset_dt[, replicate:=  sub(".+?_", "", file_name)]
    subset_dt[, replicate:= as.factor(gsub("[^0-9]", "", replicate))]
  }else{
    substrate_filenames<- filenames 
    substrates_list <- data.table::fread(substrate_filenames, na.strings= c(""), header = T)
    
    subset_dt <- substrates_list[`Modification Name` == ptm_key[ptm_id == 
                                ptm_type][[1]] & `Amino Acid Letter` == 
                                ptm_key[ptm_id == ptm_type][[2]]]
    
    
    subset_dt[, class:= ifelse(`File Name` %like% samp_key, "sample", 
                               ifelse(`File Name` %like% cntrl_key, "control", "unclassified"))]
    subset_dt[, replicate:= str_extract(`File Name`, "_R[0-9]")]
    subset_dt[, replicate:= as.factor(gsub("[^0-9]", "", replicate))]
    subset_dt <- subset_dt[AScore > ascore_cutoff]
    setnames(subset_dt, "UniProt Identifier", "uniprot_id")
    setnames(subset_dt, "File Name", "file_name")
  }
  subset_dt[, (aa_cols):= lapply(.SD, function(x) ifelse(is.na(x), "-", 
                                paste(x))), .SDcols= aa_cols]
  return(subset_dt)
}

parse_sample_name <- function(sample_paths){
  sample_subpaths <- read.table(
    text=sample_paths,sep="/", colClasses="character"
  )
  sample_names <- sample_subpaths[,ncol(sample_subpaths)]
}

add_barcode <- function(substrates_dt, freq = FALSE, legacy = FALSE){
  # aa_cols <- c(paste0("-",rev(seq(1:7))),"0", paste0(seq(1:7)))
  
  substrate_aa<- substrates_dt[,.SD, .SDcols = aa_cols]
  substrate_barcode <- apply(substrate_aa, 1, paste0, collapse = "")
  substrate_sampleid <- paste0(substrates_dt[,1], substrate_barcode)
  substrates_dt<- unique(cbind(substrates_dt, substrate_barcode))
  
  if (isTRUE(freq)){
    if (isTRUE(legacy)){
      unique_substrate_dt <- unique(substrates_dt[, -c("Reference")])
      barcode_freq <- dcast(unique_substrate_dt, substrate_barcode~class, fun.aggregate =
                              length, value.var = "substrate_barcode") 
    }else{
      unique_substrate_dt <- unique(substrates_dt[, substrate_barcode, by = .(
                                   replicate, class)])
      barcode_freq <- dcast(unique_substrate_dt, substrate_barcode ~ class, 
                            fun.aggregate = length, value.var = "substrate_barcode")
      if (!"control" %in% colnames(barcode_freq)){
        barcode_freq[, control:= 0]
      }
    }
    
    substrates_dt <- merge(substrates_dt, barcode_freq, by = "substrate_barcode")    
  }
  return(substrates_dt)
}

add_uniprotid <- function(file_dt, col_name){
  file_expand <- data.table(do.call(rbind, apply(file_dt, 1, function(x) { 
    do.call(expand.grid, strsplit(x, ";"))})))
  
  if (all(grepl("|", file_expand[, ..col_name]))){
    file_expand[, uniprot_id:= str_extract_all(unlist(.SD),"\\|(.*?)\\|"), .SDcols = col_name]
    
    dropped <- nrow(file_expand[uniprot_id == "character(0)"])
    warning(paste0("No Uniprot ID found for ", dropped, " substrates."))
    
    file_expand[, uniprot_id:= ifelse(uniprot_id == "character(0)", NA, str_sub(
               unlist(uniprot_id), start = 2, end = str_length(unlist(
               uniprot_id)) - 1))]
    file_expand[, uniprot_id:= str_trim(uniprot_id, "both")]
  }else{
    file_expand[, uniprot_id:= eval(parse(text=paste(col_name)))]
  }
  
  cols <- c("replicate", "control", "sample")
  if(all(cols %in% colnames(file_dt))){
    file_expand[, (cols):= lapply(.SD, function(x) as.numeric(
                                  as.character(x))), .SDcols = cols]
  }
  return(file_expand)
}

bkgrnd_corr <- function(substrates_dt){
  cntrl_limit <- min(substrates_dt[,"control"])
  sample_limit <- max(substrates_dt[,"sample"])
  corr_dt <- substrates_dt[control == cntrl_limit & sample == sample_limit]
  return(corr_dt)
}