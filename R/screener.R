multi_screener <- function(screener_raw, screener_uniprot, path, method = 
                          c("prod", "log2_sum", "w_prod"), pval_corr = FALSE, 
                          type = "aa", property = NULL, constrain = NULL){
  method <- match.arg(method)
  
  if (is.null(constrain)){
    constrain = 0.9
  }
  n_kinases <- unique(screener_raw[, kinase])
  
  message("Generating scoring matrix...\n")
  
  score_matrix <- 
    lapply(
      n_kinases,
      function(x) {
        make_scorematrix(screener_raw[kinase == x],
          screener_uniprot[uniprot_id %in% 
            screener_raw[kinase == x][, uniprot_id]],
          method,
          pval_corr = FALSE,
          type,
          property
        )
      }
    )
  
  all_pssm <- lapply(score_matrix, "[[", 1)
  all_pssm <- rbindlist(Map(cbind, all_pssm, kinase = n_kinases))

  all_score <- lapply(score_matrix, "[[", 2)
  all_score <- rbindlist(Map(cbind, all_score, kinase = n_kinases))
  
  message("Calculating cutpoints...\n")
  cp<- cutpointr::cutpointr(data = all_score, x = score, class = type, 
                            subgroup = kinase, 
                            pos_class = "sample", neg_class = "bkgrnd",
                            direction = ">=", 
                            use_midpoints = TRUE, 
                            method = maximize_boot_metric, 
                            boot_stratify = TRUE, 
                            boot_cut = 10, 
                            boot_runs = 50,
                            metric = sens_constrain, 
                            constrain_metric = specificity,
                            min_constrain = constrain)
  print(summary(cp))
  
  table_cp <- data.table(kinase=n_kinases, cutpoint = cp$optimal_cutpoint)
  
  cutpoint_plot <-
    ggplot2::ggplot(
      all_score[score > 0],
      aes(x = score, colour = type, fill = type)
    ) +
    geom_density(size = 1, alpha = 0.1) +
    scale_x_log10() +
    geom_vline(data = table_cp, aes(xintercept = cutpoint)) +
    facet_wrap(~kinase, scales = "free") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.8))
                        
  all_score <- merge(all_score,table_cp, by = "kinase")
  
  score_quantile <-
    data.table::data.table(t(rbind(
      mapply(function(x) {
        all_score[!is.infinite(score) & kinase == x][, c(max(.SD), min(.SD)),
          .SDcols = "score"
        ]
      }, n_kinases),
      mapply(function(x) {
        all_score[!is.infinite(score) & score >
          cutpoint & kinase == x][, quantile(.SD,
          probs = c(0.1, 0.5, 0.9),
          na.rm = TRUE
        ), .SDcols = "score"]
      }, n_kinases)
    )))
  
  colnames(score_quantile) <- c("max", "min", "Q10", "Q50","Q90")
  score_quantile$kinase <- n_kinases
  score_quantile <- merge(table_cp, score_quantile, by = "kinase")
  
  all_pssm <- merge(all_pssm, table_cp, by = "kinase")
  
  output_dir = file.path(path, "output")
  
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  
  ggsave("screener_cutpoints.jpg", cutpoint_plot, width = 8, height= 5.5, 
         units = c("in"), path = output_dir)
  
  output <- list(cutpoint_plot, all_pssm, score_quantile)
  return(output)
}

screener_cutpoint <- function(kinase_dt, kinase_uniprot, method = c("prod", 
                             "log2_sum", "w_prod"), pval_corr = FALSE, 
                             type = "aa", property = NULL, constrain = NULL){
  method <- match.arg(method)
  
  if (is.null(constrain)){
    constrain = 0.9
  }
  
  kinase_cd <- unique(kinase_dt[, kinase])
  score_matrix <- make_scorematrix(kinase_dt, kinase_uniprot, method, pval_corr, 
                                   type, property)
  kinase_score <- data.table(score_matrix[[2]], kinase_cd)
  kinase_fisher_melt <- score_matrix[[1]]
  
  cp<- cutpointr(kinase_score, score, type, 
                 method = maximize_boot_metric, 
                 metric = sens_constrain, 
                 constrain_metric = specificity, 
                 min_constrain = constrain, 
                 pos_class = "sample", direction = ">=", 
                 use_midpoints = TRUE, silent = TRUE,
                 boot_stratify = TRUE, boot_cut = 100)
  
  print(summary(cp))
  cutpoint <- cp$optimal_cutpoint 
  cutpoint_plot <- ggplot(kinase_score[score > 0], 
                          aes(x = score, colour = type)) +
                   geom_density(size = 1) + scale_x_log10() + theme_bw() +
                   geom_vline(xintercept = cutpoint) + 
                   xlab(paste0(kinase_cd, "\nScore"))
  
  score_quantile <- data.table(kinase = kinase_cd, 
                               quantiles=c("Q10", "Q50","Q90"), 
                               values = kinase_score[score > cutpoint][, 
                                        quantile(.SD, probs = c(0.1, 0.5, 0.9), 
                                        na.rm = T), .SDcols = "score"], 
                               threshold = cutpoint,
                               max = max(kinase_score[, score]),
                               min = min(kinase_score[, score]))
  
  score_quantile <- dcast(score_quantile, kinase + threshold  + max + min ~ 
                          quantiles, value.var = "values")
  output_data <- cbind(kinase = kinase_cd, kinase_fisher_melt, cutpoint)
  output <- list(cutpoint_plot, output_data, score_quantile)
  return(output)
}

make_scorematrix <- function(kinase_dt, kinase_uniprot, 
                             method = c("prod", "log2_sum", "w_prod"), 
                             pval_corr = FALSE, 
                             type = "aa",
                             property = NULL){
  method <- match.arg(method)
  # aa_cols <- c(paste0("-",rev(seq(1:7))),"0", paste0(seq(1:7)))
  
  kinase_cd <- unique(kinase_dt[, kinase])
  
  kinase_fisher <- substrate_fisher_test(kinase_dt, kinase_uniprot, type, property)
  valid_cols <- which(apply(kinase_fisher[[2]], 2, function(x) 
                     length(unique(x)) != 1)) #remove columns where there is no data
  kinase_fisher <- lapply(kinase_fisher, function(x) x[, valid_cols])
  kinase_fisher_melt <- fisher_long(kinase_fisher)
  
  kinase_score_dt <- 
      na.omit(unique(
        data.table::melt(kinase_dt[,.(uniprot_id, substrate_barcode), 
                                   by = aa_cols], 
                         id.vars = c("uniprot_id", "substrate_barcode"), 
                         variable.name = c("flank_pos"), 
                         value.name = "amino_acid")))
  
  kinase_score_dt[,barcode:= paste0(amino_acid, ":", flank_pos)]
  kinase_score_dt[, type:= as.factor("sample")]
  
  kinase_bkgrnd_score_dt <- 
    na.omit(unique(
      data.table::melt(kinase_uniprot[, .(uniprot_id, substrate_barcode), 
                                      by=aa_cols],
                       id.vars=c("uniprot_id", "substrate_barcode"), 
                       variable.name = c("flank_pos"), 
                       value.name = "amino_acid")))
  
  kinase_bkgrnd_score_dt[, barcode:=paste0(amino_acid, ":", flank_pos)]
  kinase_bkgrnd_score_dt[, type:= as.factor("bkgrnd")]
  
  all_score_dt <- rbind(kinase_score_dt, kinase_bkgrnd_score_dt)
  
  if (type == "aa_property"){
    all_score_dt <- merge(all_score_dt, 
                          aa_classification[, .SD, by = aa, .SDcols = property], 
                          by.x= "amino_acid", by.y = "aa")
    
    all_score_dt[, barcode:= lapply(.SD, function(x) paste0(x, ":", flank_pos)),
                 .SDcols=property]
  }
  
  all_score_dt <- 
    unique(merge(kinase_fisher_melt[, .(fisher_odds, fisher_pval), by = barcode], 
                 all_score_dt, 
                 by = "barcode")[, -c("uniprot_id")])
  
  if (isTRUE(pval_corr)){
    all_score_dt[, fisher_odds:= ifelse(fisher_pval > 0.05, 1, fisher_odds)]
    kinase_fisher_melt[, fisher_odds:= ifelse(fisher_pval > 0.05, 1, fisher_odds)]
  }

  if(method == "w_prod"){
    all_score_dt[, score:= get_score(fisher_odds, method, fisher_pval), 
                 by = substrate_barcode]
  }else{
    all_score_dt[, score:= get_score(fisher_odds, method), 
                 by = substrate_barcode]
  }
  
  unique_all_score_dt <- unique(all_score_dt[,score, by = .(
                                substrate_barcode, type)])
  
  output <- list(kinase_fisher_melt, unique_all_score_dt)
  return(output)
}

make_screener <- function(){
  # aa_cols <- c(paste0("-",rev(seq(1:7))),"0", paste0(seq(1:7)))
  path <- choose.dir()
  filenames<- toupper(list.files(path = path, full.names=T))
  screener_list <- lapply(filenames, function(x) fread(x, 
                          na.strings = c(""), header = TRUE))
  for (i in 1:length(screener_list)){
    screener_list[[i]][, file_name:= parse_sample_name(filenames[i])]
  }
  screener_dt <- rbindlist(screener_list, use.names = TRUE, fill = TRUE)
  screener_dt[, Species:= tolower(Species)]
  screener_dt[, kinase:= as.factor(sub("_.*","",file_name))]
  screener_dt[, (aa_cols):= lapply(.SD, function(x) 
                            toupper(x)), .SDcols = aa_cols]
  screener_dt[, uniprot_id:= as.factor(str_trim(uniprot_id, "both"))]
  
  id_cols <- setdiff(colnames(screener_dt), 
                     c("file_name", "AScore", "Row", "replicate"))
  valid_cols <- c("substrate_barcode", "file_name", "kinase", "uniprot_id", 
                  "class", "Species")
  match_cols <- colnames(screener_dt)[colnames(screener_dt) %in% c(valid_cols, 
                                     aa_cols)]
  
  unique_screener_dt <- unique(screener_dt[, ..match_cols])
  unique_screener_dt[, substrate_barcode:= ifelse(is.na(substrate_barcode), 
                                           gen_barcode(.SD), 
                                           substrate_barcode), .SDcols = 
                                           aa_cols]
  unique_screener_dt[, class:= as.factor("sample")]
  unique_screener_dt <- unique_screener_dt[Species == "human" | is.na(Species)]
  return(unique_screener_dt)
}

gen_barcode <- function(aa_seq){
  substrate_barcode <- apply(aa_seq, 1, paste0, collapse = "")
  substrate_barcode <- str_replace_all(substrate_barcode, "NA", "-")
  return(substrate_barcode)
}
