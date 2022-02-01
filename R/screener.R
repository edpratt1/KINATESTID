#' Title
#'
#' @param screener_input 
#' @param uniprot_input 
#' @param path 
#' @param method 
#' @param pval_corr 
#' @param type 
#' @param norm_method 
#' @param property 
#' @param constrain 
#'
#' @return
#' @export
#'
#' @examples
multi_screener <- function(screener_input, 
                           uniprot_input, 
                           path, 
                           method = c("prod", "log2_sum", "w_prod"), 
                           pval_corr = FALSE, 
                           type = "aa", 
                           norm_method = c("none", "bkgrnd"),
                           property = NULL, 
                           constrain = NULL){
  
  file_check <- tryCatch(
    {
      check_screener(screener_input)
      check_uniprot(uniprot_input)
    },
    error = function(e){
      message(e)
      return(FALSE)
    })
  
  if(!isTRUE(file_check)){
    return()
  }
  
  method <- match.arg(method)
  norm_method <- match.arg(norm_method)
  
  if (is.null(constrain)){
    constrain = 0.9
  }
  n_kinases <- unique(screener_input[, kinase])
  pb <- txtProgressBar(min = 0, 
                       max = length(n_kinases), 
                       initial = 0, 
                       style = 3)
  message("\nGenerating scoring matrix...\n")
  
  score_matrix <- 
    lapply(
      seq_along(n_kinases),
      function(x) {
        setTxtProgressBar(pb, x)
        output <- make_scorematrix(screener_input[kinase == n_kinases[x]],
          uniprot_input[uniprot_id %in% screener_input[kinase == n_kinases[x]][, uniprot_id]],
          method,
          pval_corr,
          type,
          norm_method,
          property
        )
        return(output)
      }
    )
  
  all_pssm <- lapply(score_matrix, "[[", 1)
  all_pssm <- rbindlist(Map(cbind, all_pssm, kinase = n_kinases))

  all_score <- lapply(score_matrix, "[[", 2)
  all_score <- rbindlist(Map(cbind, all_score, kinase = n_kinases))
  
  message("\nCalculating cutpoints...\n")
  cp<- cutpointr::cutpointr(data = all_score, 
                            x = score, 
                            class = type, 
                            subgroup = kinase, 
                            pos_class = "sample", 
                            neg_class = "bkgrnd",
                            direction = ">=", 
                            use_midpoints = TRUE, 
                            method = maximize_boot_metric, 
                            boot_stratify = TRUE, 
                            boot_cut = 30, 
                            boot_runs = 20,
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
  
  cp_point <- data.table(subgroup = n_kinases, 
                         sens = cp$sensitivity, 
                         spec = 1- cp$specificity,
                         auc = round(cp$AUC, 2))
  
  roc_plot <- plot_roc(cp, display_cutpoint = FALSE) + 
    geom_point(data = cp_point, aes(x = spec, y = sens, group = subgroup), 
               col = "black") +
    geom_text(data = cp_point, aes(x = 0.75, y = 0.1, label = auc)) +
    theme_minimal() + 
    theme(axis.text = element_text(size = 10, face = "bold"),
          axis.text.x = element_text(angle = 45),
          axis.title = element_text(size = 12, face = "bold"),
          legend.position = "none") +
    facet_wrap(~subgroup) 
  
  all_score <- merge(all_score,table_cp, by = "kinase")
  
  score_quantile <-
    data.table::data.table(t(rbind(
      mapply(function(x) {
        all_score[!is.infinite(score) & 
                    kinase == x][, c(max(.SD), min(.SD)), .SDcols = "score"]
      }, 
      n_kinases),
      mapply(function(x) {
        all_score[!is.infinite(score) & 
                    score > cutpoint & 
                    kinase == x][, quantile(.SD, probs = c(0.1, 0.5, 0.9),
                                            na.rm = TRUE
                                            ), .SDcols = "score"]
        }, 
        n_kinases)
    )))
  
  colnames(score_quantile) <- c("max", "min", "Q10", "Q50","Q90")
  score_quantile$kinase <- n_kinases
  score_quantile <- merge(table_cp, score_quantile, by = "kinase")
  
  #This needs to happen BEFORE cp selection
  if (norm_method == "bkgrnd"){
    score_quantile <- merge(score_quantile,
                            unique(all_score[, .(bkgrnd_mean, bkgrnd_sd), 
                                             by = kinase]),
                            by = "kinase")
  }
  
  all_pssm <- merge(all_pssm, table_cp, by = "kinase")
  all_pssm <- merge(all_pssm, unique(screener_input[, source, by = kinase]), 
                    by = "kinase")
  
  output_dir = file.path(path, "output")
  
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  
  ggsave("screener_cutpoints.jpg", cutpoint_plot, width = 8, height= 5.5, 
         units = c("in"), path = output_dir)

  ggsave("roc.jpg", roc_plot, width = 8, height= 5.5, 
         units = c("in"), path = output_dir)
  
  norm_tbl <- 
    data.table::data.table(reshape2::dcast(
      all_score[, sum(str_count(substrate_barcode, "Y")), 
                by = .(kinase, type)], 
      kinase ~ type, 
      value.var = "V1"))
  norm_tbl[, norm:= bkgrnd/sample]
  
  settings <- list(method = method, 
                   pval_corr = pval_corr, 
                   type = type, 
                   norm_method = norm_method, 
                   property = property)
  
  output <- list(cutpoint_plot, all_pssm, score_quantile, settings)
  return(output)
}

#' Title
#'
#' @param kinase_dt 
#' @param kinase_uniprot 
#' @param method 
#' @param pval_corr 
#' @param type 
#' @param property 
#' @param constrain 
#'
#' @return
#' @export
#'
#' @examples
screener_cutpoint <- function(kinase_dt, 
                              kinase_uniprot, 
                              method = c("prod", "log2_sum", "w_prod"), 
                              pval_corr = FALSE, 
                              type = "aa",
                              property = NULL, 
                              constrain = NULL){
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

#' Title
#'
#' @param kinase_dt 
#' @param kinase_uniprot 
#' @param method 
#' @param pval_corr 
#' @param type 
#' @param norm_method 
#' @param property 
#' @param verbose 
#'
#' @return
#' @export
#'
#' @examples
make_scorematrix <- function(kinase_dt, kinase_uniprot, 
                             method = c("prod", "log2_sum", "w_prod"), 
                             pval_corr = FALSE, 
                             type = "aa",
                             norm_method = c("none", "bkgrnd"),
                             property = NULL,
                             verbose = FALSE){
  method <- match.arg(method)
  norm_method <- match.arg(norm_method)

  kinase_cd <- unique(kinase_dt[, kinase])
  
  kinase_fisher <- substrate_fisher_test(kinase_dt, kinase_uniprot, type, property)
  valid_cols <- which(apply(kinase_fisher[[2]], 
                            2, 
                            function(x){ 
                              length(unique(x)) != 1
                              }
                            )
                      ) #remove columns where there is no data
  kinase_fisher <- lapply(kinase_fisher, function(x) x[, valid_cols])
  kinase_fisher_melt <- fisher_long(kinase_fisher)
  
  kinase_score_dt <- 
      na.omit(unique(
        data.table::melt(kinase_dt[, .(uniprot_id, substrate_barcode), 
                                   by = aa_cols], 
                         id.vars = c("uniprot_id", "substrate_barcode"), 
                         variable.name = c("flank_pos"), 
                         value.name = "amino_acid")))
  
  kinase_score_dt[, barcode:= paste0(amino_acid, ":", flank_pos)]
  kinase_score_dt[, type:= as.factor("sample")]
  
  kinase_bkgrnd_score_dt <- 
    na.omit(unique(
      data.table::melt(kinase_uniprot[, .(uniprot_id, substrate_barcode), 
                                      by=aa_cols],
                       id.vars = c("uniprot_id", "substrate_barcode"), 
                       variable.name = c("flank_pos"), 
                       value.name = "amino_acid")))
  
  kinase_bkgrnd_score_dt[, barcode:= paste0(amino_acid, ":", flank_pos)]
  kinase_bkgrnd_score_dt[, type:= as.factor("bkgrnd")]
  
  all_score_dt <- rbind(kinase_score_dt, kinase_bkgrnd_score_dt)
  all_score_dt <- all_score_dt[flank_pos %in% core_aa_cols]
  
  if (type == "aa_property"){
    all_score_dt <- merge(all_score_dt, 
                          aa_classification[, .SD, by = aa, .SDcols = property], 
                          by.x= "amino_acid", by.y = "aa")
    
    all_score_dt[, barcode:= lapply(.SD, function(x) paste0(x, ":", flank_pos)),
                 .SDcols = property]
  }
  
  all_score_dt <- 
    unique(merge(kinase_fisher_melt[, .(fisher_odds, fisher_pval), by = barcode], 
                 all_score_dt, 
                 by = "barcode")[, -c("uniprot_id")])
  
  if (isTRUE(pval_corr)){
    all_score_dt[, fisher_odds:= ifelse(fisher_pval > 0.05, 1, fisher_odds)]
    kinase_fisher_melt[, fisher_odds:= ifelse(fisher_pval > 0.05, 1, fisher_odds)]
  }

  if (method == "w_prod"){
    all_score_dt[, score:= get_score(fisher_odds, method, fisher_pval), 
                 by = substrate_barcode]
  }else{
    all_score_dt[, score:= get_score(fisher_odds, method), 
                 by = substrate_barcode]
  }
  
  if (norm_method == "bkgrnd"){
    all_score_dt <- cbind(all_score_dt, 
                          unique(all_score_dt[type == "bkgrnd"], 
                                 by = c("substrate_barcode", 
                                        "score"))[, .(
                            bkgrnd_mean = mean(score), 
                            bkgrnd_sd = sd(score))])
    
    all_score_dt[, score:= (score - bkgrnd_mean) / bkgrnd_sd]
    
    unique_all_score_dt <- unique(all_score_dt[, score, 
                                               by = .(substrate_barcode, 
                                                      type, 
                                                      bkgrnd_mean, 
                                                      bkgrnd_sd)]
                                  )
  }else{
    unique_all_score_dt <- unique(all_score_dt[,score, 
                                               by = .(substrate_barcode, 
                                                      type)]) 
    }

  if (isTRUE(verbose)){
    output <- all_score_dt
    return(output)
  }else{
    output <- list(kinase_fisher_melt, unique_all_score_dt)
    return(output)
  }
}

#' Title
#'
#' @return
#' @export
#'
#' @examples
make_screener <- function(){
  expt_file_key <- "KALIP"
  path <- choose.dir()
  filenames<- toupper(list.files(path = path, full.names=T))
  
  screener_list <- lapply(filenames, 
                          function(x) {
                            fread(x, na.strings = c(""), header = TRUE)
                          }
                         )
  
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
  unique_screener_dt[, source:= as.factor(ifelse(file_name %like% "KALIP", "EXPT", "LIT"))]
  return(unique_screener_dt)
}

gen_barcode <- function(aa_seq){
  substrate_barcode <- apply(aa_seq, 1, paste0, collapse = "")
  substrate_barcode <- str_replace_all(substrate_barcode, "NA", "-")
  return(substrate_barcode)
}
