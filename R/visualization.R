substrates_qc <- function(substrates_dt, path){
  unique_substrates <- bkgrnd_corr(substrates_dt)
  bin_dt <- unique(substrates_dt[,substrate_barcode, by = .(sample, control)])
  
  sample_bin_graph <- ggplot(bin_dt, aes(x = sample)) + 
    geom_histogram(stat = "count") + theme_bw() + 
    xlab("No of Samples Present") +
    ylab("No of Unique Substrates")
  
  venn_data <- peptide_intersect(substrates_dt)
  
  sub_heatmap <- substrates_freq(unique_substrates, T)
  
  counts_dt <- unique(substrates_dt[,substrate_barcode,by=.(file_name, class)])
  counts_graph <- ggplot(counts_dt, aes(x = file_name, fill = class)) + 
                  geom_histogram(stat = "count") + theme_bw() + 
                  facet_wrap(~class, scales = "free") +
                  theme(axis.text = element_text(angle = 45, hjust = 1), 
                  plot.margin = unit(c(0,1,0,2), "cm"), legend.position =
                  "none") + xlab("")
  
  top_row <- plot_grid(counts_graph)
  
  first_col <- plot_grid(sample_bin_graph, venn_data[[1]], ncol=1)
  second_col <- plot_grid(sub_heatmap)
  bottom_row <- plot_grid(first_col, second_col, ncol = 2, rel_widths = 
                         c(1, 1.5))
  
  
  gg_all <- plot_grid(top_row, bottom_row, ncol = 1, rel_heights = c(1,1.5))
  
  output_dir = file.path(path,"output")
  
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  
  ggsave("substrates_qc_plot.jpg", gg_all, width = 8.5, height= 11, units = 
        c("in"), path = output_dir)
  
  qc_output <- list(gg_all, venn_data[[2]])
  return(qc_output)
}

peptide_intersect <- function(substrates_dt){
  intersect_var = "replicate"
  sample_key <- unique(substrates_dt[, ..intersect_var])[[1]]
  no_samples <- length(sample_key)
  sample_perm <- combn(sample_key, 2, simplify=F)
  
  comparisons <- list()
  barcode <-  unique(substrates_dt[, substrate_barcode, by = intersect_var])
  barcode_lengths <- barcode[,length(unique(substrate_barcode)), 
                             by= intersect_var]
  
  for (i in 1:no_samples){
    comparisons[[i]] <-  unique(
      barcode[replicate == sample_key[i]][,substrate_barcode]
    )
    names(comparisons)[i] <- sample_key[i]
  }
  
  peptide_overlaps <- list()
  for (j in 1:length(sample_perm)){
    index<-paste(sample_perm[[j]])
    peptide_overlaps[[j]] <- Reduce(intersect, 
                                    comparisons[names(comparisons) %in% dput(index)])
    names(peptide_overlaps)[j] <- paste0(index, collapse = ":")
  }
  
  if (no_samples >=3){
    peptide_overlaps[[length(peptide_overlaps) + 1]] <- Reduce(intersect, 
                                                               peptide_overlaps)
    names(peptide_overlaps)[length(peptide_overlaps)] <- paste0(sample_key, 
                                                                collapse = ":")
  }
  
  if (no_samples == 1){
    venn_plot <- ggplot() + theme_void()
  }
  
  if (no_samples == 2){
    col_pal <- brewer.pal(no_samples, "Dark2")
    col_pal <- col_pal[1:no_samples] #workaround for required minimum colors
    area1 <- barcode_lengths[1][[2]]
    area2 <- barcode_lengths[2][[2]]
    n12 <- length(peptide_overlaps[[1]])
    labels <- paste("Replicate",sample_key, sep="\t ")
    venn_plot<- draw.pairwise.venn(area1, area2, n12, 
                                  scaled = T, fill = col_pal, col = 
                                  rep("white", 2), category = labels, 
                                  fontface = rep("bold", 3), cat.fontfamily = 
                                  rep("Arial", 2), cat.fontface = 
                                  rep("bold", 2), cat.col = col_pal, cex = 
                                  rep(1.5, 3), margin=0.08)
    venn_plot <- as.ggplot(grobTree(venn_plot))
  }
  
  if (no_samples == 3){
    col_pal <- brewer.pal(3, "Dark2")
    area1 <- barcode_lengths[1][[2]]
    area2 <- barcode_lengths[2][[2]]
    area3 <- barcode_lengths[3][[2]]
    n12 <- length(peptide_overlaps[[1]])
    n13 <- length(peptide_overlaps[[2]])
    n23 <- length(peptide_overlaps[[3]])
    n123 <- length(peptide_overlaps[[4]])
    labels <- paste("Replicate",sample_key, sep="\t ")
    venn_plot<- draw.triple.venn(area1, area2, area3, n12, n23, n13, n123, 
                                scaled = T, fill = col_pal, col =
                                rep("white", 3), category = labels, fontface =
                                rep("bold", 7), cat.fontfamily = 
                                rep("Arial", 3), cat.fontface = rep("bold", 3),
                                cat.col = col_pal, cex = rep(1.5, 7), 
                                margin=0.08)
    venn_plot <- as.ggplot(grobTree(venn_plot))
  }
  overlap_output <- list(venn_plot, peptide_overlaps)
  return(overlap_output)
  
}

pssm_volcano <- function(pssm, odds_cutoff = NULL, pval_cutoff = NULL){
  if (is.null(odds_cutoff)){
    odds_cutoff <- 0.5
  }
  if(is.null(pval_cutoff)){
    pval_cutoff <- 0.05
  }
  
  pssm_long <- fisher_long(pssm)
  pssm_long[, log2Odds:= log2(fisher_odds)]
  
  volcano_plot <- EnhancedVolcano(pssm_long, lab = pssm_long[, barcode], 
                                  x='log2Odds', y='fisher_pval', 
                                  FCcutoff = odds_cutoff, 
                                  pCutoff = pval_cutoff, 
                                  boxedLabels = TRUE, 
                                  drawConnectors = TRUE, 
                                  title = NULL, 
                                  subtitle = NULL,
                                  xlab = bquote(~Log[2]~"fisher odds"))
  
}

substrates_heatmap <- function(substrates_dt, scramble = FALSE, seed = NULL){
  if (is.null(seed)){
    seed = 42
  }
  
  set.seed(seed)
  
  aa_prop_table <- aa_freq(substrates_dt)
  
  if (isTRUE(scramble)){
    ptm_col <- aa_prop_table[,`0`]
    order<-sample(1:nrow(aa_prop_table))
    aa_prop_table <- aa_prop_table[order, ]
    aa_prop_table$`0` <- ptm_col
    aa_clean <- Filter(function(x)!all(is.na(x)), aa_prop_table)
    aa_heatmap <- data.matrix(aa_clean[,-1])
    aa_heatmap_plot <- pheatmap(aa_heatmap, cluster_cols = F, cluster_rows = F, 
                                scale= "column")
    aa_heatmap_plot <- as.ggplot(aa_heatmap_plot)
  }else{
    aa_clean <- Filter(function(x)!all(is.na(x)), aa_prop_table)
    aa_heatmap <- data.matrix(aa_clean[,-1])
    rownames(aa_heatmap)<- aa_clean[[1]]
    aa_heatmap_plot <- pheatmap(aa_heatmap, cluster_cols = F, cluster_rows = F, 
                                scale= "column")
    aa_heatmap_plot <- as.ggplot(aa_heatmap_plot)
  }
  return(aa_heatmap_plot)  
}

silico_heatmap <- function(output_dt){
  dm_palette <- colorRampPalette(c("#FFFFFF", "#FFB7B4", "#FF6361", "#FF3127"))
  key <- data.table(perf = c("inactive", "low", "medium", "high"), color = c(0, 1, 2, 3))
  
  output_dt$perf <- factor(output_dt$perf, levels = c("inactive", "low", "medium", "high"))
  output_merge <- merge(output_dt, key, by = "perf")
  dm_rownames <- dcast(output_merge, substrate_barcode ~ kinase, value.var = "color")[,1]
  dm <- data.matrix(dcast(output_merge, substrate_barcode ~ kinase, value.var = "color", value.factor=TRUE)[,-1])
  rownames(dm) <- dm_rownames
  dm_plot <- pheatmap(t(dm), cluster_rows = FALSE, cluster_cols = FALSE, color = dm_palette(5))
  dm_plot <- as.ggplot(dm_plot)
}


plot_activity<- function(output_dt, kinase_set){
  output_dt$perf <- as.factor(output_dt$perf)
  output_dt$perf <- factor(output_dt$perf, levels = c("inactive", "low", 
                          "medium", "high"))
  output_dt$kinase <- droplevels(output_dt$kinase)
  output_dt$kinase <- factor(output_dt$kinase, levels = 
                            kinase_set[kinase_set %in% output_dt[,kinase]])
  output_plot <- ggplot(output_dt, aes (x = kinase, y = perf, color = perf)) + 
                 geom_point(size = 5)+ theme_bw() + 
                 facet_wrap(~substrate_barcode) + theme(axis.text = 
                 element_text(face = "bold", size = "12"), axis.title =
                 element_text(face = "bold", size = 12), legend.position = "none",
                 axis.text.x = element_text(angle = 45, vjust = 0.6), 
                 strip.text = element_text(face="bold", size = 12)) 
}

substrates_property <- function(substrates_dt, uniprot_dt, path){
  property = c("property_chemical", "property_hydropathy", "property_volume", 
              "property_polar")
  
  property_plots <- list()
  property_tables <- list()
  
  for (i in 1:length(property)){
    fisher_tables <- substrate_fisher_test(substrates_dt, uniprot_dt, 
                                           "aa_property", property[i])
    fisher_melt <- melt(fisher_tables[[2]], varnames = 
                       c("property","flank_pos"), value.name = "fisher_odds")
    fisher_plot <- ggplot(fisher_melt, aes(x = flank_pos, y = fisher_odds, 
                         colour = property)) + geom_point(size = 2) + 
                         facet_wrap(~property, scales = "free_x", ncol = 2) + 
                         geom_line(size = 1) + geom_hline(yintercept = 1, 
                         linetype = "dashed") + theme_bw() + xlab("") + ylab(
                         "Fisher Odds") + theme(legend.position = "none")
    property_plots[[i]] <- fisher_plot
    property_tables[[i]] <- fisher_tables
  }
  
  top_row <- plot_grid(property_plots[[1]], property_plots[[3]])
  bottom_row <- plot_grid(property_plots[[2]], property_plots[[4]])
  substrates_property_plot <- plot_grid(top_row, bottom_row, rel_heights = 
                                          c(1, 0.5),ncol=1)
  
  output_dir = file.path(path, "output")
  
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  
  ggsave("substrates_properties_plot.jpg", substrates_property_plot, width = 11,
         height= 8.5, units = c("in"), path = output_dir)
  
  property_pval_combn <- rbind(property_tables[[1]][[1]], 
                               property_tables[[2]][[1]], 
                               property_tables[[3]][[1]], 
                               property_tables[[4]][[1]])
  
  property_odds_combn <- rbind(property_tables[[1]][[2]], 
                               property_tables[[2]][[2]], 
                               property_tables[[3]][[2]], 
                               property_tables[[4]][[2]])
  
  property_tables_combn <- list(property_pval_combn, property_odds_combn)
  return(property_tables_combn)
}
