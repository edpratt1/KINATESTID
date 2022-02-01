## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  library(devtools)
#  install_local("path_to_source/KINATESTID2.tar.gz", dependencies = TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages(path_to_source, repos = NULL, type="source")

## ----results="hide", warning=FALSE, message=FALSE-----------------------------
library(KINATESTID)

## ----results="hide", warning=FALSE, message=FALSE, include=FALSE--------------
library(kableExtra)
library(magrittr)

## ----eval=TRUE----------------------------------------------------------------
path <- system.file("extdata", package = "KINATESTID")
substrates <- import_substrates(path, ptm_type = "Phosphorylation (STY)-Y",
                                legacy =  FALSE, freq = TRUE)

## ----eval=TRUE, echo=FALSE----------------------------------------------------
kbl(substrates[1:10,]) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", 
                                      "responsive", font_size = 7), 
                                       full_width = FALSE) %>%
  scroll_box(width = "100%", height = "300px")

## ----eval=FALSE---------------------------------------------------------------
#  uniq_subs <- bkgrnd_corr(substrates)

## ----eval=FALSE, message=FALSE, results="hide"--------------------------------
#  output_dir <- getwd()
#  uniprot <- import_uniprot(uniq_subs, uniprot_tryptic, path = output_dir)

## ----eval=FALSE, message=FALSE, results="hide"--------------------------------
#  pssm <- substrate_fisher_test(substrates_dt = uniq_subs, uniprot_dt = uniprot,
#                                type = "aa")

## ----eval=FALSE---------------------------------------------------------------
#  pssm <- substrate_fisher_test(substrates_dt = uniq_subs, uniprot_dt = uniprot,
#                                type = "aa_property",
#                                property = "property_chemical")

## ----eval=FALSE, warning=FALSE, message=FALSE, results='hide', error=FALSE----
#  screener <- multi_screener(screener_raw, screener_uniprot,
#                             path = output_dir,
#                             method = "prod",
#                             pval_corr = FALSE,
#                             type = "aa",
#                             norm_method = "none",
#                             constrain = 0.90)

## ---- echo=FALSE, out.width="100%"--------------------------------------------
knitr::include_graphics("screener_cutpoints.jpg")

## ----eval=FALSE, warning=FALSE, message=FALSE, error=FALSE, results='hide'----
#  candidates <- generate_substrates(pssm, uniprot, screener,
#                                    target_kinase = "BTK",
#                                    screening_kinase = "ALL",
#                                    n_hits = 10)
#  

## ----eval=FALSE, warning=FALSE, message=FALSE, results='hide'-----------------
#  venn <- peptide_intersect(substrates)[[1]]

## ---- echo = FALSE------------------------------------------------------------
knitr::include_graphics("venn_plot.jpeg")

## ----eval=FALSE, warning=FALSE, message=FALSE, results='hide'-----------------
#  heatmap <- substrates_heatmap(uniq_subs, scramble = FALSE, seed = 123)

## ---- echo = FALSE------------------------------------------------------------
knitr::include_graphics("heatmap.jpeg")

## ----eval=FALSE, warning=FALSE, message=FALSE, results='hide'-----------------
#  volcano <- pssm_volcano(pssm, odds_cutoff = 0.5, pval_cutoff = 0.05)
#  volcano

## ---- echo = FALSE------------------------------------------------------------
knitr::include_graphics("volcano_plot.jpeg")

## ----eval=FALSE---------------------------------------------------------------
#  qc <- substrates_qc(substrates, uniprot, output_dir)

## ----eval=FALSE, warning=FALSE, message=FALSE, results='hide'-----------------
#  library(motifStack)
#  motif <- new("psam", mat= pssm[[2]], name = "Affinity Logo",
#               color = colorset(alphabet = "AA", colorScheme = "chemistry"))
#  plot(motif)

## ---- echo = FALSE------------------------------------------------------------
knitr::include_graphics("motifStack_plot.jpeg")

## ----eval=FALSE---------------------------------------------------------------
#  old_peptides <- fread("path_to_source/data_file.csv")

## ----eval=FALSE---------------------------------------------------------------
#  old_peptides <- data.table(reshape2::melt(old_peptides,
#                                            id.vars = c("substrate_barcode", "kinase"),
#                                            value.name = "amino_acid",
#                                            variable.name = c("flank_pos")))
#  old_peptides[, barcode:= paste0(amino_acid, ":", flank_pos)]
#  old_peptides <- na.omit(old_peptides)

## ----eval=FALSE---------------------------------------------------------------
#  scores <- multi_candidate_screener(screener, old_peptides, "ALL", FALSE)
#  head(scores)
#  #        kinase active     perf      score   cutpoint substrate_barcode n_active
#  #     1:    ABL   TRUE   medium   4.686755   2.696951                U6       11
#  #     2:    ARG   TRUE   medium   7.352788   2.928542                U6       11
#  #     3:    AXL   TRUE     high  17.898059   2.401886                U6       11
#  #     4:    BTK   TRUE     high 162.656945   3.214976                U6       11
#  #     5:    CSK  FALSE inactive  96.024821 599.754951                U6       11

## ----eval=FALSE---------------------------------------------------------------
#  normscores <- norm_scores(scores, screener)
#  lapply(normscores, head)
#  # [[1]]
#  #   substrate_barcode ABL ARG AXL BTK CSK FES FLT3 FYN HCK JAK2 LCK LYN PYK2 SRC SYK TYRO3 YES
#  # 1                U1  75  75  83  71  51  59   96  53  54   86  46  68   68  79  70    82  83
#  # 2                U2  77  79  86  84  75  76   90  63  62   74  38  79   81  79  75    82  82
#  # 3                U3  89  91  81  86  73  65   84  69  60   58  37  81   83  84  79    78  91
#  # 4                U4  82  85  86  88  76  74   76  70  60   54  45  79   84  84  79    78  90
#  # 5                U5  88  90  87  91  62  84   97  72  64   42  44  84   72  67  70    79  74
#  # 6                U6  78  82  88  91  72  86   79  74  60   54  36  80   82  83  78    76  90
#  #
#  # [[2]]
#  # kinase norm_thresholds
#  # 1:    ABL              74
#  # 2:    ARG              77
#  # 3:    AXL              73
#  # 4:    BTK              74
#  # 5:    CSK              76
#  # 6:    FES              69

## -----------------------------------------------------------------------------
sessionInfo()

