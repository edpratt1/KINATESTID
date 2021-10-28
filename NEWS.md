# KINATESTID (development version)

# KINATESTID 2.1.0

### New features
* All analysis and scoring functions only use the -4 to +4 flanking positions 
  instead of -7 to +7.

* There are now file checks for data inputs for main functions to ensure all 
  files are formatted properly prior to analysis.  

# KINATESTID 2.0.2

### New features
* `multi_screener()` now outputs a list containing relevant analysis settings
  (e.g. `method`, `pval_corr`, etc) which will be automatically used by all
  downstream functions that take screener as an input. It is no longer necessary
  to manually specify these in any function that requires screener. Example
  functions this applies to:
    + `generate_substrates()`
    + `score_candidates()`
    + `multi_candidate_screener()`

* `multi_candidate_screener` now accepts an boolean(T/F) value to screen against
  a kinase family instead of the full screener library.

### Minor improvements and bug fixes
* `multi_candidate_screener()` and `multi_screener()` now include percent progress
  bars.
  
* `multi_screener()` now provides ROC curves labeled with the optimal cutpoint
  for every kinase in the screener library.
  
* `multi_screener()` now accepts a user-defined specificity constraint using
  variable `constrain`. The default remains 90% specificity.
  
* `generate_substrates()`now allows user to manually curate amino acid residues
  when a maximum number of permutations (currently 250,000) is exceeded.

* Manual installation of Bioconductor is no longer necessary.

 
# KINATESTID 2.0.1
## Minor improvements and bug fixes
* `multi_screener()` now accepts any user-generated screener and uniprot files
  rather than only the package-provided ones.

* `substrate_fisher_test()` now works with amino acid properties.

* Vignette no longer requires the `motifStack` package to render properly.
  
