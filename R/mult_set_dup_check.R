#!/usr/bin/env Rscript

#' mult_set_dup_check
#'
#' Multiple sets and Analysis Duplication Check
#' @param Loci.Set Numeric setting allowable missing data for running analysis (may use "ignore").
#' @param All.Pairwise Logical indicating whether all pairwise loci should be analyzed in haplotype analysis.
#' @param Run Which tests are being run
#' @param Output Logical Should analysis results be written to output directory.
#' @note This function is for internal use only.
mult_set_dup_check <- function(Loci.Set,
                               All.Pairwise,
                               Run,
                               Output){
  if(!missing(Loci.Set)){
    if(length(Loci.Set)>1 &&
      (All.Pairwise | "L" %in% Run | "A" %in% Run)){
      Err.Log(Output,"MultipleSets")
      stop("Analysis Stopped, multiple sets / analysis duplication error",
        call. = F)
    }
  }
}
