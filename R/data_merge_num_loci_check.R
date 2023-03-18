#!/usr/bin/env/Rscript

#' data_merge_num_loci_check
#'
#' Check for data merging and number of loci
#' @param Output Logical Should analysis results be written to output directory.
#' @param Merge.Output Logical Should analysis results be merged into a single file for easy access.
#' @param All.Pairwise Logical indicating whether all pairwise loci should be analyzed in haplotype analysis.
#' @param Tab The dataset read in as a table
#' @note This function is for internal use only.
data_merge_num_loci_check <- function(Output,
                                      Merge.Output,
                                      All.Pairwise,
                                      Tab)
if(Output && Merge.Output && All.Pairwise){
  if(ncol(Tab)>52){Err.Log(Output,"AllPairwise.Merge")}
}
