#!/usr/bin/env Rscript

#' bad_data_def
#'
#' Check for bad data defintions
#' @param Tab Logical indicating whether data is HLA class I/II genotyping data only.
#' @param Data.Col Input list defining which loci to use for analyses (combinations permitted).
#' @param Output Numeric Exon(s) for targeted amino acid analysis.
#' @note This function is for internal use only.
bad_data_def <- function(Tab,
                         Data.Col,
                         Output)
if(length(which(Tab[,Data.Col]==0))>0 ||
   length(which(Tab[,Data.Col]==0))>1){
  Err.Log(Output,"Bad.Data")
  stop("Analysis Stopped, Bad Data Definitions", call. = F)
}
