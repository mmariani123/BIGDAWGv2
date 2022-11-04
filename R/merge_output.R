#' merge_output
#'
#' merge_output function
#' @param Run Run param
#' @param BD.out BD.out param
#' @param OutDir OutDir param

merge_output <- function(Run,
                         BD.out,
                         OutDir){

  cat("\nMerging data files ...\n")
  if("HWE" %in% Run){
    Run <- Run[-which(Run=="HWE")]
  }
  if(length(Run)>=1){
    MergeData_Output(BD.out,Run,OutDir)
  }

}
