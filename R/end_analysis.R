#' end_analysis
#'
#' end_analysis function
#' @param Output Logical Should analysis results be written to output directory.
#' @param OutDir OutDir parameter
#' @param BD.out BD.out parameter
#' @param Return Return parameter

end_analysis <- function(Output,
                         OutDir,
                         BD.out,
                         Return){

  # ======================================================================== ####
  cat("\n>>>>>>>>>>>>>>>>>>>>>>>>>> End Analysis <<<<<<<<<<<<<<<<<<<<<<<<<<\n")

  if(Output){
    setwd(OutDir)
    save(BD.out, file="Analysis.RData")
  }

  options(warn=0)

  if(Return){
    return(BD.out)
  }

}
