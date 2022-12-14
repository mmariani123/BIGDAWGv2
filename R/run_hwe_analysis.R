#' run_hwe_analysis
#'
#' Function to run the Hardy-Weinberg analysis
#' @param HLA HLA parameter
#' @param Trim Trim parameter
#' @param Tab Tab parameter
#' @param Output Output parameter
#' @param Verbose Verbose parameter
#' @param BD.out BD.out parameter

# ===================================================================================================================================== ####
# Hardy Weignberg Equilibrium 'HWE' ___________________________________________________________________________________________________ ####

run_hwe_analysis <- function(HLA,
                             Trim,
                             Tab,
                             Output,
                             Verbose,
                             BD.out){

  cat("\n>>>> STARTING HARDY-WEINBERG ANALYSIS...\n")
  #cat(paste(rep("_",50),collapse=""),"\n")
  if(HLA && Trim){
    cat("HWE performed at user defined resolution.\n")
  }else if(HLA){
    cat("HWE performed at maximum available resolution.\n")
  }

  HWE <- HWE.wrapper(Tab,Output,Verbose)
  BD.out[['HWE']] <- HWE
  rm(HWE)

  return(BD.out)

}
