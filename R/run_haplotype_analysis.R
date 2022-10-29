#' run_aplotype_analysis
#'
#' Function to run the haplotype analysis
#' @param nloci nloci parameter
#' @param All.Pairwise nloci parameter
#' @param SID SID parameter
#' @param Tabsub Tabsub parameter
#' @param loci loci parameter
#' @param loci.ColNames loci.colNames parameter
#' @param genos genos parameter
#' @param grp grp parameter
#' @param Strict.Bin Strict.Bin parameter
#' @param Output Output parameter
#' @param Verbose Verbose parameter
#' @param Cores Cores parameter
#' @param BD.out BD.out paramter
#' @param SAFE SAFE parameter

# ======================================================================= ####
# Haplotype Analysis 'H' ________________________________________________ ####

run_haplotype_analysis <- function(nloci,
                                   All.Pairwise,
                                   SID,
                                   Tabsub,
                                   loci,
                                   loci.ColNames,
                                   genos,
                                   grp,
                                   Strict.Bin,
                                   Output,
                                   Verbose,
                                   Cores,
                                   BD.out,
                                   SAFE){

  #cat(paste(rep("_",50),collapse="","\n"))

  # Sanity check for set length and All.Pairwise=T
  if(nloci<2){
    Err.Log(Output,"Loci.No")
    stop("Analysis Stopped.", call. = F)
  }else if(All.Pairwise & nloci<=2){
    Err.Log(Output,"Loci.No.AP")
    stop("Analysis Stopped.", call. = F)
  }

  Haps.list <- H.MC.wrapper(SID,
                            Tabsub,
                            loci,
                            loci.ColNames,
                            genos,
                            grp,
                            All.Pairwise,
                            Strict.Bin,
                            Output,
                            Verbose,
                            Cores)

  if(All.Pairwise){
    if(length(BD.out[['H']])>0){
      BD.out[['H']] <- c(BD.out[['H']],Haps.list)
    }else{
      BD.out[['H']] <- Haps.list
    }
  }else{
    BD.out[['H']][[SetName]] <- Haps.list
  }

  rm(list=ls()[!(ls() %in% SAFE)])

  return(BD.out)

}
