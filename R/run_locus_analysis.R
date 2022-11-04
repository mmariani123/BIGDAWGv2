#' run_locus_analysis
#'
#' Function to run the locus analysis
#' @param nloci nloci parameter
#' @param loci loci parameter
#' @param loci.ColNames loci.ColNames parameter
#' @param genos genos parameter
#' @param grp parameter
#' @param Strict.Bin Strict.Bin parameter
#' @param Output Output parameter
#' @param Verbose Verbose parameter
#' @param SetName SetName parameter
#' @param SAFE SAFE parameter
#' @param BD.out BD.out parameter

# ======================================================================== ####
# Locus Level 'L' ________________________________________________________ ####

run_locus_analysis <- function(nloci,
                               loci,
                               loci.ColNames,
                               genos,
                               grp,
                               Strict.Bin,
                               Output,
                               Verbose,
                               SetName,
                               SAFE,
                               BD.out){

    #cat(paste(rep("_",50),collapse=""))

    L.list <- L.wrapper(nloci,
                        loci,
                        loci.ColNames,
                        genos,
                        grp,
                        Strict.Bin,
                        Output,
                        Verbose)

    BD.out[['L']][[SetName]] <- list(binned=L.list[['AB']],
                                   freq=L.list[['AF']],
                                   OR=L.list[['OR']],
                                   chisq=L.list[['CS']],
                                   table=L.list[['FB']])

    rm(list=ls()[!(ls() %in% SAFE)])

    return(list(BD.out,SAFE))

}
