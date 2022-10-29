#' run_amino_acid_analysis
#'
#' Function to run the amino acid analysis
#' @param UPL.flag UPL.flag parameter
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
#' @param BD.out BD.out parameter
#' @param SAFE SAFE parameter
#' @param Release Release parameter

# ======================================================================== ####
# Amino Acid Level 'A' ___________________________________________________ ####

run_amino_acid_analysis <- function(UPL.flag,
                                    nloci,
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
                                    SAFE,
                                    Release){

    #cat(paste(rep("_",50),collapse=""))

    if(UPL.flag){
      cat(paste0("Using updated protein exon alignments ",
                 "for amino acid analysis.\n"))
    }

    A.list <- A.wrapper(loci,
                        loci.ColNames,
                        genos,grp,
                        Exon,EPL,
                        Cores,
                        Strict.Bin,
                        Output,
                        Verbose)

    if(Output){
      ## write to file
      write.table(Release,
                  file = "Set_Parameters.txt",
                  sep="\t",
                  row.names = F,
                  col.names=F,
                  quote = F,
                  append=T)
    }

    BD.out[['A']][[SetName]] <- list(log=A.list[['AL']],
                                     binned=A.list[['AB']],
                                     freq=A.list[['AF']],
                                     OR=A.list[['OR']],
                                     chisq=A.list[['CS']],
                                     table=A.list[['FB']])

    rm(list=ls()[!(ls() %in% SAFE)])

}
