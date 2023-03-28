#!/usr/bin/env Rscript

#' run_dla_checks
#'
#' Check format and quality of the HLA data
#' @param Tab
#' @param Data.Col
#' @param Loci.Set
#' @param Trim
#' @param EVS.rm
#' @param Run
#' #DRB345.test
#' @param Output
#' @param Cores
#' @param Res
#' @param EPL
#' @note This function is for internal BIGDAWG use only.
run_dla_checks <- function(Tab,
                           Data.Col,
                           Loci.Set,
                           Trim,
                           EVS.rm,
                           Run,
                           #DRB345.test,
                           Output,
                           Cores,
                           Res,
                           EPL){

  if(Trim | EVS.rm | "A" %in% Run ){#| #DRB345.test){
    cat("Running HLA specific check functions...\n")
  }

  # Check Locus*Allele Formatting across all loci
  CheckCol <-
    sum(unlist(apply(Tab[,Data.Col],
                     MARGIN=c(1,2),
                     FUN = function(x){grepl("\\*",na.omit(x))}
  )
  )
  )

  TotalCol <- (dim(Tab[,Data.Col])[1] *
                 dim(Tab[,Data.Col])[2]) -
    (length(which(Tab[,Data.Col]=="^")) +
       sum(is.na(Tab[,Data.Col])))

  if(CheckCol>0 && CheckCol!=TotalCol){
    Err.Log(Output,"Bad.Format.HLA")
    stop("Analysis Stopped, HLA improperly formatted ...",
         call. = F)
  }

  # Separate DRB345 if exists as single column pair
  # and check zygosity

  #if(DRB345.test){
  #
  #  cat("Processing DRB345 column data.\n")

  #  DRBFLAG <- T

  #  # Expand DRB3/4/5 to separate column pairs
  #  Tab <- DRB345.parser(Tab)
  #  colnames(Tab) <- sapply(colnames(Tab),
  #                          FUN=gsub,
  #                          pattern="\\.1",
  #                          replacement="")

  #  # Redefine Data Columns
  #  Data.Col <- seq(3,ncol(Tab))

  #  # Define DR Loci to Process
  #  getCol <- grep("DRB",colnames(Tab))
  #  Loci.DR <- unique(colnames(Tab)[getCol])
  #
  #  # Process Loci
  #  Tab.list     <- lapply(seq_len(nrow(Tab)),
  #                         FUN=function(z){Tab[z,getCol]}
  #  )
  #
  #  Tab.tmp      <- mclapply(Tab.list,
  #                           FUN=DRB345.Check.Wrapper,
  #                           Loci.DR=Loci.DR,
  #                           mc.cores=Cores)
  #  Tab.tmp      <- do.call(rbind,Tab.tmp)
  #  Tab[,getCol] <- Tab.tmp[,grep("DRB",colnames(Tab.tmp))]
  #  Tab          <- cbind(Tab,Tab.tmp[,'DR.HapFlag'])
  #
  #  colnames(Tab)[ncol(Tab)] <- "DR.HapFlag"
  #
  #  #Identify DR345 flagged haplotypes and Write to File
  #  DR.Flags <- Tab[which(Tab[,'DR.HapFlag']!=""),
  #                  c(1,2,getCol,ncol(Tab))]
  #
  #  row.names(DR.Flags) <- NULL
  #
  #  if(Output){
  #    if(!is.null(DR.Flags)){
  #      Err.Log(Output,"Bad.DRB345.hap") ; cat("\n")
  #      write.table(DR.Flags,
  #                  file="Flagged_DRB345_Haplotypes.txt",
  #                  sep="\t",
  #                  quote=F,
  #                  row.names=F,
  #                  col.names=T)
  #    }
  #  }
  #
  #  cat("\n")
  #
  #}else{DRBFLAG <- F}

  # Separate locus and allele names if data is
  # formatted as Loci*Allele
  Tab[,Data.Col] <- apply(Tab[,Data.Col],
                          MARGIN=c(1,2),
                          FUN=Stripper)

  # Sanity Check for Resolution if Trim="T" and Trim Data
  if(Trim & CheckHLA(Tab[,Data.Col])) {
    cat("--Trimming Data.\n")
    #Tab.untrim <- Tab
    Tab[,Data.Col] <- apply(Tab[,Data.Col],
                            MARGIN=c(1,2),
                            GetField,
                            Res=Res)
    rownames(Tab) <- NULL
  }else if(Trim){
    Err.Log(Output,"Bad.Format.Trim")
    stop("Analysis Stopped, Trimming error ...",call. = F)
  }

  # Sanity Check for Expression Variant Suffix Stripping
  #if(EVS.rm & CheckHLA(Tab[,Data.Col])){
  #  cat("--Stripping Expression Variants Suffixes.\n")
  #  Tab[,Data.Col] <- apply(Tab[,Data.Col],
  #                          MARGIN=c(1,2),
  #                          gsub,
  #                          pattern="[[:alpha:]]",
  #                          replacement="")
  #  EVS.loci <- as.list(names(EPL))
  #  EPL <- lapply(EVS.loci,
  #                EVSremoval,
  #                EPList=EPL)
  #  names(EPL) <- EVS.loci
  #  rm(EVS.loci)
  #
  #}else if(EVS.rm){
  #
  #  Err.Log(Output,"Bad.Format.EVS")
  #  stop("Analysis Stopped. EVS formatting error ...",
  #       call. = F)
  #
  #}

  #EPL <- BIGDAWG::ExonPtnList

  # Sanity Check for Amino Acid Test Feasibility
  #if("A" %in% Run){

  #  cat(paste0("Running Amino Acid Analysis specific",
  #             " checks functions...\n"))

  #  Release <- EPL$Release.Version

  #  # Sanity Check for Known HLA loci in Bundled
  #  # Database Release
  #  cat(paste0("--Checking loci against database version",
  #             Release,
  #             ".\n"))
  #  test <- CheckLoci(names(EPL),
  #                    unique(colnames(Tab)[Data.Col])
  #  )
  #  if(test$Flag){
  #    Err.Log(Output,
  #            "Bad.Locus.HLA",
  #            test$Loci)
  #    stop("Analysis stopped.",
  #         call. = F)
  #  }

  #  # Sanity Check for Known HLA alleles
  #  # in Bundled Database Release
  #  #cat(paste0("--Checking alleles against database version",
  #  #           Release,
  #  #           ".\n"))
  #  #test <- CheckAlleles(EPL, Tab[,Data.Col])
  #  #if(test$Flag){
  #  #  Err.Log(Output,
  #  #          "Bad.Allele.HLA",
  #  #          test$Alleles)
  #  #  stop("Analysis stopped. Bad Allele HLA ...",call. = F)
  #  #}

  #  # Sanity Check for Analysis and HLA Allele
  #  # Resolution (MUST perform THIS STEP AFTER TRIM!!!!)
  #  #if(Res<2 | !CheckHLA(Tab[,Data.Col])){
  #  #  Err.Log(Output,"Low.Res")
  #  #  cat("You have opted to run the amino acid analysis.\n")
  #  #  stop("Analysis stopped. Low res error ...",call. = F)
  #  #}

  #} # End A if statement

  # LOCI SET COLUMN DEFINITIONS
  # This section MUST follow DRB345 processing (above)
  # on the chance that DRB345 is formatted as single column
  # and DRB3, DRB4, or DRB5 is defined in Loci.Set.

  if(missing(Loci.Set)){
    Set <- list(Data.Col)
  }else{
    Loci.Set <- lapply(Loci.Set,
                       FUN=function(x){sapply(x,toupper)})
    Set <- lapply(Loci.Set,
                  FUN=function(x){
                    seq(1,ncol(Tab))[colnames(Tab) %in% x]
                  })
  }

  # LOCUS SET DEFINED DOES NOT EXIST IN DATA
  if(!missing(Loci.Set)){
    Loci.Set <- unique(unlist(Loci.Set))
    Loci.Data <- colnames(Tab)[Data.Col]
    if(sum(Loci.Set %in% Loci.Data) != length(Loci.Set)){
      Err.Log(Output,"PhantomSets")
      stop("Analysis Stopped. PhantomSets error ...",
           call. = F) }
  }

  Release <- 'current user-defined dog/dla'
  DRBFLAG <- FALSE
  return(list(Set,Release,DRBFLAG))

}
