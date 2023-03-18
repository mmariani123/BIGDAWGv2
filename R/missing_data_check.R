#!/usr/bin/env Rscript

#' missing_data_check
#'
#' Check for missing dataV
#' @param Missing Numeric setting allowable missing data for running analysis (may use "ignore").
#' @param Output Logical Should analysis results be written to output directory.
#' @param Tab The dataset read in as a table
#' @param Data.Col The data column names
#' @param NAStrings Vector of various character strings defining NA values
#' @note This function is for internal use only.
missing_data_check <- function(Missing,
                               Output,
                               Tab,
                               Data.Col,
                               NAstrings){
  if(Missing == "ignore"){
    cat("Ignoring any missing data...\n")
    Err.Log(Output,"Ignore.Missing")
    rows.rm <- NULL
  }else{
    if(Missing > 2){
      if("H" %in% Run){Err.Log(Output,"Big.Missing")}
    }
    cat(paste0("Removing any missing data.",
              " This will affect Hardy-Weinberg Equilibrium test.\n"))
    geno.desc <- summaryGeno.2(Tab[,Data.Col], miss.val=NAstrings)
    test <- geno.desc[,2] + 2*geno.desc[,3]
    rows.rm <- which(test > Missing)
    if(length(rows.rm) > 0){
      rows.rm <- which(test > Missing)
      ID.rm <- Tab[rows.rm,1]
      Tab <- Tab[-rows.rm,]
      if(Output){write.table(ID.rm,
                            file="Removed_SampleIDs.txt",
                            sep="\t",
                            row.names=F,
                            col.names=F,
                            quote=F)
      }
      rm(ID.rm)
    }
    rm(geno.desc,test)
    if(nrow(Tab)==0){
      Err.Log(Output,"TooMany.Missing")
      stop("Analysis Stopped, too many missing data points",
          call. = F)
    }
  }
}
