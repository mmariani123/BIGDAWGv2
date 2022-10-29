#' case_control_summary
#'
#' case_control_summary function
#' @param HLA Logical Indicating whether data is HLA class I/II genotyping data only.
#' @param Trim Logical indicating if HLA alleles should be trimmed to a set resolution.
#' @param Tab Input data to be worked on
#' @param Res Numeric setting what desired resolution to trim HLA alleles.
#' @param Output Logical Should analysis results be written to output directory.
#' @param Verbose Logical Should a summary of each analysis be displayed in console.

# =============================================================== ####
# Case-Control Summary __________________________________________ ####

case_control_summary <- function(Trim,
                                 Res,
                                 Tab,
                                 HLA,
                                 Verbose,
                                 Output
                                ){

  cat("\n>>>> CASE - CONTROL SUMMARY STATISTICS\n")
  #cat(paste(rep("_",50),collapse=""),"\n")
  if(Trim){
    rescall <- paste(Res,"-Field",sep="")
  }else{
    rescall <- "Not Defined"
  }
  Check <- PreCheck(Tab,
                    colnames(Tab),
                    rescall,
                    HLA,
                    Verbose,
                    Output)
  if(Output){
    write.table(Check,
                file="Data_Summary.txt",
                sep=": ",
                col.names=F,
                row.names=T,
                quote=F)
    rm(Check,rescall)
  }
}
