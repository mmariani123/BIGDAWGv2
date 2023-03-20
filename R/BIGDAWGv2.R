#' BIGDAWGv2 Main Wrapper Function
#'
#' This is the main wrapper function for each analysis.
#' @param Data Name of the genotype data file.
#' @param HLA Logical Indicating whether data is HLA class I/II genotyping data only.
#' @param Species Select species string: 'hla','dla','gla','cla'
#' @param Run.Tests Specifics which tests to run.
#' @param Loci.Set Input list defining which loci to use for analyses (combinations permitted).
#' @param Exon Numeric Exon(s) for targeted amino acid analysis.
#' @param All.Pairwise Logical indicating whether all pairwise loci should be analyzed in haplotype analysis.
#' @param Trim Logical indicating if HLA alleles should be trimmed to a set resolution.
#' @param Res Numeric setting what desired resolution to trim HLA alleles.
#' @param EVS.rm Logical indicating if expression variant suffixes should be removed.
#' @param Missing Numeric setting allowable missing data for running analysis (may use "ignore").
#' @param Strict.Bin Logical specify if strict rare cell binning should be used in ChiSq test.
#' @param Cores.Lim Integer setting the number of cores accessible to BIGDAWG (Windows limit is 1 core).
#' @param Results.Dir Optional, string of full path directory name for BIGDAWG output.
#' @param Return Logical Should analysis results be returned as list.
#' @param Output Logical Should analysis results be written to output directory.
#' @param Merge.Output Logical Should analysis results be merged into a single file for easy access.
#' @param Verbose Logical Should a summary of each analysis be displayed in console.
#' @examples
#' \dontrun{
#' ### The following examples use the synthetic data set bundled with BIGDAWG
#'
#' # Haplotype analysis with no missing genotypes for two loci sets
#' # Significant haplotype association with phenotype
#' # BIGDAWG(Data="HLA_data", Run.Tests="H", Missing=0, Loci.Set=list(c("DRB1","DQB1")))
#'
#' # Hardy-Weinberg and Locus analysis ignoring missing data
#' # Significant locus associations with phenotype at all but DQB1
#' # BIGDAWG(Data="HLA_data", Run.Tests="L", Missing="ignore")
#'
#' # Hardy-Weinberg analysis trimming data to 2-Field resolution with no output to files (console only)
#' # Significant locus deviation at DQB1
#' BIGDAWG(Data="HLA_data", Run.Tests="HWE", Trim=TRUE, Res=2, Output=FALSE)
#' }
BIGDAWGv2 <- function(Data,
         HLA=TRUE,
         Species='hla',
         Run.Tests,
         Loci.Set,
         Exon,
         All.Pairwise=FALSE,
         Trim=FALSE,
         Res=2,
         EVS.rm=FALSE,
         Missing=2,
         Strict.Bin=FALSE,
         Cores.Lim=1L,
         Results.Dir,
         Return=FALSE,
         Output=TRUE,
         Merge.Output=FALSE,
         Verbose=TRUE){

  options(warn=-1)

  MainDir <- getwd()

  on.exit(setwd(MainDir), add = TRUE)

  ################## Check if data present: #############

  if(missing(Data)){
    Err.Log("P.Missing","Data") ;
    stop("Analysis Stopped. Missing Data",call.=FALSE)
  }

  ################ Check Parameters #####################

  HLA <- as.logical(HLA)

  Check.Params(HLA,
               Loci.Set,
               Exon,
               All.Pairwise,
               Trim,
               Res,
               EVS.rm,
               Missing,
               Cores.Lim,
               Return,
               Output,
               Merge.Output,
               Verbose)

  ######## CHECK MULTICORE LIMITATIONS #################

  Cores <- BIGDAWG::Check.Cores(Cores.Lim,Output)

  cat(rep("=",40))
  cat("\n       BIGDAWG: Bridging ImmunoGenomic Data Analysis Workflow Gaps\n")
  cat(rep("=",40),"\n")
  cat("\n>>>>>>>>>>>>>>>>>>>>>>>> BEGIN Analysis <<<<<<<<<<<<<<<<<<<<<<<<\n\n")

  # Define Output object
  BD.out <- list()

  ###########################################################
  #### ================================================= ####
  #### ______________ Read in Data _____________________ ####
  #### ================================================= ####
  ###########################################################

  NAstrings=c("NA","","****","-","na","Na")

  if(is.character(Data)){

    if(Data=="HLA_data"){

      # Using internal synthetic set
      Tab <- BIGDAWG::HLA_data
      Data.Flag <- "Internal Synthetic Data Set"

    }else{

      # Read in data file entered as string
      if(!file.exists(Data)){
        Err.Log(Output,"Bad.Filename", Data)
        stop("Analysis stopped.",call.=F)
      }
      Tab <- read.table(Data,
                        header = T,
                        sep="\t",
                        stringsAsFactors = F,
                        na.strings=NAstrings,
                        fill=T,
                        comment.char = "#",
                        strip.white=T,
                        blank.lines.skip=T,
                        colClasses="character")
      Data.Flag <- Data

    }

  }else{

    # Using R object
    Tab <- Data
    Data.Flag <- deparse(substitute(Data))

    # Convert Empty Cells to NA
    for(i in 3:ncol(Tab)){
      putCell <- which(sapply(Tab[,i],nchar)==0)
      if(length(putCell) > 0 ){Tab[putCell,i] <- NA}
    }

  }

  ###########################################################
  #### ================================================= ####
  #### _________ Declare Data Input Parameter __________ ####
  #### ================================================= ####
  ###########################################################

  cat("Data Input: ",Data.Flag,"\n\n\n")

  # Convert GLS data
  if(ncol(Tab) == 3){GLSFLAG = TRUE} else {GLSFLAG = FALSE}
  if(GLSFLAG && !HLA){Err.Log(Output,"notHLA.GLS")}
  if(GLSFLAG && HLA){
    cat("Converting Gene List Strings to Tabular Format...\n\n")
    Tab <- GLSconvert(Tab,Convert="GL2Tab",
                      System="HLA",
                      File.Output="R",
                      Strip.Prefix=T,
                      Abs.Fill=T,
                      Cores.Lim=Cores)
  }

  ###########################################################
  #### ================================================= ####
  #### ___________________ Prep Data ___________________ ####
  #### ================================================= ####
  ###########################################################

  # Prep Data for processing and checks
  Tab <- prepData(Tab)

  ###########################################################
  #### ================================================= ####
  #### ___________ Set Output Directory ________________ ####
  #### ================================================= ####
  ###########################################################

  # Define and Change to the required output directory

  if(Output){
    if(missing(Results.Dir)){
      OutDir <- paste0(MainDir,
                       .Platform$file.sep,
                      "output ",
                      format(Sys.time(), "%d%m%y %H%M%S")
                      )
      dir.create(OutDir)
    } else {
      OutDir <- Results.Dir
    }
  }
  if(Output){setwd(OutDir)}

###########################################################
#### ================================================= ####
#### _____ Data Processing and Sanity Checks _________ ####
#### ================================================= ####
###########################################################

  cat(">>>> DATA PROCESSING AND CHECKS.\n")

#### General processing and checks for all data

  # Define Data Columns
  Data.Col <- seq(3,ncol(Tab))

  # RUN TESTS DEFINITIONS
  if(missing(Run.Tests)){

    Run <- c("HWE","H","L","A")

  }else{

    Run <- Run.Tests

  }

  if("A" %in% Run){
    if(Species=='hla'){

    }else if(Species=='dla'){

    }else if(Species=='gla'){

    }else if(Species=='cla'){

    }else{
      stop(paste0('Innapropriate <Species> parameter, ',
                  'program terminating ...'))
    }
  }else{
    cat(paste0('Amino Acid test not selected, ',
               'test will not be run.\n'))
  }

############## BAD DATA DEFINITIONS - No 1's or 0's######

BIGDAWGv2::bad_data_def(Tab=Tab,
                        Data.Col=Data.Col,
                        Output=Output)

###################### MISSING DATA ##############################

rows.rm <- BIGDAWGv2::missing_data_check(Missing=Missing,
                                         Run=Run,
                                         Output=Output,
                                         Tab=Tab,
                                         Data.Col=Data.Col,
                                         NAstrings=NAstrings)

######### MULTIPLE SETS AND ANALYSIS DUPLICATION ###########
############################################################
############################################################
############################################################
############################################################

BIGDAWGv2::mult_set_dup_check(Loci.Set,
                              All.Pairwise,
                              Run,
                              Output)

############ DATA MERGE AND NUMBER OF LOCI ###############
##########################################################
##########################################################
##########################################################
##########################################################

BIGDAWGv2::data_merge_num_loci_check(Output,
                                     Merge.Output,
                                     All.Pairwise,
                                     Tab)

#################### HLA Specific Checks #################
##########################################################
##########################################################
##########################################################
##########################################################

if(Species=='hla'){
  #THe position list will have to be specific to the species
  #ie Update ptn list
  #and the position list has to do with the core exons

  #Check for the updated ExonPtnList 'UpdatePtnList' and
  #use if found.
  UpdatePtnList <- NULL
  UPL <- paste0(path.package('BIGDAWG'),
              "/data/UpdatePtnAlign.RData")
  if(file.exists(UPL)) {
    load(UPL)
    EPL <- UpdatePtnList
    rm(UpdatePtnList)
    UPL.flag=T
  }else{
    rm(UpdatePtnList,UPL)
    EPL <- BIGDAWG::ExonPtnList
    UPL.flag=F
  }

  if(Trim & !HLA){Err.Log(Output, "NotHLA.Trim")}
  if(EVS.rm & !HLA){Err.Log(Output, "NotHLA.EVS.rm")}
  if(!HLA){
    DRBFLAG <- NULL
  }else{
    DRB345.test <- length(grep("DRB345",colnames(Tab)))>0
  }
}

################## What is below for ? #####################
############################################################
############################################################
############################################################
############################################################

if(HLA){

  runHlaCheckOutput <- run_hla_checks(Trim=Trim,
                            EVS.rm=EVS.rm,
                            Run=Run,
                            DRB345.test=DRB345.test,
                            Output=Output,
                            Cores=Cores,
                            Res=Res)

  Set <- runHlaCheckOutput[[1]]
  Release <- runHlaCheckOutput[[2]]

}

###########################################################
#### ================================================= ####
#### _____________ Case-Control Summary ______________ ####
#### ================================================= ####
###########################################################

case_control_summary(Trim=Trim,
                     Tab=Tab,
                     Res=Res,
                     HLA=HLA,
                     Verbose=Verbose,
                     Output=Output)

###########################################################
#### ================================================= ####
#### ____________ Write to Parameter File ____________ ####
#### ================================================= ####
###########################################################

if(Output){

    if(HLA && !is.null(DRBFLAG)){
      DRB345.tmp <- DRBFLAG
    }else{
      DRB345.tmp <- NULL
    }
    if(HLA){
      Trim.tmp <- Trim
    }else{
      Trim.tmp <- NULL
    }
    if(HLA && Trim){
      Res.tmp <- Res
    }else{
      Res.tmp <- NULL
    }
    if(HLA){
      EVS.rm.tmp <- EVS.rm
    }else{EVS.rm.tmp <- NULL
    }

    if(!missing(Exon)){
      Exon.tmp <- paste(unique(unlist(Exon)),collapse=",")
    }else{
      Exon.tmp <- NULL
    }

    Params.Run <- list(Time = format(
                        Sys.time(), "%a %b %d %X %Y"),
                       BD.Version = as.character(
                                      packageVersion("BIGDAWG")),
                       Cores.Used = Cores,
                       File = Data.Flag,
                       Output.Results = Output,
                       Merge = Merge.Output,
                       Return.Object = Return,
                       Display.Results = Verbose,
                       HLA.Data = HLA,
                       Exon = Exon.tmp,
                       DRB345.Parsed = DRB345.tmp,
                       Tests = paste(Run,collapse=","),
                       All.Pairwise = All.Pairwise,
                       Trim = Trim.tmp,
                       Resolution = Res.tmp,
                       Suffix.Stripping = EVS.rm.tmp,
                       Missing.Allowed = Missing,
                       Strict.Binning = Strict.Bin,
                       Samples.Removed = length(rows.rm))

    Params.Run <- do.call(rbind,Params.Run)
    write.table(Params.Run,
                file="Run_Parameters.txt",
                sep=": ",
                row.names=T,
                col.names=F,
                quote=F)
  }

###########################################################
#### ================================================= ####
#### ________ Hardy Weinberg Equilibrium 'HWE' _______ ####
#### ================================================= ####
###########################################################

if("HWE" %in% Run){

  BD.out <- run_hwe_analysis(HLA=HLA,
                             Trim=Trim,
                             Tab=Tab,
                             Output=Output,
                             Verbose=Verbose,
                             BD.out=BD.out)

} #END HARDY-WEINBERG

###########################################################
#### ================================================= ####
#### ________________ Set Loop Begin _________________ ####
#### __ (loop through each defined locus/loci set) ___ ####
#### ================================================= ####
###########################################################

if(sum(c("H","L","A") %in% Run) > 0){

  cat("\n>>>>>>>>>>>> Begin Locus Sets <<<<<<<<<<<\n\n")

  if(length(Set)==1){

    cat("Your analysis has 1 set to analyze.\n")

  }else{

    cat(paste("Your analysis has ",
              length(Set),
              " sets to analyze.",
              sep=""),
        "\n")
  }

  for(k in 1:length(Set)){
    cat("\n")
    cat(paste(rep(">",35),collapse=""),"Running Set",k,"\n")

    cols <- Set[[k]]
    Tabsub <- Tab[,c(1,2,cols)]

    #Set Specific Global Variables
    SID <- Tabsub[,1] # sample IDs
    genos <- Tabsub[,3:ncol(Tabsub)] # genotypes
    genos[genos==""] <- NA
    grp <- Tabsub[, 2] # phenotype
    #nGrp0 <- length(which(grp==0))*2 #nalleles
    #nGrp1 <- length(which(grp==1))*2 #nalleles
    loci <- unique(gsub(".1","",colnames(genos),fixed=T)) # name of loci
    loci.ColNames <- gsub(".1","",colnames(genos),fixed=T) # column names
    nloci <- as.numeric(length(loci)) # number of loci
    SetName <- paste('Set',k,sep="")

    if(HLA==T){genos[genos=='^'] <- "00:00"}

    if(Output){

      OutSetDir <- paste0(OutDir,
                          .Platform$file.sep,
                          "Set",
                          k)
      dir.create(OutSetDir)
      setwd(OutSetDir)

      Params.set <- list(Set = paste("Set",k),
                         Loci.Run = paste(loci,collapse=",")
                         )

      Params.set <- do.call(rbind,Params.set)
      write.table(Params.set,
                  file="Set_Parameters.txt",
                  sep=": ",
                  row.names=T,
                  col.names=F,
                  quote=F)
      }

    SAFE <- c(ls(),"SAFE")


###########################################################
#### ================================================= ####
#### __________ Haplotype Analysis 'H'  ______________ ####
#### ================================================= ####
###########################################################

if("H" %in% Run){

  output <- run_haplotype_analysis(
              nloci=nloci,
              All.Pairwise=All.Pairwise,
              SID=SID,
              Tabsub=Tabsub,
              loci=loci,
              loci.ColNames=loci.ColNames,
              genos=genos,
              grp=grp,
              Strict.Bin=Strict.Bin,
              Output=Output,
              Verbose=Verbose,
              Cores=Cores,
              SAFE=SAFE,
              BD.out=BD.out)

  BD.out <- output[[1]]
  SAFE <- output[[2]]

}

###########################################################
#### ================================================= ####
#### ________________ Locus Level 'L' ________________ ####
#### ================================================= ####
###########################################################

if("L" %in% Run){

BD.out <- run_locus_analysis(
            nloci=nloci,
            loci=loci,
            loci.ColNames=loci.ColNames,
            genos=genos,
            grp=grp,
            Strict.Bin=Strict.Bin,
            Output=Output,
            Verbose=Verbose,
            SetName=SetName,
            SAFE=SAFE,
            BD.out=BD.out)

  BD.out <- output[[1]]
  SAFE <- output[[2]]

}

###########################################################
#### ================================================= ####
#### _____________ Amino Acid Level 'A' ______________ ####
#### ================================================= ####
###########################################################

if("A" %in% Run){

  if(Species=='hla'){
    BD.out <- run_amino_acid_analysis_hla(
                UPL.flag=UPL.flag,
                nloci=nloci,
                All.Pairwise=All.Pairwise,
                SID=SID,
                Tabsub=Tabsub,
                loci=loci,
                loci.ColNames,
                genos=genos,
                grp=grp,
                Exon=Exon,
                EPL=EPL,
                Cores=Cores,
                Strict.Bin=Strict.Bin,
                Output=Output,
                Verbose=Verbose,
                Release=Release,
                SetName=SetName,
                SAFE=SAFE,
                BD.out=BD.out)
    BD.out <- output[[1]]
    SAFE <- output[[2]]

  }else if(Species=='dla'){

  }else if(Species=='cla'){

  }else if(Species=='gla'){

  }else{
    stop("Must select an appropriate species, program terminating ...")
  }

} #END AMINO ACID

###########################################################
#### ================================================= ####
#### ___________________ End Analysis ________________ ####
#### ================================================= ####
###########################################################

}

rm(k)

}# END SET LOOP

if(Output){

  if(Merge.Output){

    merge_output(BD.out=BD.out,
                 Run=Run,
                 OutDir=OutDir)

  }

}

BD.out <- end_analysis(Output=Output,
                       OutDir=OutDir,
                       BD.out=BD.out,
                       Return=Return)

} # END FUNCTION
