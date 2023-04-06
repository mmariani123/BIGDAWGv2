#' File Fetcher
#'
#' Download Protein Alignment and Accessory Files
#' @param Loci HLA Loci to be fetched. Limited Loci available.
#' @note This function is for internal BIGDAWG use only.
GetFiles <- function(Loci){
  #downloads *_prot.txt alignment files
  #downloads hla_nom_p.txt file

  # Get P-Groups Files
  #download.file(
  #  "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/wmda/hla_nom_p.txt",
  #  destfile="hla_nom_p.txt",
  #  method="libcurl")

  # Get Release Version
  df <- readLines(con=paste0("C:/Users/Mike_2/Documents",
                             "/GitHub/BIGDAWGv2/canid_DQB1_prot.txt"),
        n=3)

  df <- readLines(con="canid_DQB1_prot.txt", n=3)
  #df <- readLines(con="dla_nom_p.txt", n=3)
  RD <- unlist(strsplit(df[2],split=" "))[3]
  #RV <- paste(unlist(strsplit(df[3], split=" "))[3:4],collapse=" ")
  RV <- df[3]
  write.table(c(RD,RV),
              file="Dog_Release.txt",
              quote=F,
              col.names=F,
              row.names=F)

  # Get Locus Based Alignments
  for(i in 1:length(Loci)){
    Locus <- Loci[i]
    FileName <- paste0(Locus,"_prot.txt")
    URL <- paste0(
      "https://raw.githubusercontent.com/mmariani123/BIGDAWGv2/latest/alignments",
      FileName)
    #"https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/alignments is 404"
    #But individual files are found, e.g.:
    #https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/alignments/B_prot.txt
    download.file(URL,
                  destfile = FileName,
                  method="libcurl")
  }

}

#' HLA P group File Formatter
#'
#' Format the hla_nom_p.txt read table object for a specific locus.
#' @param x P group object from read.table command.
#' @param Locus Locus to be filtered on.
#' @note This function is for internal BIGDAWG use only.
PgrpFormat_dog <- function(x,Locus){

  # Identify for Locus ... change necessary if DRB
  x.sub <- x[which(x[,1]==Locus),]
  rownames(x.sub) <- NULL
  x.sub[,2] <- sapply(x.sub[,2],
                      function(i){
                        paste(
                          paste0(Locus,
                                 "*",
                                 unlist(strsplit(i,"/"))),
                          collapse="/")
                      }
  )
  colnames(x.sub) <- c("Locus","Allele","P.Group")

  #Expand
  x.list <- list()
  for(i in 1:nrow(x.sub)) {
    if(grepl("/",x.sub[i,'Allele'],fixed=T)){
      tmp <- unlist(strsplit(x.sub[i,'Allele'],"/"))
      tmp <- cbind(rep(Locus,length(tmp)),
                   tmp,
                   rep(x.sub[i,'P.Group'],
                       length(tmp))
      )
      colnames(tmp) <- colnames(x.sub)
      x.list[[i]] <- tmp
    }else{
      x.list[[i]] <- x.sub[i,]
    }
  }
  x.list <- do.call(rbind,x.list)
  x.list <- x.list[order(x.list[,'Allele']),]
  rownames(x.list) <- NULL
  colnames(x.list) <- c("Locus","Allele","P.Group")
  return(x.list)
}

#' HLA P group Finder
#'
#' Identify P group for a given allele if exists.
#' @param x Allele of interest.
#' @param y Formatted P groups.
#' @note This function is for internal BIGDAWG use only.
PgrpExtract <- function(x,y) {
  getRow <- grep(x,y[,'Allele'],fixed=T)
  if(length(getRow)>=1){
    if(length(getRow)>1){
      getRow <- getRow[which(sapply(as.character(y[getRow,'Allele']),
                                    nchar)==nchar(x))]
    }
    return(as.character(y[getRow,'P.Group']))
  }else{
    return("")
  }
}

#' Protein Exon Alignment Formatter
#'
#' Dynamically creates an alignmnet of Allele exons for Analysis.
#' @param Locus Locus alignment to be formatted.
#' @param RefTab Reference exon protein information for alignment formatting.
#' @param Species The species under consideration
#' @note This function is for internal BIGDAWG use only.
ExonPtnAlign.Create <- function(Locus,RefTab,Species){

  if(Species=='hla'){
  #########################################################################
  # Need to remove if DRB split into single locus files
  if(grepl("DRB",Locus)){
    Locus.get <- "DRB"
  }else{
    Locus.get <- Locus
  }
  #########################################################################

  AlignMatrix <- NULL; rm(AlignMatrix)

  #Read in P-Groups
  Pgrps <- read.table("hla_nom_p.txt",
                      fill=T,
                      header=F,
                      sep=";",
                      stringsAsFactors=F,
                      strip.white=T,
                      colClasses="character")
  Pgrps[,1] <- gsub("\\*", "", Pgrps[,1])
  Pgrps <- PgrpFormat(Pgrps,Locus)

  }else if(Species=="dla"){
  Locus.get <- Locus
  browser()
  Pgrps <- read.table(system.file(paste0('extdata/',
                                  Species,
                                  '/',
                                  Species,
                                  '_nom_p.txt'),
                                  package = 'BIGDAWGv2'),
                      fill=T,
                      header=F,
                      sep=";",
                      stringsAsFactors=F,
                      strip.white=T,
                      colClasses="character")
  #browser()
  Pgrps[,1] <- gsub("DLA-|\\*", "", Pgrps[,1])
  Pgrps <- PgrpFormat(Pgrps,Locus)
  }

  #Read in Alignment
  Name <- system.file(paste0('extdata/',
                             Species,
                             '/',
                             Species,
                             '-',
                             Locus.get,
                             '_prot.txt'),
                      package = 'BIGDAWGv2')
  Align <- read.table(Name,
                      fill=T,
                      header=F,
                      sep="\t",
                      stringsAsFactors=F,
                      strip.white=T,
                      colClasses="character")

  #Trim
  Align <- as.matrix(Align[-nrow(Align),]) #Remove Footer

  #Begin Formatting
  Align <- as.matrix(Align[-grep("\\|",Align[,1]),1]) #Remove Pipes
  Align[,1] <- sapply(Align[,1],FUN=sub,pattern=" ",replacement="~")
  Align[,1] <- sapply(Align[,1],FUN=gsub,pattern=" ",replacement="")
  Align <- strsplit(Align[,1],"~")
  Align <- as.matrix(do.call(rbind,Align))

  #browser()

  #Adjust rows to blank where Sequence column == Allele Name
  Align[which(Align[,1]==Align[,2]),2] <- ""

  # Remove Prot Numbering Headers
  Align <- Align[-which(Align[,1]=="Prot"),]

  # Get Unique Alleles
  Alleles <- unique(Align[,1])

  # Loop Through and Build Alignment Block
  Block <- list()
  for(i in Alleles) {
    getRows <- which(Align[,1]==i)
    Block[[i]] <- paste(Align[getRows,2],collapse="")
  }
  Block <- cbind(Alleles,do.call(rbind,Block))

  #Fill end gaps with * to make char lengths even
  Block.len <- max(as.numeric(sapply(Block[,2],FUN=nchar)))
  for( i in 1:nrow(Block) ) {
    Block.miss <- Block.len - nchar(Block[i,2])
    if( Block.miss > 0 ) {
      Block[i,2] <- paste0(as.character(Block[i,2]),
                           paste(rep(".",Block.miss),collapse=""))
    }
  }; rm(i)

  #Split Allele name into separate Locus and Allele, Send Back to Align object
  AlignAlleles <- do.call(rbind,strsplit(Block[,1],"[*]"))
  AlignAlleles <-
    cbind(AlignAlleles,
          apply(AlignAlleles,MARGIN=c(1,2),FUN=GetField,Res=2)[,2])
  rownames(AlignAlleles) <- NULL

  Align <- cbind(AlignAlleles,Block)
  colnames(Align) <- c("Locus",
                       "Allele",
                       "Trimmed",
                       "FullName",
                       "Sequence")

  #Split Sub Alignment into composite elements and
  #Extract relevant positions
  Align.split <- strsplit(Align[,'Sequence'],"")
  Align.split <- do.call(rbind,Align.split)
  AlignMatrix <- cbind(Align[,1:4],Align.split)
  rownames(AlignMatrix) <- NULL

  #browser()
  if(Species=='dla'){
    #truncate ref alleles to 2 positions for now for testing
    for(z in 1:nrow(AlignMatrix)){
      splittedAllele <- unlist(
        strsplit(AlignMatrix[z,'FullName'],
                split='\\*|:'))
      AlignMatrix[z,'FullName'] <-
        paste0(splittedAllele[1],
            '*',
            paste0(splittedAllele[2:3],collapse=":"))
    }
    RefAllele <- paste(RefTab[which(RefTab[,'Locus']==Locus),'Reference.Locus'],
                       RefTab[which(RefTab[,'Locus']==Locus),'Reference.Allele'],
                       sep="*")
  }else{

  #Ensure Reference in Row 1
  RefAllele <- paste(RefTab[which(RefTab[,'Locus']==Locus),'Reference.Locus'],
                     RefTab[which(RefTab[,'Locus']==Locus),'Reference.Allele'],
                     sep="*")
  if(!AlignMatrix[1,'FullName']==RefAllele){

    Align.tmp <- rbind(AlignMatrix[
      which(AlignMatrix[1,'Allele']==RefAllele),],
      AlignMatrix[
        -which(AlignMatrix[1,'Allele']==RefAllele),])
    AlignMatrix <- Align.tmp
    rm(Align.tmp)

  }

  }

  #browser()
  #Save Reference Row
  RefSeq <- AlignMatrix[1,]

  #Ensure Locus Specific Rows
  AlignMatrix <- AlignMatrix[which(AlignMatrix[,'Locus']==Locus),]

  #Rebind Reference if removed (only for DRB3, DRB4, and DRB5)
  if(!AlignMatrix[1,'FullName']==RefAllele){
    AlignMatrix <- rbind(RefSeq,AlignMatrix)
  }

  #Remove columns with no amino acids positions
  #(only for DRB3, DRB4, and DRB5)
  #Count occurence of "." and compare to nrow
  #of AlignMatrix
  rmCol <- which(apply(AlignMatrix,
                       MARGIN=2,
                       {FUN=function(x){
                         length(which(x=="."))
                       }}) == nrow(AlignMatrix)
  )

  if(length(rmCol)>0){
    AlignMatrix <- AlignMatrix[,-rmCol]
  }

  #Propagate Consensus Positions
  for(i in 5:ncol(AlignMatrix)){
    x <- AlignMatrix[,i]
    x[which(x=="-")] <- x[1]
    AlignMatrix[,i] <- x
  }

  #Rename amino acid positions based on reference numbering
  #Deletions are named according to the preceding position with a .1,.2,etc.
  RefStart <- as.numeric(RefTab[which(RefTab[,'Locus']==Locus),
                                'Reference.Start'])
  RefArray <- AlignMatrix[1,5:ncol(AlignMatrix)]
  Names <- NULL ; RefPos <- RefStart
  for(i in 1:length(RefArray) ) {

    if(RefArray[i]==".") {
      Names <- c(Names,
                 paste0("Pos.",
                        RefPos-1,
                        ".",
                        Iteration) )
      Iteration = Iteration + 1
    }else{
      Iteration=1
      Names <- c(Names, paste0("Pos.",RefPos))
      RefPos <- RefPos + 1
      if(RefPos==0){RefPos <- 1}
    }

  }
  colnames(AlignMatrix)[5:ncol(AlignMatrix)] <- Names
  rownames(AlignMatrix) <- NULL

  #Add Absent Allele (Absence due to lack of allele
  #and not lack of typing information)
  AlignMatrix <- rbind(c(Locus,
                         "00:00:00:00",
                         "00:00",
                         paste0(Locus,"*00:00:00:00"),
                         rep("^",ncol(AlignMatrix)-4)),
                       AlignMatrix)

  #Assign P groups
  AlignMatrix <- cbind(AlignMatrix,
                       sapply(AlignMatrix[,'FullName'],
                              PgrpExtract,y=Pgrps) )

  colnames(AlignMatrix)[ncol(AlignMatrix)] <- "P.group"

  #Tally Unknowns as separate column
  AlignMatrix <- cbind(AlignMatrix,
                       apply(AlignMatrix[,5:(ncol(AlignMatrix)-1)],
                             MARGIN=1,
                             FUN=function(x){length(which(unlist(grep("*",x,fixed=T))>0))}
                       )
  )

  colnames(AlignMatrix)[ncol(AlignMatrix)] <- "Unknowns"

  #Tally Null Positions as separate column
  AlignMatrix <- cbind(AlignMatrix,
                       apply(AlignMatrix[,5:(ncol(AlignMatrix)-2)],
                             MARGIN=1,
                             FUN=function(x){
                               length(which(unlist(grep("-",x,fixed=T))>0))}
                       )
  )

  colnames(AlignMatrix)[ncol(AlignMatrix)] <- "NullPositions"

  #Tally InDels
  AlignMatrix <- cbind(AlignMatrix,
                       apply(AlignMatrix[,5:(ncol(AlignMatrix)-3)],
                             MARGIN=1,
                             FUN=function(x){
                               length(which(unlist(grep(".",x,fixed=T))>0))}
                       )
  )

  colnames(AlignMatrix)[ncol(AlignMatrix)] <- "InDels"

  rownames(AlignMatrix) <- NULL
  FileName <-
    paste0(
      devtools::package_file(),
        '/inst/extdata/',
        Species,
        '/ExonPtnAlign_',
        Locus,
        '.obj')
  save(AlignMatrix,
       file=FileName)

}

#' Alignment Object Creator
#'
#' Create Object for Exon Protein Alignments.
#' @param Loci Loci to be bundled.
#' @param Release IMGT/HLA database release version.
#' @param RefTab Data of reference exons used for protein alignment creation.
#' @param Species The species being considered
#' @note This function is for internal BIGDAWG use only.
AlignObj.Create <- function(Loci,Release,RefTab,Species){

  AlignMatrix <- NULL; rm(AlignMatrix)

  ExonPtnList <- list()
  for(i in 1:length(Loci)) {
    Locus <- Loci[i]
    FileName <- system.file(
      paste0('extdata/',
             Species,
             "/ExonPtnAlign_",
             Locus,
             ".obj"),
      package = 'BIDAWGv2')
    load(FileName) #Loads AlignMatrix
    ExonPtnList[[Locus]] <- AlignMatrix
  }

  ExonPtnList[['Release.Version']] <- as.character(Release[2,])
  ExonPtnList[['Release.Date']] <- as.character(Release[1,])
  ExonPtnList[['RefExons']] <- RefTab
  browser()
  save(ExonPtnList,
       file=paste0(devtools::package_file(),
         '/inst/extdata/',
         Species,
         '/ExonPtnAlign.obj'))

}

#' Updated Alignment Object Creator
#'
#' Synthesize Object for Exon Protein Alignments.
#' @param Loci Loci to be bundled.
#' @param Release IMGT/HLA database release version.
#' @param RefTab Data of reference exons used for protein alignment creation.
#' @param Species The Species under consideration
#' @note This function is for internal BIGDAWG use only.
AlignObj.Update <- function(Loci,Release,RefTab,Species){

  AlignMatrix <- NULL; rm(AlignMatrix)

  browser()

  UpdatePtnList <- list()
  for(i in 1:length(Loci)){
    Locus <- Loci[i]
    FileName <- system.file(
      paste0('extdata/',
             Species,
             "/ExonPtnAlign_",
             Locus,
             ".obj"),
      package = 'BIGDAWGv2')
    load(FileName) #Loads AlignMatrix
    UpdatePtnList[[Locus]] <- AlignMatrix
  }
  UpdatePtnList[['Release.Version']] <- as.character(Release[2,])
  UpdatePtnList[['Release.Date']] <- as.character(Release[1,])
  UpdatePtnList[['RefExons']] <- RefTab
  #save(UpdatePtnList,file="UpdatePtnAlign.RData")

  ##MM 04/03/2023
  #I will add the ExonPtnMap field to ExonPtnListbecause
  #I can't seem to find the code where this is implemented
  #For reference look agt BIGDAWG::ExonPtnList$ExonPtnMap

  `12` <- data.frame(Locus=c('12',
                                   '12',
                                   '12',
                                   '12',
                                   '12',
                                   '12',
                                   '12'),
                           Exon=c(1,2,3,4,5,6,7),
                           Start=c(34,34,34,34,34,34,34),
                           Stop=c(80,80,80,80,80,80,80),
                           stringsAsFactors=FALSE)

  `64` <- data.frame(Locus=c('64',
                                    '64',
                                    '64',
                                    '64',
                                    '64',
                                    '64',
                                    '64'),
                            Exon=c(1,2,3,4,5,6,7),
                            Start=c(34,34,34,34,34,34,34),
                            Stop=c(80,80,80,80,80,80,80),
                            stringsAsFactors=FALSE)

  `79` <- data.frame(Locus=c('79',
                                    '79',
                                    '79',
                                    '79',
                                    '79',
                                    '79',
                                    '79'),
                            Exon=c(1,2,3,4,5,6,7),
                            Start=c(34,34,34,34,34,34,34),
                            Stop=c(80,80,80,80,80,80,80),
                            stringsAsFactors=FALSE)

  `88` <- data.frame(Locus=c('88',
                                    '88',
                                    '88',
                                    '88',
                                    '88',
                                    '88',
                                    '88'),
                            Exon=c(1,2,3,4,5,6,7),
                            Start=c(34,34,34,34,34,34,34),
                            Stop=c(80,80,80,80,80,80,80),
                            stringsAsFactors=FALSE)

  DQA1 <- data.frame(Locus=c('DQA1',
                                    'DQA1',
                                    'DQA1',
                                    'DQA1',
                                    'DQA1',
                                    'DQA1',
                                    'DQA1'),
                            Exon=c(1,2,3,4,5,6,7),
                            Start=c(34,34,34,34,34,34,34),
                            Stop=c(80,80,80,80,80,80,80),
                            stringsAsFactors=FALSE)

  DQB1 <- data.frame(Locus=c('DQB1',
                                    'DQB1',
                                    'DQB1',
                                    'DQB1',
                                    'DQB1',
                                    'DQB1',
                                    'DQB1'),
                            Exon=c(1,2,3,4,5,6,7),
                            Start=c(34,34,34,34,34,34,34),
                            Stop=c(80,80,80,80,80,80,80),
                            stringsAsFactors=FALSE)

  DRB1 <- data.frame(Locus=c('DRB1',
                                    'DRB1',
                                    'DRB1',
                                    'DRB1',
                                    'DRB1',
                                    'DRB1',
                                    'DRB1'),
                            Exon=c(1,2,3,4,5,6,7),
                            Start=c(34,34,34,34,34,34,34),
                            Stop=c(80,80,80,80,80,80,80),
                            stringsAsFactors=FALSE)

  UpdatePtnList[['ExonPtnMap']] <- list(
    '12'=`12`,
    '64'=`64`,
    '79'=`79`,
    '88'=`88`,
    'DQA1'=DQA1,
    'DQB1'=DQB1,
    'DRB1'=DRB1
  )

  browser()
  save(UpdatePtnList,
       file=paste0(
         devtools::package_file(),
         '/inst/extdata/',
         Species,
         '/UpdatePtnAlign.RData'))

  ##Example, look at human:
  ##BIGDAWG::ExonPtnList$ExonPtnMap

  #$A
  #Locus Exon Start Stop
  #[1,] "A"   "1"  "1"   "24"
  #[2,] "A"   "2"  "25"  "114"
  #[3,] "A"   "3"  "115" "206"
  #[4,] "A"   "4"  "207" "298"
  #[5,] "A"   "5"  "299" "337"
  #[6,] "A"   "6"  "338" "348"
  #[7,] "A"   "7"  "349" "364"
  #[8,] "A"   "8"  "365" "365"

  #$B
  #Locus Exon Start Stop
  #[1,] "B"   "1"  "1"   "24"
  #[2,] "B"   "2"  "25"  "114"
  #[3,] "B"   "3"  "115" "206"
  #[4,] "B"   "4"  "207" "298"
  #[5,] "B"   "5"  "299" "337"
  #[6,] "B"   "6"  "338" "348"
  #[7,] "B"   "7"  "349" "362"

  #$C
  #Locus Exon Start Stop
  #[1,] "C"   "1"  "1"   "24"
  #[2,] "C"   "2"  "25"  "114"
  #[3,] "C"   "3"  "115" "206"
  #[4,] "C"   "4"  "207" "298"
  #[5,] "C"   "5"  "299" "338"
  #[6,] "C"   "6"  "339" "349"
  #[7,] "C"   "7"  "350" "365"
  #[8,] "C"   "8"  "366" "366"

  #$DPA1
  #Locus  Exon Start Stop
  #[1,] "DPA1" "1"  "1"   "33"
  #[2,] "DPA1" "2"  "34"  "115"
  #[3,] "DPA1" "3"  "116" "209"
  #[4,] "DPA1" "4"  "210" "260"

  #$DPB1
  #Locus  Exon Start Stop
  #[1,] "DPB1" "1"  "1"   "33"
  #[2,] "DPB1" "2"  "34"  "121"
  #[3,] "DPB1" "3"  "122" "215"
  #[4,] "DPB1" "4"  "216" "252"
  #[5,] "DPB1" "5"  "253" "258"

  #$DQA1
  #Locus  Exon Start Stop
  #[1,] "DQA1" "1"  "1"   "27"
  #[2,] "DQA1" "2"  "28"  "110"
  #[3,] "DQA1" "3"  "111" "204"
  #[4,] "DQA1" "4"  "205" "255"

  #$DQB1
  #Locus  Exon Start Stop
  #[1,] "DQB1" "1"  "1"   "36"
  #[2,] "DQB1" "2"  "37"  "126"
  #[3,] "DQB1" "3"  "127" "220"
  #[4,] "DQB1" "4"  "221" "257"
  #[5,] "DQB1" "5"  "258" "261"

  #$DRB1
  #Locus  Exon Start Stop
  #[1,] "DRB1" "1"  "1"   "33"
  #[2,] "DRB1" "2"  "34"  "123"
  #[3,] "DRB1" "3"  "124" "217"
  #[4,] "DRB1" "4"  "218" "254"
  #[5,] "DRB1" "5"  "255" "262"
  #[6,] "DRB1" "6"  "263" "266"

  #$DRB3
  #Locus  Exon Start Stop
  #[1,] "DRB3" "1"  "1"   "33"
  #[2,] "DRB3" "2"  "34"  "123"
  #[3,] "DRB3" "3"  "124" "217"
  #[4,] "DRB3" "4"  "218" "254"
  #[5,] "DRB3" "5"  "255" "262"
  #[6,] "DRB3" "6"  "263" "266"

  #$DRB4
  #Locus  Exon Start Stop
  #[1,] "DRB4" "1"  "1"   "33"
  #[2,] "DRB4" "2"  "34"  "123"
  #[3,] "DRB4" "3"  "124" "217"
  #[4,] "DRB4" "4"  "218" "254"
  #[5,] "DRB4" "5"  "255" "262"
  #[6,] "DRB4" "6"  "263" "266"

  #$DRB5
  #Locus  Exon Start Stop
  #[1,] "DRB5" "1"  "1"   "33"
  #3[2,] "DRB5" "2"  "34"  "123"
  #[3,] "DRB5" "3"  "124" "217"
  #[4,] "DRB5" "4"  "218" "254"
  #[5,] "DRB5" "5"  "255" "262"
  #[6,] "DRB5" "6"  "263" "266"

}
