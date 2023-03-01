#!/usr/bin/env Rscript

#' get_files
#'
#' functions to download files, used to create alnprot.txt _nom_p.txt files
#' @param fileUrl the Url of the file, not including the file name
#' @param fileDlNames the file name
#' @param fileOutNames the desired names for the downloaded files
#' @param species the desired specia, e.g., "dla", will affect which
#' folder the files are downloaded to.
#' @note This function is for general use
#' @examples
#' \dontrun{
#' # The following is an example:
#' # get_files(fileUrl=fileUrl,
#'   #fileDlNames=fileDlNames,
#'  fileOutNames=fileOutNames)
#' }
get_files <- function(fileUrl,
                      fileDlNames,
                      fileOutNames,
                      species){
  #downloads files such as *_prot.txt alignment files,
  #_nom_p.txt files, .fasta files, etc. for a select species.
  #setwd("./data")
  lapply(fileDlNames,
        FUN = function(x){
              download.file(
              url = paste0(fileUrl, x),
              destfile=paste0("./data/",species,"/",x),
              method="libcurl")
          }
        )
}

#' make_prot_file
#'
#' Function to create the .alnprot.txt files from the
#' files from the '-protein.aln' files hosted on GitHub
#' @param fileIn Input file paths
#' @param groupSize the size of the chunks of amino acid strings
#' the default is 10, and they are output as tab-separated columns
#' the output files(s)
#' @param allBases if set to FALSE, then only mismatched bases
#' will be included and the remainder will be set to '-' .
#' @param species human: hla , dog: dla, cow: bla, chicken: cla
#' @param pos the offset position to account for leader peptides for
#' the gene being anaylyzed.
#' @note This function is for general use
#' calc_spaces
#'
#' This function is for positioning the pos markers and '|' indicators
#' in the _nom_p.txt files
#' @param x
#' @param groupSize
#' @param allBases
#' @param species
#' @param pos
#' @param lengths The length of the AA sequences corresponding to the
#' third column in the download -protein align files.
#' @param ncharName the length in characters of the allele nanmes
#' corresponding to the first field in the downloaded -protein.aln files
#' @note This function is for general use
#' @examples
#' \dontrun{
#' # The following is an example:
#' # get_files(fileUrl=fileUrl,
#'   #fileDlNames=fileDlNames,
#'  fileOutNames=fileOutNames)
#' DRB1 starts at -29,
#' DQA1 starts at -23,
#' DQB1 starts at -32,
#' DPA1 starts at -31
#' and DPB1 starts at -29,
#' so the protein alignments
#' all start at -30 or -40
#'
#' posVec <- c(
#'  -30,
#'  1,
#'  60,
#'  120,
#'  180,
#'  240,
#'  300,
#'  360,
#'  381)
#'
#' dla.prot.files <- list.files(path = "./data/dla",
#'                             pattern = "dla",
#'                             full.names = TRUE)
#'
#' make_prot_file(lengths = posVec,
#'            nCharName = 13)
#'}

make_prot_file <- function(x = fileOutNamesProt[3],
                        groupSize = 10,
                        allBases = FALSE,
                        species = "dla",
                        pos=29,
                        lengths,
                        nCharName){

  library(data.table)
  library(dplyr)
  library(magrittr)
  library(plyr)
  library(abind)

  #outFile <- "C:\\Users\\mmari\\OneDrive\\Desktop\\testOut.txt"
  #if(file.exists(outFile)){
  #  unlink(outFile)
  #}

  fullPath <- paste0("./data/",species,"/",x)
  print(paste0("Processing file ",fullPath))
  df <- read.table(fullPath,
                   sep="",
                   header = FALSE,
                   stringsAsFactors = FALSE)
  groupNames   <- df[,1]
  strings      <- df[,2]
  protLengths  <- df[,3]
  uniqProtLengths <- as.numeric(unique(protLengths))
  groupNamesFinal <- c()
  #print(colnames(df))
  casted <- aggregate(df, V2 ~ V1, FUN=paste, collapse = "")

  if(length(unique(nchar(casted$V2)))!=1){
    dupEntries <- which(unique(nchar(casted$V2)!=max(uniqProtLengths)))
    stop(paste0('There are duplicated entries in the file!',
                ' near lines ',
                paste(
                  as.character(
                    which(nchar(casted$V2)!=max(uniqProtLengths))
                  ),
                  collapse=",")
    )
    )
  }
  #print(casted)
  #stop()
  dfMat <- data.frame()
  #for(i in 1:length(uniqProtLengths)){
  #dfNow <- df[which(df[,3]==uniqProtLengths[i]),]
  if(allBases==FALSE){
    for(j in 1:nrow(casted)){
      if(j!=1){
        #print(dfNow[1,2])
        #print(dfNow[j,2])
        matchInd <- mapply(function(x, y) which(x == y),
                           strsplit(casted[1,2], split=""),
                           strsplit(casted[j,2], split=""))
        if(length(matchInd[[1]])!=0){
          for(k in seq_along(matchInd)){
            substring(casted[j,2],
                      matchInd[k],
                      matchInd[k]) <- "-"
          }
        }
      }
    }
  }
  #account for leader peptide here:
  leadPadding <- plyr::round_any(pos, 10, f=ceiling)
  leaderPad <- leadPadding - pos
  repDots <- paste(rep(" ", times=leaderPad), collapse="")
  pepNow <- paste0(repDots, casted[,2])
  #print(pepNow)

  choppedStrings <- array(rep("NA", times=10),
                          c(nrow=nrow(casted),
                            ncol=ceiling((max(uniqProtLengths)+leaderPad)/groupSize),
                            1)
  )

  pepValues <-
    unlist(
      strsplit(
        pepNow,
        split = paste0("(?<=.{",groupSize,"})"),
        perl = TRUE
      )
    )

  nCol <- ceiling((max(uniqProtLengths)+leaderPad)/groupSize)
  nRow <- nrow(casted)

  #Don't forget to add the leader values to the pepValues below:
  pepValues <- paste0(pepValues,rep("",nCol*nRow-length(pepValues)))
  pepValues <- c(
    paste(rep(' ',
              times=leaderPad,
              collapse = "")),
    pepValues)

  for(i in 1:length(pepValues)){
    choppedStrings[ceiling(i/nCol),ifelse(i%%nCol!=0,i%%nCol,nCol),] <-
      pepValues[i]
    #print(i)
  }

  dfMat <- choppedStrings
  dfMat <- abind::abind(
    array(
      rep(unique(df$V1)),
      c(nrow=nrow(dfMat),ncol=1,1)
    ),
    dfMat,
    along = 2)

  #lengths <- c(as.integer(-pos)-leaderPad, as.integer(-pos), 1, unique(df[,3]))
  #lengths <- c(as.integer(-pos)-leaderPad, 1, unique(df[,3]))
  lengths = posVec

  outFile <- paste0("./data/",species,"/",x,"prot.txt")

  #Check if file exists, if it does erase it:
  if(file.exists(outFile)){
    unlink(outFile)
  }

  dfMat1 <- dfMat[,1,1]
  dfMat2 <- dfMat[,-1,1]

  if(file.exists(outFile)){
    unlink(outFile)
  }

  #Create file header (first lines)
  protHeader <- paste0(
    '# file: ',outFile,'\n',
    '# date: ','2023-02-06\n',
    '# version: BIGDAWG 2.0\n',
    '# origin: github.com/mmariani123/bigdawgv2\n',
    '# repository: github.com/mmariani123/bigdawgv2\n',
    '# author: Mariani Systems LLC, Michael P. Mariani ',
    '(m.mariani123@gmail.com)\n'
  )

  #Create output file and add header:
  fileConn <- file(outFile)
  writeLines(protHeader, fileConn)
  close(fileConn)

  setsList <- list()
  negStart <- lengths[1]
  #for(i in seq_along(1:ceiling(nCol/10))){
  for(i in seq_along(1:ceiling(max(lengths)/100))){
    fileConn <- file(outFile,open="a")
    print(paste0('current iter is: ',as.character(i)))
    lengthsNow <-
      lengths[
        which(
          ((lengths-negStart)<=i*100 & (lengths-negStart>(i-1)*100)) |
            ((lengths-negStart)<=i*100 & i==1 & (lengths-negStart==0))
        )
      ]
    print(paste("The current positions are: ",
                paste(as.character(lengthsNow),
                      collapse=" "), collapse=" "))

    setsList[[i]] <- lengthsNow

    posLine  <- list()
    vertLine <- list()
    for(j in 1:length(lengthsNow)){
      if(j==1 & i ==1){
        posLine[[j]] <-
          paste('prot',
                paste(rep(" ",
                          times=nCharName-nchar('prot')-1),
                      collapse=""),
                lengths[1],
                collapse="")

        vertLine[[j]] <-
          paste('prot',
                paste(rep(" ",
                          times=nCharName-nchar('prot')-1),
                      collapse=""),
                '|',
                collapse="")

      }
      else if(i==1 & j!=1){
        #consider the neg nature of the first peptide seq:
        interDist <- abs(lengths[j]-lengths[j-1])-nchar(lengths[j-1])-ifelse(j==2,3,0)
        print(paste0('interDist = ',as.character(interDist)))
        #Calculate the between column spaces becasue each column is grouped
        #into 10 AA currently
        bcs <- ceiling((abs(lengths[j]-lengths[j-1])/10))
        print(paste0('between column spaces = ',as.character(bcs)))

        posLine[[j]] <- paste(
          paste(rep(" ",times=interDist+bcs), collapse=""),
          lengths[j], collapse="")

        vertLine[[j]] <- paste(
          paste(rep(" ",times=interDist+bcs+ifelse(j==2,2,0)), collapse=""),
          '|', collapse="")

      }
      else if(i!=1 & j==1){
        #remaining from last line:
        interDist <-
          setsList[[i]][1] -
          lengths[length(unlist(setsList[1:i-1]))] -
          ((i-1)*100 - (lengths[length(unlist(setsList[1:i-1]))] + abs(negStart))) +
          #0
          -3
        #(nchar('prot')-1) -
        #6
        #(nCharName)
        #4
        #14
        #4 -
        #10
        #(nCharName-nchar('prot')) -
        #nchar('prot') #+
        #nchar(lengths[length(unlist(setsList[1:i-1]))]) #+
        #(nCharName - nchar('prot')-1)
        print(paste0('interDist = ',as.character(interDist)))
        #Calculate the between column spaces becasue each column is grouped
        #into 10 AA currently
        bcs <- ceiling((((setsList[[i]][1] -
                            lengths[length(unlist(setsList[1:i-1]))]) -
                           ((i-1)*100 - (lengths[length(unlist(setsList[1:i-1]))] + abs(negStart)))
        )/10))
        print(paste0('between column spaces = ',as.character(bcs)))

        posLine[[j]] <- paste('prot',
                              paste(rep(" ",times=interDist+bcs),collapse=""),
                              setsList[[i]][1],collapse="")

        vertLine[[j]] <- paste('prot',
                               paste(rep(" ",times=interDist+bcs),collapse=""),
                               '|',collapse="")

      }else{
        interDist <- setsList[[i]][j] -
          setsList[[i]][j-1] -
          7
        #nchar(setsList[[i]][j-1])
        print(paste0('interDist = ',as.character(interDist)))
        #Calculate the between column spaces becasue each column is grouped
        #into 10 AA currently
        bcs <- ceiling((abs(setsList[[i]][j] - setsList[[i]][j-1] + abs(negStart))/10))
        print(paste0('between column spaces = ',as.character(bcs)))

        posLine[[j]] <- paste(
          paste(rep(" ",times=interDist+bcs),collapse=""),
          setsList[[i]][j],collapse="")

        vertLine[[j]] <- paste(
          paste(rep(" ",times=interDist+bcs+2),collapse=""),
          '|',collapse="")

      }
    }
    #print(paste(posLine,collapse = ""))
    writeLines(paste(posLine,collapse = ""), fileConn)
    writeLines(paste(vertLine,collapse = ""), fileConn)
    close(fileConn)

    write.table(
      cbind(dfMat1,
            dfMat2[,c((i*10-9):(ifelse(i*10>nCol,nCol,i*10)))]),
      file=outFile,
      append=TRUE,
      col.names = FALSE,
      row.names = FALSE,
      quote = FALSE,
      sep=" ")

    fileConn <- file(outFile, open='a')
    writeLines('',fileConn)
    close(fileConn)
  }
}

#' make_prot_files
#'
#' Function to create the .alnprot.txt files from the
#' files from the '-protein.aln' files hosted on GitHub
#' @param filesIn list of full file paths
#' @param groupSize the size of the chunks of amino acid strings
#' the default is 10, and they are output as tab-separated columns
#' the output files(s)
#' @param allBases if set to FALSE, then only mismatched bases
#' will be included and the remainder will be set to '-' .
#' @param species human: hla , dog: dla, cow: bla, chicken: cla
#' @param pos the offset position to account for leader peptides for
#' the gene being anaylyzed.
#' @note This function is for general use
make_prot_files <- function(filesIn=NULL,
                           groupSize=10,
                           allBases=TRUE,
                           species='hla',
                           pos=29){
  lapply(filesIn,
         FUN = function(x){
           make_prot(x,
                     groupSize = groupSize,
                     allBases = allBases,
                     species = species,
                     pos = pos)
         })
}

#' make_p_group_file
#'
#' Function to create the _nom_p.txt files from the 'alnprot.txt'
#' files created with 'make_prot_file()'
#' @param filesInP These are the protein files created with
#' @param species human: hla , dog: dla, cow: bla, chicken: cla
#' make_prot_file(),
#' @note This function is for general use
make_p_group_file <- function(filesInP,
                              species='hla'){
  dfs <- lapply(filesInP,
         FUN=function(x){
           print(x)
           dfIn <- readLines(paste0("./data/",species,"/",x))
           comments  <- dfIn[which(grepl(pattern = "^#", dfIn))]
           protHead  <- dfIn[which(grepl(pattern = "^Prot", dfIn))]
           sequences <- dfIn[which(!grepl(pattern = "^Prot|^#", dfIn))]
           allSeqs   <- unlist(strsplit(sequences, split="\t"))
           pGroups   <- allSeqs[which(grepl(pattern="\\*", allSeqs))]
           pGroupsSplit  <- strsplit(pGroups, split="\\*|:")
           pGroupsSplit2 <- sapply(pGroupsSplit, "[[", 2)
           pGroupsSplit3 <- sapply(pGroupsSplit, "[[", 3)
           pGroupsSplit4 <- paste(pGroupsSplit2, pGroupsSplit3, sep=":")
           pIds <- pGroupsSplit4 %>% unique()
           #print(length(pIds))
           geneTmp <- gsub("(.*)\\*.*","\\1", pGroups[1], perl=TRUE)
           fileOutP <- paste0("./data/",species,"/",x,"_nom_p.txt")
           if(file.exists(fileOutP)){
             print(paste0(fileOutP, "aleady exists, removing ..."))
             file.remove(fileOutP)
           }
           print(paste0("Creating ", fileOutP, "..."))
           sink(fileOutP,append=TRUE)
           for(i in 1:length(pIds)){
            groupsTmp <-
              pGroups[which(grepl(pGroups,
                pattern = paste0(geneTmp,"\\*",pIds[i])))] %>%
                  unique() %>% gsub(".*\\*","",.)
            suffixGroup <- paste0(pIds[i],"P")
            outGroups  <- paste0(geneTmp,"*;",
                                 paste(groupsTmp, collapse="/"),
                                 ";",suffixGroup)
            writeLines(outGroups)
           }
           sink()
         })
}

#' make_prot_file
#'
#' Function to create the _nom_p.txt files from the 'alnprot.txt'
#' files created with 'make_prot_file()'
#' @param x x param
#' @param groupSize groupSize param
#' @param allBases allBases param
#' @param species species param
#' @param pos pos param
#' make_prot_file(),
#' @note This function is for general use
make_prot <- function(x = fileOutNamesProt[1],
                      groupSize = 10,
                      allBases = FALSE,
                      species = "dla",
                      pos=29){
  fullPath <- paste0("./data/",species,"/",x)
  print(paste0("Processing file ",fullPath))
  df <- read.table(fullPath,
                   sep="",
                   header = FALSE,
                   stringsAsFactors = FALSE)
  groupNames   <- df[,1]
  strings      <- df[,2]
  protLengths  <- df[,3]
  uniqProtLengths <- as.numeric(unique(protLengths))
  groupNamesFinal <- c()
  print(colnames(df))
  casted <- aggregate(df, V2 ~ V1, FUN=paste, collapse = "")
  #print(casted)
  #stop()
  dfMat <- data.frame()
  #for(i in 1:length(uniqProtLengths)){
  #dfNow <- df[which(df[,3]==uniqProtLengths[i]),]
  if(allBases==FALSE){
    for(j in 1:nrow(casted)){
      if(j!=1){
        #print(dfNow[1,2])
        #print(dfNow[j,2])
        matchInd <- mapply(function(x, y) which(x == y),
                           strsplit(casted[1,2], split=""),
                           strsplit(casted[j,2], split=""))
        if(length(matchInd[[1]])!=0){
          for(k in seq_along(matchInd)){
            substring(casted[j,2],
                      matchInd[k],
                      matchInd[k]) <- "-"
          }
        }
      }
    }
  }
  #account for leader peptide here:
  finalLength <- uniqProtLengths + pos
  repDots <- paste(rep(".",times=pos),collapse="")
  pepNow <- paste0(repDots,casted[,2])
  print(pepNow)

  choppedStrings <- array(rep("NA",times=10),
                          c(nrow=nrow(casted),
                            ncol=ceiling((max(uniqProtLengths)+pos)/groupSize),
                            1)
  )

  pepValues <- unlist(
    strsplit(
      gsub(paste0("(.{",groupSize,"})"), "\\1 ",
           pepNow,
           perl = TRUE),
      split=" "),
  )

  nCol<- ceiling((max(uniqProtLengths)+pos)/groupSize)
  nRow <- nrow(casted)

  for(i in 1:length(pepValues)){
    choppedStrings[ceiling(i/nCol),ifelse(i%%nCol!=0,i%%nCol,nCol),] <- pepValues[i]
    print(i)
  }

  dfMat <- choppedStrings

  lengths <- c("Prot",as.integer(-pos),1,unique(df[,3]))
  outFile <- paste0("./data/",species,"/",x,"prot.txt")
  fileConn <- file(outFile)
  writeLines(paste(lengths,collapse = "\t"), fileConn)
  close(fileConn)
  write.table(dfMat,
              file=outFile,
              append=TRUE,
              col.names = FALSE,
              row.names = FALSE,
              quote = FALSE,
              sep="\t")
  #stop()
}
