#!/usr/bin/env Rscript

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
#' dla.prot.files <- list.files(path = "./data/dla",
#'                             pattern = "dla",
#'                             full.names = TRUE)
#'calc_spaces(lengths = posVec,
#'            nCharName = 13)
#'}

calc_spaces <- function(x = fileOutNamesProt[3],
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
  leaderPad <- plyr::round_any(pos, 10, f=ceiling) - pos
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
  pepValues <- paste0(pepValues,rep("",nCol*nRow-length(pepValues)))

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
  lengths <- c(as.integer(-pos)-leaderPad, 1, unique(df[,3]))

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
  for(i in seq_along(1:ceiling(sum(abs(range(lengths)))/100))){
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
