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

# fasta_to_aln
#'
#' Function to create a .aln file from a .fasta e.g. using
#' https://github.com/ANHIG/IMGTHLA/blob/Latest/fasta/DRB1_prot.fasta
#' as input.
#' @param fileIn Input file paths
#' @param appendName Boolean, append the 1st and 2nd allele names
#' Otherwise it will just use the 2nd column as the output name
#' @param species human: hla , dog: dla, cow: bla, chicken: cla
#' @note This function is for general use
#' @examples
#' \dontrun{
#' # The following is an example:
#'}
fasta_to_aln <- function(fileIn='',
                         appendName=FALSE,
                         species='hla'){
  fullPath <- paste0("./data/",species,"/",fileIn)
  print(paste0("Converting .fasta file ",
               fullPath,
               ' to .aln file'))
  conn = file(fullPath, "r")
  fileOut <- gsub(".fasta$",".aln",fileIn,perl=TRUE)
  fullOut <- paste0("./data/",species,"/",fileOut)
  if(exists(fullOut)){
    unlink(fullOut)
  }
  connOut = file(fullOut, 'a')
  #Check number of lines in file (doesnt include final blank line):
  numLines <- sapply(fullPath,R.utils::countLines)
  head1Check <- TRUE
  seqList <- list()
  headerLine <- TRUE
  i <- 1
  while(TRUE){
    lineIn = readLines(conn, n = 1)
    if(length(lineIn) == 0){
      break
    }else if(i==numLines[[1]]){
      lineOut <- paste(
        headerLine[2],
        paste(seqList, collapse=''),
        headerLine[3],
        collapse = "\t")
      writeLines(lineOut,
                 con = connOut)
      seqList <- list()
      #break
    }else if(grepl("^>", lineIn, perl=TRUE) & i==1){
      headerLine <- unlist(strsplit(lineIn,split=' +'))[2:3]
        #next
    }else if(grepl("^>", lineIn, perl=TRUE) & i!=1){
      if(head1Check <- TRUE){
        lineOut <- paste(
          headerLine[1],
          paste(seqList, collapse=''),
          headerLine[2],
          collapse = "\t")
        writeLines(lineOut, con=connOut)
        headerLine <- unlist(strsplit(lineIn,split=' +'))[2:3]
        seqList <- list()
        head1Check <- FALSE
        #next
      }else{
        lineOut <- paste(
          headerLine[2],
          paste(seqList, collapse=''),
          headerLine[3],
          collapse = "\t")
        headerLine <- unlist(strsplit(lineIn,split=' +'))[2:3]
        seqList <- list()
      #next
      }
    }else if(grepl("^[[:alnum:]]", lineIn)==TRUE){
      seqList <- rlist::list.append(seqList,lineIn)
      #next
    }#else{
    #  stop(paste0('Error: there is an indexing error
    #              outputing the .aln file'))
    #}
    print(lineIn)
    i <- i+1
  }
  close(conn)
}

# fasta_to_aln_2
#'
#' Function to create a .aln file from a .fasta e.g. using
#' https://github.com/ANHIG/IMGTHLA/blob/Latest/fasta/DRB1_prot.fasta
#' as input.
#' @param fileIn Input file paths
#' @param appendName Boolean, append the 1st and 2nd allele names
#' Otherwise it will just use the 2nd column as the output name
#' @param species human: hla , dog: dla, cow: bla, chicken: cla
#' @note This function is for general use
#' @examples
#' \dontrun{
#' # The following is an example:
#'}
fasta_to_aln_2 <- function(fileIn='',
                         appendName=FALSE,
                         species='hla'){
  fullPath <- paste0("./data/",species,"/",fileIn)
  print(paste0("Converting .fasta file ",
               fullPath,
               ' to .aln file'))
  conn = file(fullPath, "r")
  fileOut <- gsub(".fasta$",".aln",fileIn,perl=TRUE)
  fullOut <- paste0("./data/",species,"/",fileOut)
  if(exists(fullOut)){
    unlink(fullOut)
  }
  connOut = file(fullOut, 'a')
  #Check number of lines in file (doesnt include final blank line):
  numLines <- sapply(fullPath,R.utils::countLines)
  head1Check <- TRUE
  seqList <- list()
  headerLine <- TRUE
  i <- 1
  while(TRUE){
    lineIn = readLines(conn, n = 1)
    if(length(lineIn) == 0){
      break
    }else if(grepl("^>", lineIn, perl=TRUE) & i==1){
      headerLineFirst <-
        unlist(strsplit(lineIn,split=' +'))[2:3]
      #next
    }else if(grepl("^>", lineIn, perl=TRUE)){
        headerLine <- paste(
          headerLineFirst[1],
          paste(seqList[[1]], collapse=''),
          headerLineFirst[2],
          collapse = "\t")
        headerLine <-
          unlist(strsplit(lineIn,split=' +'))[2:3]
        lineOut2 <- paste(
          headerLine[1],
          paste(seqList[2:length(seqList)], collapse=''),
          headerLine[2],
          collapse = "\t")
        writeLines(lineOut1, con=connOut)
        writeLines(lineOut2, con=connOut)
        head1Check <- FALSE
        seqList <- list()
        writeLines(lineOut2, con=connOut)
    }else{
      seqList <- rlist::list.append(seqList,lineIn)
      #next
    }#else{
    #  stop(paste0('Error: there is an indexing error
    #              outputing the .aln file'))
    #}
    print(lineIn)
    i <- i+1
  }
  close(conn)
}

#' calc_spaces
#'
#' This function is for positioning the pos markers and '|' indicators
#' in the _nom_p.txt files
#' @param fileIn The full input file path
#' @param groupSize The number of AA's in each column
#' @param allBases Show all bases or just mutations
#' @param species Specify the species; e.g., human, dog, chicken, cow
#' @param pos The offset position (length) of the leader peptide
#' @param lengths The length of the AA sequences corresponding to the
#' third column in the download -protein align files.
#' @param ncharName the length in characters of the allele nanmes
#' corresponding to the first field in the downloaded -protein.aln files
#' @param protHeader First header lines of the output file fiven as a single
#' string with '\n' to denote newlines.
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

make_prot_file <- function(fileIn = "",
                        groupSize = 10,
                        allBases = FALSE,
                        species = "dla",
                        pos=23,
                        lengths,
                        nCharName,
                        protHeader=""){

  #library(data.table)
  #library(dplyr)
  #library(magrittr)
  #library(plyr)
  #library(abind)

  #fileOut <- "C:\\Users\\mmari\\OneDrive\\Desktop\\testOut.txt"
  #if(file.exists(fileOut)){
  #  unlink(fileOut)
  #}

  df <- read.table(fileIn,
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
  repLeadSpace <- paste(rep(" ", times=leaderPad), collapse="")
  pepNow <- paste0(repLeadSpace, casted[,2])
  #print(pepNow)

  choppedStrings <- array(rep("NA", times=10),
                          c(nrow=nrow(casted),
                            ncol=ceiling((max(uniqProtLengths)+leaderPad)/groupSize),
                            1)
  )

  pepVals <-
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
  #pepValues <- paste0(pepValues,rep("",nCol*nRow-length(pepValues)))
  pepValues <- paste0(pepVals,rep("",nCol*nRow))
  #pepValues <- c(
  #  paste(rep(' ',
  #            times=leaderPad,
  #            collapse = "")),
  #  pepVals2)

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
  #lengths = posVec

  fileOut <-
    paste0("./inst/extdata/",
           species,"/",
           paste0(gsub(".*/|-protein.aln$","",fileIn),
                  sep="_prot.txt"))

  #Check if file exists, if it does erase it:
  if(file.exists(fileOut)){
    unlink(fileOut)
  }

  dfMat1 <- dfMat[,1,1]
  dfMat2 <- dfMat[,-1,1]

  #Create output file and add header:
  fileConn <- file(fileOut)
  writeLines(protHeader, fileConn)
  close(fileConn)

  setsList <- list()
  negStart <- lengths[1]
  #for(i in seq_along(1:ceiling(nCol/10))){
  iterNums <- seq_along(1:ceiling(max(lengths)/100))
  for(i in iterNums){
    fileConn <- file(fileOut, open="a")
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
          paste(
                paste(rep(" ",
                          times=nCharName-nchar('prot')+4),
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
        browser()
        #prevLineRemaining <- ((i-1)*100 - ((max(unlist(setsList[[i-1]])) + abs(negStart))) + ifelse(j==1,2,0))
        prevLineRemaining <- ((i-1)*100 - ((max(unlist(setsList[[i-1]])) + abs(negStart))))
        #Subtract the inter column spaces from the previous line
        #prevLineRemaining <-
        #  prevLineRemaining +
        #  ifelse(prevLineRemaining>10,ceiling(prevLineRemaining/10)-1,0)
        #Calculate inter sequence length distance:
        interDist <-
          setsList[[i]][1] -
          (max(unlist(setsList[[i-1]])))
        finalDist <-
          interDist -
          prevLineRemaining +
          ((nCharName+1) - nchar('prot')) -
          2
          ##-3
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
        #Calculate the between column spaces because each column is grouped
        #into 10 AA currently
        #bcs <- ceiling((((setsList[[i]][1] -
        #                    lengths[length(unlist(setsList[1:i-1]))]) -
        #                   ((i-1)*100 - (lengths[length(unlist(setsList[1:i-1]))] + abs(negStart)))
        #)/10))
        bcs <- ifelse(abs(finalDist>10),ceiling(abs(finalDist)/10)-1,0)
        print(paste0('between column spaces = ',as.character(bcs)))

        posLine[[j]] <- paste('prot',
                              paste(rep(" ",times=finalDist+bcs),collapse=""),
                              setsList[[i]][1],collapse="")

        vertLine[[j]] <- paste(
                               paste(rep(" ",times=finalDist+bcs+5),collapse=""),
                               '|',collapse="")

      }else{
        if(i==max(iterNums)){
          print(paste0('Assigning final index position at ',
                       as.character(max(lengths))))
          interDist <-
            max(lengths) -
            setsList[[i]][j-1] -
            nchar(setsList[[i]][j-1]) -
            1
        }else{
          interDist <-
            setsList[[i]][j] -
            setsList[[i]][j-1] -
            nchar(setsList[[i]][j-1]) -
            1
        }

          print(paste0('interDist = ',as.character(interDist)))
          #Calculate the between column spaces becasue each column is grouped
          #into 10 AA currently
          #bcs <- ceiling((abs(setsList[[i]][j] - setsList[[i]][j-1] + abs(negStart))/10))
          bcs <- ifelse(i==max(iterNums),ceiling(abs(interDist)/10)-2,ceiling(abs(interDist)/10))
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
      file=fileOut,
      append=TRUE,
      col.names = FALSE,
      row.names = FALSE,
      quote = FALSE,
      sep=" ")

    fileConn <- file(fileOut, open='a')
    writeLines('',fileConn)
    close(fileConn)
  }
}

#' make_p_group_file
#'
#' Function to create the _nom_p.txt files from the 'alnprot.txt'
#' files created with 'make_prot_file()'
#' @param fileInP These are the protein files created with
#' @param species human: hla , dog: dla, cow: bla, chicken: cla
#' make_prot_file(),
#' @note This function is for general use
make_p_group_file <- function(fileInP,
                              species='hla'){
           print(fileInP)
           dfIn <- readLines(system.file(
             paste0("extdata/",species,"/",fileInP),
             package = "BIGDAWGv2"))
           comments  <- dfIn[which(grepl(pattern = "^#", dfIn))]
           protHead  <- dfIn[which(grepl(pattern = "^Prot", dfIn))]
           sequences <- dfIn[which(!grepl(pattern = "^Prot|^#", dfIn))]
           allSeqs   <- unlist(strsplit(sequences, split=" "))
           pGroups   <- allSeqs[which(grepl(pattern="\\*", allSeqs))]
           pGroupsSplit  <- strsplit(pGroups, split="\\*|:")
           pGroupsSplit2 <- sapply(pGroupsSplit, "[[", 2)
           pGroupsSplit3 <- sapply(pGroupsSplit, "[[", 3)
           pGroupsSplit4 <- paste(pGroupsSplit2, pGroupsSplit3, sep=":")
           pIds <- pGroupsSplit4 %>% unique()

           #pIds <- gsub(" .*","",pIds)

           #print(length(pIds))
           geneTmp <- gsub("(.*)\\*.*","\\1", pGroups[1], perl=TRUE)
           fileOutP <- paste0("./inst/extdata/",species,"/",fileInP,"_nom_p.txt")
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
}
