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

#DRB1 starts at -29,
#DQA1 starts at -23,
#DQB1 starts at -32,
#DPA1 starts at -31
#and DPB1 starts at -29,
#so the protein alignments
#all start at -30 or -40

#' make_prot_file
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
make_prot_file <- function(filesIn=NULL,
                           groupSize=10,
                           allBases=TRUE,
                           species='hla',
                           pos=29){
  lapply(filesIn,
         FUN = function(x){
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

              choppedStrings <- array(NA,
                                      c(nrow=nrow(casted),
                                      ncol=ceiling((max(uniqProtLengths)+pos)/groupSize),
                                      groupSize)
              )

              pepValues <- unlist(
                strsplit(
                  gsub(paste0("(.{",groupSize,"})"), "\\1 ",
                       pepNow,
                       perl = TRUE),
                  split=" "),
              )

              nCol<- 1:ceiling((max(uniqProtLengths)+pos)/groupSize)
              nRow <- nrow(casted)

              for(i in 1:length(pepValues)){
                  choppedStrings[ceiling(i/13),i,] <- pepValues[i]
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
        })

}

#' make_prot_file
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
