#' get_files
#'
#' functions to download files, used to create alnprot.txt _nom_p.txt files
#' @param fileUrl the Url of the file, not including the file name
#' @param fileDlNames the file name
#' @param fileOutNames the desired names for the downloaded files
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
                      fileOutNames){
  #downloads files such as *_prot.txt alignment files,
  #_nom_p.txt files, .fasta files, etc. for a select species.
  setwd("./data")
  lapply(fileDlNames,
        FUN = function(x){
              download.file(
              url = paste0(fileUrl, x),
              destfile=x,
              method="libcurl")
          }
        )
}

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
#' @note This function is for general use
make_prot_file <- function(filesIn=NULL,
                           groupSize=10,
                           allBases=TRUE,
                           species='hla'){
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
           dfMat <- data.frame()
           for(i in 1:length(uniqProtLengths)){
              dfNow <- df[which(df[,3]==uniqProtLengths[i]),]
              if(allBases==FALSE){
                for(j in 1:nrow(dfNow)){
                  if(j!=1){
                    #print(dfNow[1,2])
                    #print(dfNow[j,2])
                  matchInd <- mapply(function(x, y) which(x == y),
                                     strsplit(dfNow[1,2], split=""),
                                     strsplit(dfNow[j,2], split=""))
                    if(length(matchInd[[1]])!=0){
                      for(k in seq_along(matchInd)){
                        substring(dfNow[j,2],
                                  matchInd[k],
                                  matchInd[k]) <- "-"
                      }
                    }
                  }
                }
              }
              choppedStrings <- matrix(
                                  unlist(
                                    strsplit(
                                      gsub(paste0("(.{",groupSize,"})"), "\\1 ", dfNow[,2], perl = FALSE),
                                      split=" "
                                    ),
                                  ),
                                    ncol=ifelse(i==1,
                                                ceiling(uniqProtLengths[1]/groupSize),
                                                ceiling(length((uniqProtLengths[i-1]+1):uniqProtLengths[i])/groupSize)
                                                ),
                                    byrow = TRUE
                                  )
             groupNamesFinal <- dfNow[,1]
             dfCol <- cbind(groupNamesFinal,
                            choppedStrings)
             dfMat <- plyr::rbind.fill(dfMat,data.frame(dfCol))
            }
           lengths <- c("Prot",1,unique(df[,3]))
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
           dfIn <- readLines(x)
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
