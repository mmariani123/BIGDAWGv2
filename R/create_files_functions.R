#!/usr/bin/env/Rscript

#' BIGDAWGv2 Main Wrapper Function
#'
#' This is the main wrapper function for each analysis.
#' @param Data Name of the genotype data file.
#' @param HLA Logical Indicating whether data is HLA class I/II genotyping data only.
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
#'
# Script to create dla_nom_p.txt for DLA

library(data.table)
library(dplyr)
library(magrittr)

dla.prot.files <- list.files(path = ".",
                             pattern = "canid")

dla.prots.list <- lapply(
       dla.prot.files,
       read.table,
       sep="",
       header=FALSE,
       stringsAsFactors=FALSE,
       skip=9,
       strip.white=TRUE,
       colClasses="character")

dla.bind <- do.call(rbind, dla.prots.list)

dla.prots.protein.list <- lapply(dla.prot.files,
                          FUN=function(x){
                            lines.in=readLines(x)
                            lines.out=lines.in[which(grepl("Prot",lines.in))]
                            return(lines.out)
                          })

dog.loci <- dla.bind[,1] %>%
  gsub("\\*.*","",.) %>%
    unique()

dog.loci
#For dog:
"DLA-DQA1"
"DLA-DQB1"
"DLA-DRB1"

#for(i in 1:length(dog.loci)){
#  select.loci <- dla.bind[which(grepl(dog.loci[i],dla.bind[,1])),1]
#  #print(select.loci)
#  #collapse.loci <- paste(gsub(".*\\*","",select.loci),collapse="/")
#
#}
#select.loci
#collapse.loci

nom.p <- dla.bind[,1]
nom.p <- gsub("\\*",";",nom.p) %>% paste0(.,";")

write.table(nom.p,
            file="dla.nom.p.txt",
            col.names = FALSE,
            row.names = FALSE,
            sep = " ",
            quote=FALSE)

# file: dla_nom_p.txt
# date: 2022-11-06
# version: Copyright © EMBL 2022
# origin: github.com/m.mariani123/BIGDAWGv2/dla_nom_p.txt
# repository: github.com/m.mariani123/BIGDAWGv2/dla_nom_p.txt
# author: Michael P Mariani <m.mariani123@gmail.com>

dog.fasta.files <- list("dla-12-protein.fasta",
                  "dla-64-protein.fasta",
                  "dla-79-protein.fasta",
                  "dla-88-protein.fasta",
                  "dla-dqa1-protein.fasta",
                  "dla-dqb1-protein.fasta",
                  "dla-drb1-protein.fasta")

dog.aln.files <- list("dla-12-protein.aln",
                      "dla-64-protein.aln",
                      "dla-79-protein.aln",
                      "dla-88-protein.aln",
                      "dla-dqa1-protein.aln",
                      "dla-dqb1-protein.aln",
                      "dla-drb1-protein.aln")

fileUrl      <- "https://raw.githubusercontent.com/ANHIG/IPDMHC/Latest/alignments/"
fileDlNames  <- dog.aln.files

fileOutNames <- list("dla-12-protein.fasta",
                     "dla-64-protein.fasta",
                     "dla-79-protein.fasta",
                     "dla-88-protein.fasta",
                     "dla-dqa1-protein.fasta",
                     "dla-dqb1-protein.fasta",
                     "dla-drb1-protein.fasta")

fileOutNames <- list("dla-12-protein.aln",
                     "dla-64-protein.aln",
                     "dla-79-protein.aln",
                     "dla-88-protein.aln",
                     "dla-dqa1-protein.aln",
                     "dla-dqb1-protein.aln",
                     "dla-drb1-protein.aln")

get_files <- function(fileUrl,
                      fileDlNames,
                      fileOutNames){
  #downloads *_prot.txt alignment files
  #downloads hla_nom_p.txt file

  # Get P-Groups Files
  lapply(fileDlNames,
        FUN = function(x){
              download.file(
              url = paste0(fileUrl, x),
              destfile=x,
              method="libcurl")
          }
        )

}

get_files(fileUrl=fileUrl,
          fileDlNames=fileDlNames,
          fileOutNames=fileOutNames)

make_prot_file <- function(filesIn=NULL,
                           groupSize=10){

  lapply(filesIn,
         FUN = function(x){
           print(paste0("Processing file ",x))
           df <- read.table(x,
                      sep="",
                      header = FALSE,
                      stringsAsFactors = FALSE)

           groupNames   <- df[,1]
           strings      <- df[,2]
           protLengths  <- df[,3]
           uniqProtLengths <- as.numeric(unique(protLengths))

           groupNamesFinal <- c()
           #dfMat <- c()
           dfMat <- data.frame()
           for(i in 1:length(uniqProtLengths)){
             dfNow <- df[which(df[,3]==uniqProtLengths[i]),]
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

           outFile <- paste0(x,"prot.txt")
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

make_prot_file(fileOutNames,
               groupSize = 10)

pGroupHeader <- paste("# file: canid_DRB1_prot.txt",
                      "# date: 11/04/2022",
                      "# version: Copyright © EMBL 2022",
                      "# origin: https://www.ebi.ac.uk/ipd/mhc/alignment for canis adustis DRB1",
                      "# repository: github.com/mmariani123/BIGDAWGv2/canid_DRB1_prot.txt",
                      "# author: Michael P. Mariani (m.mariani123@gmail.com)", sep="\n")

pGroupFiles <- list("dla-12-protein.alnprot.txt",
                    "dla-64-protein.alnprot.txt",
                    "dla-79-protein.alnprot.txt",
                    "dla-88-protein.alnprot.txt",
                    "dla-dqa1-protein.alnprot.txt",
                    "dla-dqb1-protein.alnprot.txt",
                    "dla-drb1-protein.alnprot.txt")

make_p_group_file <- function(filesInP){
  dfs <- lapply(filesInP,
         FUN=function(x){
           dfIn <- readLines(x)
           #print(head(dfIn, n=10))
           comments  <- dfIn[which(grepl(pattern = "^#", dfIn))]
           protHead  <- dfIn[which(grepl(pattern = "^Prot", dfIn))]
           sequences <- dfIn[which(!grepl(pattern = "^Prot|^#", dfIn))]
           allSeqs   <- unlist(strsplit(sequences, split="\t"))
           pGroups   <- allSeqs[which(grepl(pattern="\\*", allSeqs))]
           #print(comments)
           #print(protHead)
           #print(length(pGroups))
           #print(pGroups)
           pGroupsSplit  <- strsplit(pGroups, split="\\*|:")
           pGroupsSplit2 <- sapply(pGroupsSplit, "[[", 2)
           pGroupsSplit3 <- sapply(pGroupsSplit, "[[", 3)
           pGroupsSplit4 <- paste(pGroupsSplit2, pGroupsSplit3, sep=":")
           pIds <- pGroupsSplit4 %>% unique()
           #print(length(pIds))
           geneTmp <- gsub("(.*)\\*.*","\\1", pGroups[1], perl=TRUE)
           #print(length(pIds))
           #print(paste0(x,"_nom_p.txt"))
           fileOutP <- paste0(x,"_nom_p.txt")
           if(file.exists(fileOutP)){
             print(paste0(fileOutP, "aleady exists, removing ..."))
             file.remove(fileOutP)
             #stop()
           }
           print(paste0("Creating ", fileOutP, "..."))
           sink(fileOutP,append=TRUE)
           for(i in 1:length(pIds)){
            #print(i)
            #print(pIds[i])
            #print(genesTmp[i])
            groupsTmp <-
              pGroups[which(grepl(pGroups,
                pattern = paste0(geneTmp,"\\*",pIds[i])))] %>%
                  unique() %>% gsub(".*\\*","",.)
            #groupsTmp <- gsub(".*\\*", "", groupsTmp, perl=TRUE)
            #print(groupsTmp)
            suffixGroup <- paste0(pIds[i],"P")
            outGroups  <- paste0(geneTmp,"*;",
                                 paste(groupsTmp, collapse="/"),
                                 ";",suffixGroup)
            #print(outGroups)
            writeLines(outGroups)
            #stop()
           }
           sink()
         })
         #stop()
         #return(dfs)
}

#DLA-DQA1*;005:01:1/005:01:2;005:01P

dfs <- make_p_group_file(filesInP=pGroupFiles)

dfs[1]

df <- readLines(con=fileOutName)
headers <- df[c(TRUE,FALSE)]
sequences <- df[c(FALSE,TRUE)]

RD <- unlist(strsplit(df[2],split=" "))[3]
RV <- paste(
  unlist(
    strsplit(df[3],
             split=" "))[3:4],
  collapse=" ")

# Get Locus Based Alignments
#for(i in 1:length(Loci)){
#  Locus <- Loci[i]
#  FileName <- paste0(Locus,"_prot.txt")
#  URL <- "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/alignments/"
#  download.file(URL,
#                destfile = FileName,
#                method="libcurl")
#}

#}
