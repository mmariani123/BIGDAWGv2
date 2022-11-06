#!/usr/bin/env/Rscript

# Script to create dla_nom_p.txt for DLA

library(data.table)
library(dplyr)
library(magrittr)

dla.prot.files <- list.files(path=".", pattern="canid")

dla.prots.list <- lapply(dla.prot.files,
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
