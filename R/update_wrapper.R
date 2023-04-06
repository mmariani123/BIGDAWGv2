#' Update function for protein aligment upon new IMGT HLA data release
#'
#' This updates the protein aligment used in checking HLA loci
#' and alleles as well as in the amino acid analysis.
#' @param Restore Logical specifying if the original alignment file be restored.
#' @param Force Logical specifiying if update should be forced.
#' @param Output Logical indicating if error reporting should be written to file.
#' @param CreateNew Logical indicating if error reporting should be written to file
#' @param Species Logical indicating if error reporting should be written to file
#' @param OutputDir Where the new/updated object will be placed.
UpdateRelease <- function(Force=F,
                          Restore=F,
                          Output=F,
                          CreateNew=F,
                          Species='hla',
                          OutputDir=getwd()){

##browser()

if(species=='dla' & CreateNew==TRUE){

  #Step 1 get loci

  Loci <- c("A",
            "B",
            "C",
            "DPA1",
            "DPB1",
            "DQA1",
            "DQB1",
            "DRB1",
            "DRB3",
            "DRB4",
            "DRB5")

  loci <- c('12',
            '64',
            '79',
            '88',
            'DQA1',
            'DQB1',
            'DRB1')

  Loci.get <- c("A",
                "B",
                "C",
                "DPA1",
                "DPB1",
                "DQA1",
                "DQB1",
                "DRB")

  loci.get <- c('12',
                '64',
                '79',
                '88',
                'DQA1',
                'DQB1',
                'DRB1')

  #Exon Info
  RefTab <- BIGDAWG::ExonPtnList$RefExons

  RefTab <- data.frame(Locus=c('12',
                               '64',
                               '79',
                               '88',
                               'DQA1',
                               'DQB1',
                               'DRB1'),
                       Reference.Locus=c('12',
                                         '64',
                                         '79',
                                         '88',
                                         'DQA1',
                                         'DQB1',
                                         'DRB1'),
                       Reference.Allele=c('001:02',
                                          '001:02',
                                          '001:01',
                                          '001:01',
                                          '001:01',
                                          '001:01',
                                          '001:01'),
                       Reference.Accession=c('DLA00100',
                                             'DLA00200',
                                             'DLA00300',
                                             'DLA00400',
                                             'DLA00500',
                                             'DLA00600',
                                             'DLA00700'),
                       Reference.Start=c(-7,
                                         -7,
                                         -7,
                                         -7,
                                         -7,
                                         -7,
                                         -7),
                       stringsAsFactors = FALSE)

  #Release$V1
  #Release[[1]][1]
  #Release[[1]][2]

  Release <- data.frame(
               V1=BIGDAWG::ExonPtnList$Release.Version,
               V2=BIGDAWG::ExonPtnList$Release.Date,
               stringsAsFactors = FALSE)

  Release[1,] <- "2023-03-29"
  Release[2,] <- "IPD-IMGT/HLA 3.51.0"

  cat("Formatting alignment files.\n")
  for(i in 1:length(loci)){
    locus <- loci[i]
    ExonPtnAlign.Create(locus,RefTab,Species)
  }
  AlignObj.Update(loci,Release,RefTab,Species)

  cat("Created.\n")

  #STEP 2: Download protein alignments and other ancillary files
  #cat("Creating reference object for the amino acid analysis.\n")
  #cat("Downloading alignment files from the IMGT/HLA.\n")
  #GetFiles(Loci.get)
  #Release <-
  #  read.table('Release.txt',
  #             sep="\t") # created during GetFiles download

  #STEP 3: Format alignments for exons of interest
  #cat("Formatting alignment files.\n")
  #for(i in 1:length(Loci)){
  #  Locus <- Loci[i]
  #  ExonPtnAlign.Create(Locus,RefTab)
  #}

  #STEP 4: Create ExonPtnAlign list object for BIGDAWG package
  #AlignObj.Update(Loci,Release,RefTab)

  #STEP 5: Clean up
  #cat("Cleaning up.\n")
  #invisible(file.remove(dir()[which(dir() %in% Safe!=T)]))

  #cat("Created.\n")

  #switchOutput <- switch(
  #  Species,
  #  'hla' = create_human(),
  #  'dla' = create_dog(),
  #  'cla' = create_cow(),
  #  'gla' = create_chicken(),
  #)

}else{

  if(!inherits(
    try(XML::readHTMLTable(
          "http://cran.r-project.org/web/packages/BIGDAWG/index.html",
          header=F),
        silent=T),
        "try-error")){

    MainDir <- getwd()
    on.exit(setwd(MainDir), add = TRUE)

    getDir <- path.package('BIGDAWG')
    putDir <- paste(getDir,"/data",sep="")
    if(!dir.exists(putDir)){dir.create(putDir)}

    if(!Restore) {

      #Check current version against BIGDAWG version
      if(!Force) {

        setwd(putDir)

        # Get IMGT Release Version
        invisible(
          download.file(
            paste0('ftp://ftp.ebi.ac.uk/pub/databases/ipd',
                   '/imgt/hla/release_version.txt'),
                   destfile="release_version.txt",
                   method="libcurl"))
        Release <- read.table("release_version.txt",
                              comment.char="",
                              sep="\t")
        Release <- apply(Release,
                         MARGIN=1,
                         FUN=function(x) gsub(": ",":",x)
                         )
        RV.current <- unlist(strsplit(Release[3],split=":"))[2]
        file.remove("release_version.txt")

        # Get BIGDAWG
        UPL <- paste(path.package('BIGDAWG'),
                     "/data/UpdatePtnAlign.RData",
                     sep="")
        UpdatePtnList <- NULL
        rm(UpdatePtnList)
        if(file.exists(UPL)){
          load(UPL)
          EPL <- UpdatePtnList
          rm(UpdatePtnList,UPL)
          UPL.flag=T
        }else{
          EPL <- ExonPtnList
          UPL.flag=F }

        RV.BIGDAWG <- EPL$Release.Version

        cat('Versions:\n',
            'IMGT/HLA current: ',
            RV.current,
            '\n BIGDAWG version: ',
            RV.BIGDAWG,'\n')
        if(grepl(RV.current,RV.BIGDAWG)){
          Flag <- T
        }else{
          Flag <- F
        }

      }else{

        Flag <- F

      }# End if() for setting Flag

      #Run Update if Flag = T
      if(Flag) {

        cat(paste0('\nYour database seems up to date. ',
                   'Use Force = T to force the update.'))

      }else{

        # For creating UpdatePtnAlign.RData object
        # Define download directory
        Safe <- OutputDir
        Safe <- c(Safe[!grepl(".txt",Safe)],
                  "UpdatePtnAlign.RData")

        #STEP 1: Define Loci and Read in Reference Exon Map Files
        Loci <- c("A",
                  "B",
                  "C",
                  "DPA1",
                  "DPB1",
                  "DQA1",
                  "DQB1",
                  "DRB1",
                  "DRB3",
                  "DRB4",
                  "DRB5")

        #Currently DRB1, DRB3, DRB4, DRB5 aligments in single file
        #Remove if split into individual files
        Loci.get <- c("A",
                      "B",
                      "C",
                      "DPA1",
                      "DPB1",
                      "DQA1",
                      "DQB1",
                      "DRB")

        #Exon Info
        RefTab <- BIGDAWG::ExonPtnList$RefExons

        #STEP 2: Download protein alignments and other ancillary files
        cat("Updating reference object for the amino acid analysis.\n")
        cat("Downloading alignment files from the IMGT/HLA.\n")
        GetFiles(Loci.get)
        Release <-
          read.table('Release.txt',
                     sep="\t") # created during GetFiles download

        #STEP 3: Format alignments for exons of interest
        cat("Formatting alignment files.\n")
        for(i in 1:length(Loci)){
          Locus <- Loci[i]
          ExonPtnAlign.Create(Locus,RefTab,Species = Species)
        }

        #browser()
        #STEP 4: Create ExonPtnAlign list object for BIGDAWG package
        AlignObj.Update(Loci,Release,RefTab)

        #STEP 5: Clean up
        cat("Cleaning up.\n")
        invisible(file.remove(dir()[which(dir() %in% Safe!=T)]))

        cat("Updated.\n")
      }

    }else if(Restore){

      setwd(putDir)
      if(!file.exists('UpdatePtnAlign.RData')){
        stop("No prior update to restore.",
             call.= F)
      }
      cat(paste0('Restoring original alignment reference object ',
                 'for amino acid analysis.\n'))
      invisible(file.remove('UpdatePtnAlign.RData'))
      cat('Restored.\n')

    }

  }else{

    Err.Log(Output,"No.Internet or page has moved")
    stop("Update stopped, No Internet or Page has moved.",
         call.=F)

  }

}

}
