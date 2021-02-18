# Author Neeraj kumar
# updated Feb 2021

# Description -------------------------------------------------------------

# this script makes a binary vector for genome/metagenome profile of Pfams.
#   It required the Pfam.A.clans file from Pfam data base
#  (ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/). Different version of Pfam vary in numbers. 
#   This script generates an output  .
#   PfamVector.txt file is used as an input for the next step.

#setwd("PATH") 
  getwd()

# Usage -------------------------------------------------------------------

# Rscript --vanilla script_2_pfam_vector.R -i pfam.A.clans -p pathForParsedFile -o outputfile


# install required package ------------------------------------------------

  paks <- c("tidyverse","data.table","plyr","optparse") # give required package names here
  paks<-lapply(paks, function(x) suppressMessages(require(x, character.only = TRUE)))
  rm(paks)
  # Mkaing opotions ---------------------------------------------------------
  #library("optparse")
  # maiking option using package "optparse"
  option_list = list(
    make_option(c("-i", "--Pfamfile"), type="character", default=NULL, 
                help="Pfam.A.clans file name", metavar="character"),
    make_option(c("-p", "--path"), type="character", default=NULL, 
                help="path for hmmscan pfam parsed file ", metavar="character"),
    make_option(c("-o", "--outputFileName"), type="character", default="PfamVector.txt", 
                help="output file name [default= %default]", metavar="character")) 
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
  # if (is.null(opt$file)){
  #   print_help(opt_parser)
  #   stop("At least one argument must be supplied (input file).n", call.=FALSE)
  # }
# read input files --------------------------------------------------------------
#args = commandArgs(trailingOnly=TRUE)

  pfam <- read.delim(opt$Pfamfile,
                              header=FALSE,
                              stringsAsFactors=FALSE)
  # pfam <- read.delim("Pfam-A.clans_30_8_18.tsv",
  #                            header=FALSE,
  #                            stringsAsFactors=FALSE)
  colnames(pfam) <- c("Protein_family","class","name","target.name","Description")
  Genome_data <- list.files(opt$path,full.names = TRUE)

# main part ------------------------------------------------------------

  PfamVector_2 <- NULL
  #pfamAbundance_2 <- NULL
  count <- 0
  for(filename in Genome_data)
  {
    sample <- read.csv(filename,header=FALSE, stringsAsFactors=FALSE)
    #print (sample)
    sbt = strsplit(as.character(sample[,1]),split="\\.")
    temp <-unlist(strsplit(as.character(filename),split="\\/")) # made unlist to get file name from last vector
    #sample_name <-strsplit(as.character(temp[[1]][8]),split="\\.")[[1]][1]
    print(filename)
    sample_name <-temp[length(temp)] # length(temp) return sample name.
    print(sample_name)
    #print(sbt)
    ncol = max(sapply(sbt,length))
    tempsplit <- as.data.table(lapply(1:ncol,function(i)sapply(sbt,"[",i)))
    #print(tempsplit)
    sample[,"PfamID"] <- tempsplit$V1
    #print(sample[,"PfamID"])
  sample$gid <- rep(sample_name,times = nrow(sample))
  g <- unique(sample$gid)
  p <- unique(pfam$Protein_family)
  gp <- unique(sample$PfamID)
    
 #******************** maiking Biinary Vector ******************
    
    PfamVector <- matrix(data=0,nrow=length(p),ncol=1,dimnames = list(p,g))
    imatch <- match(gp,rownames(PfamVector))
    PfamVector[imatch,] <- 1
    PfamVector_2<-cbind(PfamVector_2,PfamVector)
  }
  PfamVector_2 <- as.data.frame(PfamVector_2)
  #write table for pfam binary vector
  PfamVector_2$PfamID <- row.names(PfamVector_2)

# output table ------------------------------------------------------------
  write.table(PfamVector_2,file=opt$outputFileName,row.names= FALSE,col.names = TRUE,sep = "\t",quote =FALSE)
  
# Abundance calculations --------------------------------------------------


  #   # ****************counting pfam  abundance ############
  #   testsample <- ddply(sample,.(PfamID,gid),nrow) # main function to calculate abundance
  #   colnames(testsample)[3] <- sample_name
  #   colnames(testsample)[1] <- "Protein_family"
  #  
  #   testsample_pfam <-left_join(pfam,testsample,by="Protein_family")
  #  
  #   pfamAbundance <- as.data.frame(testsample_pfam[,c("Protein_family",sample_name)])
  #   pfamAbundance[is.na(pfamAbundance)] <- 0
  #   if(count==0){
  #     pfamAbundance_2<-pfamAbundance
  #     }else
  #     {pfamAbundance_2<- cbind(pfamAbundance_2,pfamAbundance[sample_name])
  #       }
  #   count <- count+1
  #   }
  # 
  #   
  #   PfamVector_2$description <- pfam[match(rownames(PfamVector_2),pfam$Protein_family),"Description"]
  # 
  # #write table for vector including Pfam which does not have  
  #   write.table(PfamVector_2,file="metavector.txt",row.names= TRUE,col.names = TRUE,sep = "\t",quote =FALSE)
  #   write.table(pfamAbundance_2,file = "metaPfamAbundance_2.txt",row.names= TRUE,col.names = TRUE,sep = "\t",quote =FALSE)
  # 
  # # write table excluding Pfam which does not have any match in database
  #   PfamVector_2$description <- NULL
  #   PfamVector_2 <- PfamVector_2[rowSums(PfamVector_2)!=0,]
  #   write.table(bla_2,file="metavector_ExcludeRow_zero.txt",row.names= TRUE,col.names = TRUE,sep = "\t",quote = FALSE)
  #   rownames(pfamAbundance_2) <- pfamAbundance_2$Protein_family
  #   pfamAbundance_2$Protein_family <- NULL 
  #   pfamAbundance_2 <-  pfamAbundance_2[rowSums(pfamAbundance_2)!=0,]
  #   pfamAbundance_2$Protein_family <- rownames(pfamAbundance_2) 
  #   write.table(pfamAbundance_2,file="PfamAbundance_ExcludeRow_zero.txt",row.names= TRUE,col.names = TRUE,sep = "\t",quote =FALSE)
  # 
# write the tables --------------------------------------------------------

  
# #######################  Scripts end here ############
# 
