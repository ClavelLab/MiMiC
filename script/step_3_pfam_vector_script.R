#Neeraj kumar
#pfam vector 
#last updated 16 MARCH 2020
#use-> Rscript script_2_pfam_vector.R Pfam.A.clans Path_to_parsed_files # on terminal
#metavector.txt is the output file for binary vector.

getwd()
#setwd("PATH") 
  require(c("data.table","dplyr","plyr"))
  library(data.table)
  library(dplyr)
  library(plyr)
  args = commandArgs(trailingOnly=TRUE)

  Pfam.A.clans <- read.delim(args[1], header=FALSE, stringsAsFactors=FALSE)
  pfam <- Pfam.A.clans
  Genome_data <- list.files(path=args[2],full.names = TRUE)
  colnames(pfam) <- c("Protein_family","class","name","target.name","Description")

  bla_2 <- NULL
  pfamAbundance_2 <- NULL
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
      sample[,"pfamid"] <- tempsplit$V1
  #print(sample[,"pfamid"])
      sample$gid <- rep(sample_name,times = nrow(sample))
      g <- unique(sample$gid)
      p <- unique(pfam$Protein_family)
      gp <- unique(sample$pfamid)
  
  ######### Vector #############
  #rm(Genome_to_make_vector,Pfam.A.clans)
  bla <- matrix(data=0,nrow=length(p),ncol=1,dimnames = list(p,g))
  imatch <- match(gp,rownames(bla))
  bla[imatch,] <- 1
  bla_2<-cbind(bla_2,bla)
  

  ###### abundance ############
  testsample <- ddply(sample,.(pfamid,gid),nrow) # main function to calculate abundance
  colnames(testsample)[3] <- sample_name
  colnames(testsample)[1] <- "Protein_family"
  #if(all(testsample$Protein_family!=pfam$Protein_family)) next
  testsample_pfam <-left_join(pfam,testsample,by="Protein_family")
  #testsample_pfam_temp <-left_join(pfam,testsample,by="Protein_family")
  pfamAbundance <- as.data.frame(testsample_pfam[,c("Protein_family",sample_name)])
  pfamAbundance[is.na(pfamAbundance)] <- 0
  if(count==0){
    pfamAbundance_2<-pfamAbundance
    }else
    {pfamAbundance_2<- cbind(pfamAbundance_2,pfamAbundance[sample_name])
      }
  count <- count+1
}

bla_2 <- as.data.frame(bla_2)
bla_2$description <- pfam[match(rownames(bla_2),pfam$Protein_family),"Description"]

#write table for vector including Pfam which does not have  
write.table(bla_2,file="metavector.txt",row.names= TRUE,col.names = TRUE,sep = "\t",quote =FALSE)
write.table(pfamAbundance_2,file = "metaPfamAbundance_2.txt",row.names= TRUE,col.names = TRUE,sep = "\t",quote =FALSE)

# write table excluding Pfam which does not have any match in database
bla_2$description <- NULL
bla_2 <- bla_2[rowSums(bla_2)!=0,]
write.table(bla_2,file="metavector_ExcludeRow_zero.txt",row.names= TRUE,col.names = TRUE,sep = "\t",quote = FALSE)
rownames(pfamAbundance_2) <- pfamAbundance_2$Protein_family
pfamAbundance_2$Protein_family <- NULL 
pfamAbundance_2 <-  pfamAbundance_2[rowSums(pfamAbundance_2)!=0,]
pfamAbundance_2$Protein_family <- rownames(pfamAbundance_2) 
write.table(pfamAbundance_2,file="PfamAbundance_ExcludeRow_zero.txt",row.names= TRUE,col.names = TRUE,sep = "\t",quote =FALSE)

#######################  Scripts and here ############

