
#neeraj kumar
#pfam coverage of metagenome by refernce genome
#minimum number of species covering maximum number of metagenome
#updated 15 March 2020

install.packages(c("tidyverse","plyr","data.table","ggplot2","inflection","taRifx","dplyr","readr"))
  
  library(tidyverse)
  library(plyr)
  library(data.table)
  library(ggplot2)
  library(inflection)
  library(taRifx)
  library(dplyr)
  library(readr)

#read the reference genome data set
	GenomeVector <- read.table("PiBC_Genome_Binary_Vector_111_update_3oct_2019.txt",
                             sep = "\t",header = T,stringsAsFactors = FALSE)
#read the metagenome file
	MetaGenome<- read.delim("PiBC_Vector_284_updated_3rdOct.txt",header = TRUE,
                          stringsAsFactors = FALSE,check.names = F)
	MetaGenome <-MetaGenome[,c(1:10)] #                                        ---> just for test
	myconsortia=111 # minimum consortia size
#source the function of knee point
	source("knee.R")
	rownames(MetaGenome) <- MetaGenome$PfamID # to make sure if rows name remain same 
	MetaGenome$PfamID <- NULL
#row.names(GenomeVector) <- GenomeVector$PfamID # skip if row names alreay assigned and PfamID is not available
	GenomeVector$PfamID <- NULL
#function to select the bacteria from reference set.
	scoreDiff_weighted <- function(genomes, metagenome,Minimal_consortia_size){
#row.names(GenomeVector) <- GenomeVector$PfamID
		g <- as.matrix(genomes)
		y <- as.matrix(metagenome)
		mscore <- sum(y)
		mscore1 <- sum(y)
		res <- NULL
    		#Minimal_consortia_size=23
    		for(i in 1:Minimal_consortia_size){
      			matchCount <- apply(g, 2, function(x) sum(y[x >= 1.0])) # ---> only count of pfam match
			totalpfamG <- apply(g,2,function(x) sum(x>=1)) # in case of weighted
      			misMatchCount <- totalpfamG - matchCount # number Of Pfams which are present in Bacteria but not in Metagenome.
      			matchratio <- ((matchCount)/(matchCount+misMatchCount)) # (x)/(x+y) # ---> mimic_B             
      			matchratio[matchratio=="NaN"]<- 0 # because of some values become 0, when 0/0 it gives NaN
      			print(max(matchratio))
      			print(paste ("number of itreation" ,i,sep="="))
      			a <- g[,names(sort(matchratio)[length(matchratio)])] # selected genome #for match/mismatch
      			ascore <- sort(matchratio)[length(matchratio)] # the score given to the genome while being selected.
      			pfamInGenomePerItration = unname(totalpfamG[names(sort(matchratio)[length(matchratio)])])
      			exclusiveMatch <- unname(matchCount[names(sort(matchratio)[length(matchratio)])])
      			exclusiveMismatchPerIteration <- unname(misMatchCount[names(sort(matchratio)[length(matchratio)])])
      
      			print(paste("pfam",pfamInGenomePerItration,"exclusiveMatch",exclusiveMatch,"exclusiveMismatchPerIteration",
                  		exclusiveMismatchPerIteration,sep = "="))
      #print the genome with the name
      #print (ascore)
      #unname(ascore)
      			b <- as.character(names(sort(matchratio)[length(matchratio)]))
      			print(b)
      			g <-g[,-which(colnames(g)==b),drop=F] # exclude selected genome.
      			if(i==1){
        			covered <- sum(y[a>=1])
        			rscore <- sum(y[a>=1])
        			percent <- (covered/mscore)*100
      				} else {
        			rscore <- sum(y[a>=1])
        			covered <- covered+rscore
        			percent <- (covered/mscore)*100
      				}
      			test <- names(y[a>=1,])
      			g[test,] <- 0 # making all pfam as "zero" in reference genome set ( those function were already taken ) 
      			y[a>=1]<-0 # makes all metagenome matching with highest scoring genome as '0'
      #if(rscore==0){break}
      			mscore1 <- sum(y)
      #res <- rbind(res,cbind.data.frame(bname=colnames(g)[i],percentage=percent,score=rscore)) # this line has to corrected in all scripts.
      			if(rscore==0){break}
      #if(rscore==0){next} # in case if there is 0 score then do not add to the list and go to next genome.
      			RelativeMatchPerIteration=(rscore / exclusiveMismatchPerIteration)
      
      			res <- rbind(res,cbind.data.frame(bname=names(sort(matchratio)[length(matchratio)]),
                                        percentage=percent,score=rscore,ScorePerIteration=ascore,
                                        pfamInGenomePerItration=pfamInGenomePerItration,
                                        exclusiveMismatchPerIteration=exclusiveMismatchPerIteration,
                                        RelativeMatchPerIteration=RelativeMatchPerIteration)) #for match/Mismatch
    
    						}
    return(res)
}
 
  	b_percentage_list_W <- apply(as.data.frame(MetaGenome),2,
                               function(x) scoreDiff_weighted(as.data.frame(GenomeVector),x,myconsortia))
  
  	df_W <- ldply(b_percentage_list_W, data.frame)
  	names(df_W) <- c("MetaGenome","BacterialGenome","Percentage","Score","ScorePerIteration",
                   "pfamInGenomePerItration","exclusiveMismatchPerIteration","RelativeMatchPerIteration")
  
  	df_W <- df_W %>% dplyr::group_by(MetaGenome) %>% dplyr::mutate(Index=1:length(MetaGenome))
#To plot metagenome coverage  
  #ggplot(df_W,aes(Index,Percentage,group=MetaGenome,color=MetaGenome))+geom_point()+geom_line()
  
  	df_W <- df_W %>% dplyr::group_by(MetaGenome) %>%  dplyr::mutate(kneepoint=FindKnee_uik(x=Index,y=Percentage))
  
#to plot knee point
    	ggplot(df_W,aes(Index,Percentage,group=MetaGenome,color=MetaGenome))+
    	geom_point()+geom_line()+geom_vline(xintercept = unique(df_W$kneepoint))+theme(legend.position = "none")
  
   	PiBc_knee_point_matchmismatch_2 <-  df_W %>% dplyr::group_by(MetaGenome) %>% dplyr::summarise(k=unique(kneepoint))
write.table(PiBc_knee_point_matchmismatch_2,file = "PiBc_knee_point2_mimicB.txt",sep = "\t",quote = FALSE)
  
  	df_W_keep <- df_W %>% dplyr::group_by(MetaGenome) %>%
    	dplyr::slice(1:unique(kneepoint))
  
  	df_W_keep$Index <- NULL
  	df_W_keep$kneepoint <- NULL
  #setdiff(df_W$BacterialGenome,df$BacterialGenome)  

  ###################################################################

  #total number of pfam in genome and metagenome 
  	GenomeBinary <- read.table("PiBC_Genome_Binary_Vector_111_update_3oct_2019.txt",sep = "\t", header = T)
  	row.names(GenomeBinary) <- GenomeBinary$PfamID # if coulmn name has the pfam ids and rownames does not exists.
  	GenomeBinary$PfamID <- NULL
  	total <- as.data.frame(colSums(cbind(GenomeBinary,MetaGenome)))  #skip description folder.and sums up all pfam for each genome.

  	colnames(total) <- "PfamTotal" #change name of columns # total number of pfams in each genome and metagenome.
  	total$genomes <- rownames(total) #adding column name to total score file.
  	rownames(total) <- NULL

  	temp <- cbind.data.frame(old=c(total$genomes),neu=
                             c(rep("BacterialGenome",ncol(GenomeBinary)),
                               rep("MetaGenome",ncol(MetaGenome))))#taging bacterial genome and metaGenome
  	total$name<- as.character(temp[match(total$genomes,temp$old),2]) #match two differenent dataset and add another columns in one.


##******************** Adding score to pfam table ********************************
#compare each genome from reference genome dataset to metaegenome
  	scorePfam <- function(genomes, metagenome){
    		g <- as.matrix(genomes)
    		y <- metagenome
#res <- apply(g, 2, function(x) sum(y[x==1]))
    		res <- apply(g, 2, function(x) sum(y[x>=1])) #for weighted #count of total number of pfam in bacteria and metagenome
    		return(res)
  		}
  	scorePfamTable <- as.data.frame(apply(as.data.frame(MetaGenome),2,
                                        function(x) scorePfam(as.data.frame(GenomeBinary),x))) # for apply functi
 write.table(scorePfamTable,file="PiBC_pfam_Match_to_MetaGenome.txt",sep = "\t",quote = FALSE)
  
	scorePfamTable$BacterialGenome <- rownames(scorePfamTable)

	finalStats <- left_join(df_W_keep,total[,c("PfamTotal","genomes")], by=c("BacterialGenome" = "genomes")) %>%
    	left_join(.,reshape2::melt(scorePfamTable,id.vars="BacterialGenome",value.name="AbsoluteMatch",variable.name="MetaGenome"),by=c("BacterialGenome","MetaGenome")) %>%
  	dplyr::rename(NovelMatch = Score, BactTotalPfam = PfamTotal) %>%
  	dplyr::group_by(MetaGenome) %>% 
  	dplyr::mutate(AbsoluteMismatch = BactTotalPfam - AbsoluteMatch,
                CumAbsoluteMismatch = cumsum(AbsoluteMismatch),
                CumNovelMatch = cumsum(NovelMatch),
                CumIterationRelativeMatch = cumsum(RelativeMatchPerIteration),
                AbsoluteRelativeMatch = AbsoluteMatch/AbsoluteMismatch,
                CumAbsoluteRelativeMatch = cumsum(AbsoluteRelativeMatch),
                CumExclusiveMismatchPerIteration = cumsum(exclusiveMismatchPerIteration)) %>% ungroup()

# save files using dplyr group functions. from dplyr version 0.8.0 and above the group_walk function was implemented to iterate on grouped tibbles and there is no need to use do()

	dplyr::group_by(finalStats, MetaGenome) %>% 
  	dplyr::group_walk(~ write.table(.x, paste0(unique(.y$MetaGenome), "mimic_B.txt"), row.names = F, quote = F, sep="\t"))

write.table(finalStats,file = paste("finalStats_Knee_2.txt",sep = "_"),sep = "\t",row.names = FALSE,quote = FALSE)

########### Ends here    ################################
