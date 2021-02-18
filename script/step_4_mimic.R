# Author Neeraj kumar
#update feburary 2021

# Description -------------------------------------------------------------
#This script is used to calculate the minimal microbial consrotia for a metagenome based on the reference pfam vector of microbial species

#minimum number of species covering maximum number of metagenome

# usage -------------------------------------------------------------------

#  Rscript --vanilla step_4_mimic.txt -m metagenomefileName -g genomeVectorFileName -i iterationNumber -o mimicOutputName -k kneepointbasedOutputName
 
  start_time <- Sys.time()
# install packages --------------------------------------------------------
  
  paks <- c("tidyverse","ggsignif","cowplot","RColorBrewer", "ggcorrplot","vegan",
            "ade4","colorspace","plyr","data.table","ggplot2","inflection","dplyr",
            "readr","optparse")
  paks<-lapply(paks, function(x) suppressMessages(require(x, character.only = TRUE)))
  rm(paks)
  
# Mkaing opotions ---------------------------------------------------------
  #library("optparse")
  # maiking option using package "optparse"
  option_list = list(
    make_option(c("-m", "--MetaGenomeVector"), type="character", default=NULL, 
                help="metagenome vector file name", metavar="character"),
    make_option(c("-g", "--Genome_Vector"), type="character", default=NULL, 
                help="genome vector file name", metavar="character"),
    make_option(c("-i", "--itrationNumber"), type="integer", default=10, 
                 help="iteratioin number [default= %default]", metavar="number"),
    make_option(c("-o", "--MiMiC"), type="character", default="MiMiC.txt", 
                help="MiMiC file name [default= %default]", metavar="character"),
    make_option(c("-k", "--MiMiC_KneePoint"), type="character", default="MiMiC_KneePoint.txt", 
                help="MiMiC file name at knee point [default= %default]", metavar="character")
    ) 
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
  # if (is.null(opt$file)){
  #   print_help(opt_parser)
  #   stop("At least one argument must be supplied (input file).n", call.=FALSE)
  # }

# input  tables------------------------------------------------------------------
  # reading metagenome vector file
  
  MetaGenome <- data.table::fread(opt$MetaGenomeVector,sep = "\t",
                                  header = TRUE ,
                                  stringsAsFactors = FALSE,
                                  check.names = FALSE)
  MetaGenome <- as.data.frame(MetaGenome) 
  row.names(MetaGenome) <- MetaGenome$PfamID # to make sure if rows name remain same 
  MetaGenome$PfamID <- NULL
  # reading genome vector file
	GenomeVector <- data.table::fread(opt$Genome_Vector,
	                         sep = "\t",
	                         header = T,
	                         stringsAsFactors = FALSE)
	GenomeVector <- as.data.frame(GenomeVector)
	row.names(GenomeVector) <- GenomeVector$PfamID # skip if row names alreay assigned and PfamID is not available
  GenomeVector$PfamID <- NULL
	
	  # GenomeVector <- data.table::fread("GenomeVector_outlierExclude_22807.txt",
	  #                            sep = "\t",header = T,stringsAsFactors = FALSE)
	  # GenomeVector <- as.data.frame(GenomeVector)
	  # row.names(GenomeVector) <- GenomeVector$PfamID # skip if row names alreay assigned and PfamID is not available
	  # GenomeVector$PfamID <- NULL
	  # MetaGenome <- data.table::fread("PfamVector.txt",sep = "\t",
	  #                                 header = TRUE ,stringsAsFactors = FALSE,check.names = FALSE)
	  # MetaGenome <- as.data.frame(MetaGenome)
	  # row.names(MetaGenome) <- MetaGenome$PfamID # to make sure if rows name remain same 
	  # MetaGenome$PfamID <- NULL
	  # 
	 
# consortia size ----------------------------------------------------------

#	passing the number of iteration
	Minimal_consortia_size=opt$itrationNumber

	# myconsortia<-10 # minimum consortia size
	# Minimal_consortia_size<-100
	# 

# Functions ---------------------------------------------------------------

	#knee point function uik
		FindKnee_uik <- function(x,y){
	  suppressMessages(require(inflection))
	  cut<-uik(x,y)
	  return(cut)
	}
	
	#knee point function d2uik
	FindKnee_uik2 <- function(x,y){
	  suppressMessages(require(inflection))
	  cut<-d2uik(x,y)
	  return(cut)
	}
  # MiMiC function to calculate minimal microbial consortia
	#this function compare the genome vector profile to the metagenome vector profile 
		scoreDiff_weighted <- function(genomes, metagenome,Minimal_consortia_size){
	  g <- as.matrix(genomes)
	  y <- as.matrix(metagenome)
	  mscore <- sum(y)
	  mscore1 <- sum(y)
	  res <- NULL
	 
	  for(i in 1:Minimal_consortia_size){
	    matchCount <- apply(g, 2, function(x) sum(y[x >= 1.0])) # ---> only count of Pfam match which present in genome and not in metagenome 
	    totalpfamG <- apply(g,2,function(x) sum(x>=1)) # total number of Pfam present in microbial genome 
	    misMatchCount <- totalpfamG - matchCount # number Of Pfams which are present in genome but not in Metagenome.
	    matchratio <- ((matchCount)/(matchCount+misMatchCount)) # (x)/(x+y) # ---> mimic_B             
	    matchratio[matchratio=="NaN"]<- 0 # because of some values become 0, when 0/0 it gives NaN
	    
	    bestHit <- matchCount[names(matchratio[matchratio==max(matchratio)])] #---> genome with maximum match count and maximum match mismatch ratio
	    SelectedBestHit <- names(bestHit[bestHit==max(bestHit)])[1] #---> names of the genome to be selected
	    
	    print(paste ("number of itreation" ,i,sep="=")) #---> printing number of iteration
	   
	    a<-g[,SelectedBestHit] #--> selected genome vector
	   
	    ascore <- bestHit[bestHit==max(bestHit)][1] # score
	   
	    pfamInGenomePerItration = unname(totalpfamG[SelectedBestHit]) # Total number of pfam used in the comparison 
	    exclusiveMatch <- unname(matchCount[SelectedBestHit]) # number of Pfam match to metagenome
	    exclusiveMismatchPerIteration <- unname(misMatchCount[SelectedBestHit]) # number of Pfam that did not match to the metagenome 
	    
	    print(paste("pfam",pfamInGenomePerItration,"exclusiveMatch",exclusiveMatch,"exclusiveMismatchPerIteration",
	                exclusiveMismatchPerIteration,sep = "=")) # Printing the summary on screen at each iteration
	    #print the genome with the name
	    #print (ascore)
	    #unname(ascore)
	    b <- as.character(SelectedBestHit) # name of the selected genome
	    print(b)
	    g <-g[,-which(colnames(g)==b),drop=F] # exclude selected genome.
	    if(i==1){
	      covered <- sum(y[a>=1])
	      rscore <- sum(y[a>=1])
	      percent <- (covered/mscore)*100
	      #uncovered<-(mscore-covered)
	      #uncovered_Percentage <-(uncovered/mscore)
	    } else {
	      rscore <- sum(y[a>=1])
	      covered <- covered+rscore # cumulative Pfam number covered in metagenome
	      percent <- (covered/mscore)*100 # cumulative percentage of Pfam covered in metagenome 
	      #uncovered<-(mscore-covered)
	      #uncovered_Percentage <-(uncovered/mscore)
	    }
	    
	    test <- names(y[a>=1,]) #  Pfams IDs common in genome and metagenome  
	    #print(g)
	    g[test,] <- 0 # making all Pfam as "zero" in reference genome set ( those function were already taken ) 
	    y[a>=1]<-0 # makes all metagenome matching with highest scoring genome as '0'
	   	mscore1 <- sum(y) #number of pfam in metagenome
	   
	    if(rscore==0){break} # if no pfam match found further, then stop the script
	    #RelativeMatchPerIteration=(rscore / exclusiveMismatchPerIteration)
	    
	    res <- rbind(res,cbind.data.frame(bname=SelectedBestHit,
	                                      percentage=percent,score=rscore,ScorePerIteration=ascore,
	                                      pfamInGenomePerItration=pfamInGenomePerItration,
	                                      #RelativeMatchPerIteration=RelativeMatchPerIteration,
	                                      exclusiveMismatchPerIteration=exclusiveMismatchPerIteration)) # making output dataframe
	    
	  }
	  return(res)
	}
	

# Applying MiMiC function -------------------------------------------------

 
  	b_percentage_list <- apply(as.data.frame(MetaGenome),
  	                           2,function(x) scoreDiff_weighted(as.data.frame(GenomeVector),
  	                                                            x,Minimal_consortia_size)) # making a list of the mimic calculation
  
  	df_W <- ldply(b_percentage_list, data.frame) # converting list into dataframe
  	names(df_W) <- c("MetaGenome","BacterialGenome","Percentage","Score","ScorePerIteration",
                   "pfamInGenomePerItration","exclusiveMismatchPerIteration") #naming the column in the dataframe
  
  	df_W <- df_W %>% 
  	  dplyr::group_by(MetaGenome) %>% 
  	  dplyr::mutate(Index=1:length(MetaGenome)) # adding the rank of species in metagenome
  	  #To plot metagenome coverage  
  
  	p1 <- ggplot(df_W,aes(Index,Percentage,group=MetaGenome,color=MetaGenome))+ 
  	  geom_point()+
  	  geom_line() # plot for pfam coverage of metagenome by minimal consortia
    #ggsave("Pfam_coverage_plot",plot=p1)
  	df_W <- df_W %>% 
  	  dplyr::group_by(MetaGenome) %>%  
  	  dplyr::mutate(kneepoint=FindKnee_uik(x=Index,y=Percentage)) # calculating knee point
  
  #to plot knee point
    	p2 <- ggplot(df_W,aes(Index,Percentage,group=MetaGenome,color=MetaGenome))+ 
    	  geom_point()+
    	  geom_line()+
    	  geom_vline(xintercept = unique(df_W$kneepoint))+
    	  theme(legend.position = "none") # plot the knee point on 
  #ggsave("mbarc_mimic_coverage_ncbiRef_2016_kneePoint.pdf",plot=p2)
  knee_point_summary <-  df_W %>% 
    dplyr::group_by(MetaGenome) %>% 
    dplyr::summarise(k=unique(kneepoint))
  #write.table(knee_point_summary,file = "knee_point_summary.txt",sep = "\t",quote = FALSE)

  	df_W_keep <- df_W %>% 
  	  dplyr::group_by(MetaGenome) %>%
    	dplyr::slice(1:unique(kneepoint)) # making seperate dataframe for knee point based selection
  
  	df_W_keep$Index <- NULL
  	df_W_keep$kneepoint <- NULL
  #setdiff(df_W$BacterialGenome,df$BacterialGenome)  


# Total number of pfam and match to metagenome ----------------------------

  #total number of pfam in genome and metagenome 
  #	GenomeBinary <- GenomeVector
  #	row.names(GenomeBinary) <- GenomeBinary$PfamID # if coulmn name has the pfam ids and rownames does not exists.
  #	GenomeBinary$PfamID <- NULL
  	total <- as.data.frame(colSums(cbind(GenomeVector,MetaGenome)))  #skip description folder.and sums up all pfam for each genome.

  	colnames(total) <- "PfamTotal" #change name of columns # total number of pfams in each genome and metagenome.
  	total$genomes <- rownames(total) #adding column name to total score file.
  	rownames(total) <- NULL

  	temp <- cbind.data.frame(old=c(total$genomes),neu=
                             c(rep("BacterialGenome",ncol(GenomeVector)),
                               rep("MetaGenome",ncol(MetaGenome)))) #taging bacterial genome and metaGenome
  	total$name<- as.character(temp[match(total$genomes,temp$old),2]) # number of in all genomes and metagenome


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
                                        function(x) scorePfam(as.data.frame(GenomeVector),x))) # for apply functi
# write.table(scorePfamTable,file="PiBC_pfam_Match_to_MetaGenome.txt",sep = "\t",quote = FALSE)
  
	scorePfamTable$BacterialGenome <- rownames(scorePfamTable)

# final stats without knee point ------------------------------------------


	finalStats <- dplyr::left_join(df_W,total[,c("PfamTotal","genomes")], by=c("BacterialGenome" = "genomes")) %>%
	  dplyr::left_join(.,reshape2::melt(scorePfamTable,
	                                    id.vars="BacterialGenome",
	                                    value.name="AbsoluteMatch",variable.name="MetaGenome"),by=c("BacterialGenome","MetaGenome")) %>%
	  dplyr::rename(NovelMatch = Score, BactTotalPfam = PfamTotal) %>%
  	dplyr::group_by(MetaGenome) %>% 
  	dplyr::mutate(AbsoluteMismatch = BactTotalPfam - AbsoluteMatch,
                CumAbsoluteMismatch = cumsum(AbsoluteMismatch),
                #AbsoluteRelativeMatch = AbsoluteMatch/AbsoluteMismatch,
                CumNovelMatch = cumsum(NovelMatch)
                #CumAbsoluteRelativeMatch = cumsum(AbsoluteRelativeMatch),
                #CumExclusiveMismatchPerIteration = cumsum(exclusiveMismatchPerIteration),
                #CumIterationRelativeMatch = cumsum(RelativeMatchPerIteration),
                ) %>% ungroup()
# in case if seperate file is needed for each metagenome
#saveRDS(finalStats,file = "mbarc_final_stats_2020.rds")
# save files using dplyr group functions. from dplyr version 0.8.0 and above the group_walk function was implemented to iterate on grouped tibbles and there is no need to use do()

#	dplyr::group_by(finalStats, MetaGenome) %>% 
 # 	dplyr::group_walk(~ write.table(.x, paste0(unique(.y$MetaGenome), "mimic_B.txt"), row.names = F, quote = F, sep="\t"))

write.table(finalStats,
            file = opt$MiMiC,sep = "\t",
            row.names = FALSE,
            quote = FALSE)


# final stats with knee point ---------------------------------------------
finalStats_kneePoint <- dplyr::left_join(df_W_keep,total[,c("PfamTotal","genomes")],
                                         by=c("BacterialGenome" = "genomes")) %>%
  dplyr::left_join(.,reshape2::melt(scorePfamTable,id.vars="BacterialGenome",
                                    value.name="AbsoluteMatch",variable.name="MetaGenome"),
                   by=c("BacterialGenome","MetaGenome")) %>%
  dplyr::rename(NovelMatch = Score, BactTotalPfam = PfamTotal) %>%
  dplyr::group_by(MetaGenome) %>% 
  dplyr::mutate(AbsoluteMismatch = BactTotalPfam - AbsoluteMatch,
                CumAbsoluteMismatch = cumsum(AbsoluteMismatch),
                #AbsoluteRelativeMatch = AbsoluteMatch/AbsoluteMismatch,
                CumNovelMatch = cumsum(NovelMatch)
                #CumAbsoluteRelativeMatch = cumsum(AbsoluteRelativeMatch),
                #CumExclusiveMismatchPerIteration = cumsum(exclusiveMismatchPerIteration),
                #CumIterationRelativeMatch = cumsum(RelativeMatchPerIteration),
  ) %>% ungroup()

  write.table(finalStats_kneePoint,
              file = opt$MiMiC_KneePoint,
              sep = "\t",
              row.names = FALSE,
              quote = FALSE) # 
########### Ends here    ################################
	end_time <- Sys.time()
  print(" The processing was taken..")
	end_time - start_time
	# 
