
#knee point finding functions

FindKnee_uik <- function(x,y){
  suppressMessages(require(inflection))
  cut<-uik(x,y)
  return(cut)
}


FindKnee_uik2 <- function(x,y){
  suppressMessages(require(inflection))
  cut<-d2uik(x,y)
  return(cut)
}



# load the file
#a <- read.table("~/Downloads/mbarc_mbarc_stats_mimic_B.txt",sep="\t",header = T)

#Arrange the data according to the score to consider for plotting
#b <- dplyr::arrange(a,CoveragePercentage)

#this column is added to generate numeric x-axis requried for inflection::uik function
#b$index <- as.numeric(rownames(b))

#Find the knee point using two different methods
#kneeuik <- FindKnee_uik(b$index,b$CumNovelMatch)
#kneemclust <- FindKnee_mclust(b$index,b$CumNovelMatch)

#plot the found knee points
#ggplot(b,aes(index,CumNovelMatch))+
 # geom_point()+
  #geom_vline(xintercept = kneeuik, col="red")+
  #geom_vline(xintercept = kneemclust, col="blue")
  
