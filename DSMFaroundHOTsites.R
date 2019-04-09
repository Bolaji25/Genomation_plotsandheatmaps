library(GenomicRanges)
library(genomation)
library(rtracklayer)
library(RColorBrewer)
library(BSgenome)

#setwd("Path_to_a_directory_higher_than_data")
colfunc <- colorRampPalette(c("black", "yellow","blue"))

#load data
DSMF<- import.bw("../Bigwiggle/sevinc/sm_allCs_combined_M_ucsc.bw")

#Load the feature files
hotsites<-readRDS("../featurefiles/HOT_sites_ce11_2kb.RDS")
hotsites
i<-which(seqnames(hotsites)=="chrX")

xChr<-hotsites[i]

j<-which(seqnames(hotsites)!="chrX")
aChr<-hotsites[j]

export(xChr,"../featurefiles/hotsiteschrX.gtf","gtf")
export(aChr,"../featurefiles/hotsiteschrA.gtf","gtf")
export(hotsites, "../featurefiles/hotsites_ce11_2kb.gtf", "gtf")
hotsitesA<- import("../featurefiles/hotsiteschrA.gtf", "gtf")
hotsitesX<- import("../featurefiles/hotsiteschrX.gtf", "gtf")
hotsites<- import("../featurefiles/hotsites_ce11_2kb.gtf")
#Extend gr to 1000bp on each side of the hotsites sites
flanking_hotsites_A <- resize(hotsitesA, 1000, fix="center", use.names=TRUE)
flanking_hotsites_X <- resize(hotsitesX, 2000, fix="center", use.names=TRUE)
flankingallhotsites <- resize(hotsites, 2000, fix="center", use.names=TRUE)
#Extend gr to 100bp on each side of the rex sites
#flanking_hotsites_A100 <- resize(hotsitesA, 200, fix="center", use.names=TRUE)
#flanking_hotsites_X100 <- resize(hotsitesX, 200, fix="center", use.names=TRUE)
#flankingallhotsites100 <- resize(hotsites, 200, fix="center", use.names=TRUE)


#Make score matrix for 2000bp window
sm_DSMF_hotsitesA = ScoreMatrix(target = DSMF,
                           windows = flanking_hotsites_A ,
                           strand.aware = TRUE, weight.col = "score")
sm_DSMF_hotsitesX = ScoreMatrix(target = DSMF,
                           windows = flanking_hotsites_X,
                           strand.aware = TRUE, weight.col = "score")
sm_DSMF_hotsites = ScoreMatrix(target = DSMF,
                                windows = flankingallhotsites,
                                strand.aware = TRUE, weight.col = "score")
saveRDS(sm_DSMF_hotsitesA, "../featurefiles/sm_DSMF_hotsitesA.RDS")
saveRDS(sm_DSMF_hotsitesX, "../featurefiles/sm_DSMF_hotsitesX.RDS")
saveRDS(sm_DSMF_hotsites, "../featurefiles/sm_DSMF_hotsitesall.RDS")
#Make score matrix for 200bp window
#sm_DSMF_hotsitesA100 = ScoreMatrix(target = DSMF,
 #                               windows = flanking_hotsites_A100 ,
  #                              strand.aware = TRUE, weight.col = "score")
#sm_DSMF_hotsitesX100 = ScoreMatrix(target = DSMF,
 #                               windows = flanking_hotsites_X100,
  #                              strand.aware = TRUE, weight.col = "score")
#sm_DSMF_hotsites100 = ScoreMatrix(target = DSMF,
 #                              windows = flankingallhotsites100,
  #                             strand.aware = TRUE, weight.col = "score")

#Heatmap plots for 2000bp window
heatMatrix(sm_DSMF_hotsitesA,
           xcoords = c(-1000,1000),
           winsorize = c(1,99),
           col = colfunc(10))
heatMatrix(sm_DSMF_hotsitesX,
           xcoords = c(-1000,1000),
           winsorize = c(1,99),
           col = colfunc(10))
heatMatrix(sm_DSMF_hotsites,
                xcoords = c(-1000,1000),
                winsorize = c(1,99),
                col = colfunc(10))
#Heatmap plots for 200bp window
#A100<-heatMatrix(sm_DSMF_hotsitesA100,
 #             xcoords = c(-100,100),
  #            winsorize = c(1,99),
   #           col = colfunc(10))
#X100<-heatMatrix(sm_DSMF_hotsitesX100,
 #             xcoords = c(-100,100),
  #            winsorize = c(1,99),
   #           col = colfunc(10))
#all100<-heatMatrix(sm_DSMF_hotsites100,
 #               xcoords = c(-100,100),
  #              winsorize = c(1,99),
   #             col = colfunc(10))

################################################################################################
#######Plotting with ggplot####################################################################
#data <- readRDS("../featurefiles/sm_DSMF_hotsitesA.RDS")
#data <- readRDS("../featurefiles/sm_DSMF_hotsitesX.RDS")
data <- readRDS("../featurefiles/sm_DSMF_hotsitesall.RDS")
library(ggplot2)
library(genomation)
library(tidyr)
####################################################
# Plotting function for genomation ScoreMatrix     #
####################################################
#' @importFrom ggplot2 ggplot
#' @importFrom genomation ScoreMatrix
#' @importFrom tidyr gather
#' @param data   \code {genomation} ScoreMatrix of dSMF values[0-1]
#' @param myYlab {text} Y label for graph
#' @param myXlab {text} X label for graph
#' @param feature_label {text} middle value on which windowns are aligned to
#' @param title {text} graph titleDistance between the starts of consecutive windows
#' 
plotAveragedSMF<-function(data,         #genomation ScoreMatrix
                          myYlab,       #label for Y axis (type of loci)
                          myXlab="CpG/GpC position", #label for X axis
                          feature_label="HOT sites",
                          title=NULL)
{
  #Transform genomation matrix to data frame
  df <- as.data.frame(data)
  #Remove 0 values
  df[df==0] <- NA
  #rename columns, get number of columns
  colnames(df)<-as.numeric(gsub("V","",colnames(df)))
  width_data <- max(as.numeric(colnames(df)))
  #name samples
  rownames(df)<- as.numeric(rownames(df))
  #prepare data for plotting
  d<- tidyr::gather(df,key=position, value=methylation)
  d$molecules<-rownames(df)
  d$position<-as.numeric(d$position)
  
  p <- ggplot2::ggplot(d,aes(x=position,y=molecules,width=1)) +
    ggplot2::geom_tile(aes(width=2,fill=methylation)) +
    ggplot2::scale_fill_gradientn(values = c(0,1),
                                  colors=c("yellow", "red"),
                                  na.value="black") +
    ggplot2::theme(panel.background = element_rect(fill="black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   axis.text.y = element_blank(),
                   axis.ticks.y =  element_blank(),
                   plot.title = element_text(face = "bold",hjust = 0.5)) +
    # ggplot2::guides(colour=guide_colorbar(),
    #       fill = guide_legend(title="dSMF", 
    #                          title.position="top")) + 
    ggplot2::scale_x_continuous(name=myXlab, 
                                limits = c(0,width_data),
                                expand=c(0,0),
                                breaks=c(0,width_data/2,width_data),
                                labels=c(paste0("-",width_data/2),feature_label,paste0(width_data/2))) + 
    ggplot2::ylab(myYlab) + 
    ggplot2::ggtitle(title)
  return(p)
}

#pdf("dSMF_at_allHOT_sites.pdf")
#plotAveragedSMF(data,myYlab="Genes",myXlab="Position",feature_label="HOT sites",title="dSMF at HOT sites")
ggsave("dsmfHOTsites.tiff",
       plotAveragedSMF(data,myYlab="Genes",myXlab="Position",feature_label="TSS",title="dSMF at HOT sites"), 
       units="mm", width=dim(data)[2]/10,height=dim(data)[1]/10)
dev.off()
