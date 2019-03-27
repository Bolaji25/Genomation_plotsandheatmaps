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
flanking_hotsites_A <- resize(hotsitesA, 2000, fix="center", use.names=TRUE)
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
saveRDS(heatMatrix(sm_DSMF_hotsitesA,
                   xcoords = c(-1000,1000),
                   winsorize = c(1,99),
                   col = colfunc(10)), "../featurefiles/sm_DSMF_hotsitesA.RDS", "RDS")
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
A<-heatMatrix(sm_DSMF_hotsitesA,
           xcoords = c(-1000,1000),
           winsorize = c(1,99),
           col = colfunc(10))
X<-heatMatrix(sm_DSMF_hotsitesX,
           xcoords = c(-1000,1000),
           winsorize = c(1,99),
           col = colfunc(10))
all<-heatMatrix(sm_DSMF_hotsites,
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
