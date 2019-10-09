library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Celegans.UCSC.ce11)
library(IRanges)
library("GGally")
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
#read in the public files' table, always set header so that the header is defined
#publicfiles<-read.table("./publicfiles.txt", header=TRUE, stringsAsFactors = FALSE)
#publicfiles
#for (i in 1:dim(publicfiles)[1]) {
# download.file(publicfiles$Sample[i],destfile=paste0(publicfiles$FileName[i],"_ce10.bw"))
#}
#import the bigwig files
#pub<-import("atac_wild_l3_ce10.bw")
#pub<-import(paste0(publicfiles$FileName[i],"_ce10.bw"))
#pub

#for the loops
#for (i in 1:dim(publicfiles)[1]) {
 # pub<-import(paste0(publicfiles$FileName[i],"_ce10.bw"))
  
  
  #import the downloaded chain file
  #chainFile<-import.chain("./ce10ToCe11.over.chain")
  
  #lift the files from ce10 to ce11
  #pub_ce11<-unlist(liftOver(pub, chain = chainFile))
  #to make sure all chromosomes are present in the seqlevels and seqlengths (get from BSgenome)
  #seqlevels(pub_ce11)<-seqlevels(Celegans)
  #seqlengths(pub_ce11)<-seqlengths(Celegans)
  
  
  #to export the files
  #export(pub_ce11, paste0(publicfiles$FileName[i],"_ce11.gtf"),  "gtf")
  #export(pub_ce11, paste0(publicfiles$FileName[i], "_ce11.bw"), "bw")
#}

#make GRanges object from BSgenome (accessed with "Celegans") (i remove the mito chr)
genomeGR<-GRanges(seqnames=seqnames(Celegans)[1:6],ranges=IRanges(start=1, end=seqlengths(Celegans)[1:6]), strand="*")

# get a list of all the files with the enrichment data
Datafiles<-list.files("/Users/imac/Desktop/Bolaji/PhD/analysis/Bioinformatics_scripts/Bigwiggle/","bw")
Datafiles
#create GRangesList object with different sized tiles along genome
tile1Mb<-unlist(tile(x=genomeGR,width=1000000))
tile100kb<-unlist(tile(x=genomeGR,width=100000))
#tile10kb<-unlist(tile(x=genomeGR,width=10000))
#tile1kb<-unlist(tile(x=genomeGR,width=1000))
#tile100bp<-unlist(tile(x=genomeGR,width=100))
tile50bp<-unlist(tile(x=genomeGR,width=50))
tile10bp<-unlist(tile(x=genomeGR,width=10))
tileList<-list("tile1Mb"=tile1Mb,"tile100kb"=tile100kb,"tile10kb"=tile10kb,"tile1kb"=tile1kb,"tile100bp"=tile100bp, "tile50bp"=tile50bp, "tile10bp"=tile10bp)

# read in the bedgraph files in pairs to calculate enrichment counts at different scales
for (f in 1:length(Datafiles)) {
  #get normalised counts
  Data<-import(paste0("/Users/imac/Desktop/Bolaji/PhD/analysis/Bioinformatics_scripts/Bigwiggle/",Datafiles[f]))

  #make a GRangesList from GRanges
  for (i in 1:length(tileList)) {
    tiles<-split(tileList[[i]],seqnames(tileList[[i]]))
    
    # convert GRanges back to coverage
    rle<-coverage(Data,weight="score")
    if (length(rle)>6){
      j<-grep("M",names(rle))
      rle<-coverage(Data,weight="score")[-j]
    }
    
    names(rle)<-gsub("chr","", names(rle))
    names(rle)<-paste0("chr",names(rle))
    #seqlevels(rle)<-paste0("chr", seqlevels(rle))
    
    #create views corresponding to the tiles on the coverage rle and average counts
    summedCov<-viewMeans(RleViewsList(rangesList=tiles,rleList = rle))
    
    #now save the values into tiles
    mcols(tileList[[i]])[,Datafiles[f]]<-as.data.frame(unlist(summedCov))
  }
}


#to make correlations without ordering based on similariies of correlation
#Mb1Tile<- ggcorr(as.data.frame(mcols(tileList[[1]])),label=T,hjust = 1, size = 5, color = "black", layout.exp = 7, method= c("everything","spearman"))+ggtitle("1Mbtilesize")
#bp100Tile<- ggcorr(as.data.frame(mcols(tileList[[5]])),label=T,hjust = 1, size = 5, color = "black", layout.exp = 7, method= c("everything","spearman"))+ggtitle("100bptilesize")
#bp50Tile<- ggcorr(as.data.frame(mcols(tileList[[6]])),label=T,hjust = 1, size = 5, color = "black", layout.exp = 7, method= c("everything","spearman"))+ggtitle("50ptilesize")
#bp10Tile<- ggcorr(as.data.frame(mcols(tileList[[7]])),label=T,hjust = 1, size = 5, color = "black", layout.exp = 7, method= c("everything","spearman"))+ggtitle("10bptilesize")
#bp10Tileperas<- ggcorr(as.data.frame(mcols(tileList[[7]])),label=T,hjust = 1, size = 5, color = "black", layout.exp = 7)+ggtitle("10bptilesize")
#kb10Tile<- ggcorr(as.data.frame(mcols(tileList[[3]])),label=T,hjust = 1, size = 5, color = "black", layout.exp = 7, method= c("everything","spearman"))+ggtitle("10kbtilesize")
#kb100Tile<- ggcorr(as.data.frame(mcols(tileList[[2]])),label=T,hjust = 1, size = 5, color = "black", layout.exp = 7, method= c("everything","spearman"))+ggtitle("100kbtilesize")
#kb1Tile<- ggcorr(as.data.frame(mcols(tileList[[4]])),label=T,hjust = 1, size = 5, color = "black", layout.exp = 7, method= c("everything","spearman"))+ggtitle("1kbtilesize")
#correlations<-grid.arrange(bp50Tile,bp10Tile, ncol=1, nrow=2)
#ggsave("correlationspage4.pdf", plot =correlations)




#to do the tile correlation using gapmap to be able to grop them based on similarities(i.e positives together and negatives together)
library(gapmap)
library(GGally)
library(tidyr)

# Basketball statistics provided by Nathan Yau at Flowing Data.
#dt <- read.csv("http://datasets.flowingdata.com/ppg2008.csv")

# create unordered correlation matrix
#corMat<-ggcorr(dt[, -1])
#corMat
corMat<-ggcorr(as.data.frame(mcols(tileList[[3]])))
corMat

# extract the data from the corMat object
cordata<-corMat$data
names(cordata)

#remove label column
cordata$label<-NULL

# add back missing values from other half of triangular matrix
temp<-cordata
temp$x<-cordata$y
temp$y<-cordata$x

#merge the two objects
corlong<-rbind(cordata,temp)

#rearrange the data table to a wide format
corwide<-spread(corlong,key=y,value=coefficient)

# move the row names from the x column to rownames
rowLabels<-corwide$x
corwide$x<-NULL
rownames(corwide)<-rowLabels

#fill in the diagonal identity values
corwide[is.na(corwide)]<-1

#calculate distance matrix
cordist<-dist(corwide)


cor.matrix <- as.matrix(corwide)
cor.dendro <- as.dendrogram(hclust(d = dist(x = cor.matrix)))
cor.order <- order.dendrogram(cor.dendro)

newDF<-as.data.frame(mcols(tileList[[1]]))[,cor.order]
Mb1Tile<- ggcorr(newDF,
                 label=T,hjust = 1, size = 5, color = "black",
                 layout.exp = 8, method= c("everything","spearman"))+ggtitle("1Mbtilesize")
correlations<-grid.arrange(Mb1Tile,kb100Tile, ncol=1, nrow=2)
ggsave("correlationsordered1mb100kb.pdf", plot =correlations)
# perform heirarchical clustering
#corhc<-hclust(cordist)

#extract dendogram
#dend <- as.dendrogram(corhc)

# create a colour scale generating function with the same colours as ggcorr
col_func = colorRampPalette(c("#3B9AB2","#EEEEEE","#F21A00"))
# plot the clustered correlation matrix
#gapmap(m = as.matrix(cordist), d_row= rev(dend), d_col=dend, col = col_func(11), mapping="linear")
  
