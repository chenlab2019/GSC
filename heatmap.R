args<-commandArgs(TRUE)
suppressMessages(library(DESeq2))
suppressMessages(library( "gplots" ))
library( "RColorBrewer" )
suppressMessages(library("pheatmap"))
input = args[1]
output = args[2]

data<-readRDS("vst_50_new2.rds")
type<-read.table("new.cluster.50",header=F)
classfication<-read.table("new.tcga2",header=T,sep="\t")
colnames(classfication)[1]<-"cellline"
colnames(classfication)[2]<-"subtype"
colnames(type)[1]<-"cellline"
colnames(type)[2]<-"cluster"
merged<-merge(classfication,type,by="cellline")
data2<-data[,1:50]
#rownames(data2)<-paste(data$chr,data$start,data$end,sep="_")
num<-NULL
for( i in colnames(data2)){
    order<-which(merged[,1] == i)
    num<-c(num,order)
}
merged_new<-merged[num,]
merged_new2<-data.frame(subtype=factor(merged_new[,2]),cluster=factor(merged_new[,3]))
rownames(merged_new2)<-merged_new[,1]

diff<-read.table(input,header=F)
target<-data[as.character(diff[,1]),]
mergerd_new3<-merged_new2[order(merged_new2$cluster,decreasing=F),]
mergerd_new3$cluster <- factor(mergerd_new3$cluster)
target2<-target[,order(merged_new2$cluster,decreasing=F)]

color = colorRampPalette(brewer.pal(3,"RdBu"))(21)
breaks<-seq(-1.5,1.5,by=0.15)
pdf(output)
#out<-pheatmap(target2, show_rownames=F, cluster_rows=T,cluster_cols = F,clustering_distance_cols="euclidean",scale="row",show_colnames= T,color = rev(color),annotation_col=mergerd_new3,clustering_distance_rows="euclidean", cex=1, clustering_method="complete", border_color=FALSE,fontsize_row=0.2,fontsize_col=4,breaks=breaks,annotation_colors = list(subtype=c(CL="#F8766D",MS="#00BA38",PN="#619CFF"),cluster=c("1"="#7570B3","2"="#D95F02","3"="#E6AB02")))
#cluster=c(1="#fc8d59",2="#ffffbf",3="#99d594")
pheatmap(target2, show_rownames=F, cluster_rows=F,cluster_cols = F,clustering_distance_cols="euclidean",scale="row",show_colnames= T,color = rev(color),annotation_col=mergerd_new3,clustering_distance_rows="euclidean", cex=1, clustering_method="complete", border_color=FALSE,fontsize_row=0.2,fontsize_col=4,breaks=breaks,annotation_colors = list(subtype=c(CL="#1600ff",MS="#ff0000",PN="#811788"),cluster=c("1"="#8dd3c7","2"="#fb8072","3"="#fdb462")),cellheight=0.02,cellwidth=3)
pheatmap(target2, show_rownames=F, cluster_rows=T,cluster_cols = F,clustering_distance_cols="euclidean",scale="row",show_colnames= T,color = rev(color),annotation_col=mergerd_new3,clustering_distance_rows="euclidean", cex=1, clustering_method="complete", border_color=FALSE,fontsize_row=0.2,fontsize_col=4,breaks=breaks,annotation_colors = list(subtype=c(CL="#1600ff",MS="#ff0000",PN="#811788"),cluster=c("1"="#8dd3c7","2"="#fb8072","3"="#fdb462")),cellheight=0.02,cellwidth=3)
pheatmap(target2, show_rownames=F, cluster_rows=T,cluster_cols = T,clustering_distance_cols="euclidean",scale="row",show_colnames= T,color = rev(color),annotation_col=mergerd_new3,clustering_distance_rows="euclidean", cex=1, clustering_method="complete", border_color=FALSE,fontsize_row=0.2,fontsize_col=4,breaks=breaks,annotation_colors = list(subtype=c(CL="#1600ff",MS="#ff0000",PN="#811788"),cluster=c("1"="#8dd3c7","2"="#fb8072","3"="#fdb462")))
dev.off()

