suppressMessages(library(DESeq2))
suppressMessages(library(edgeR))
library("Hmisc")
library(reshape2)
library(ggplot2)
library( "RColorBrewer" )
library(pheatmap)

args<-commandArgs(T)
input<-args[1]
a<-read.table(input,header=T,sep="\t")

dat_edgeR_norm <- cpm(a)
exprs<-log(x=(dat_edgeR_norm + 1),base=2)
time<-factor(cell104_type_mid[,1])
count=1
for ( i in unique(time) ){
  if ( dim(as.data.frame(exprs[,which(time==i)]))[2] == 1 ){
    mean_rpkm=data.frame(exprs[,which(time==i)])
  } else {
    mean_rpkm=data.frame(rowMeans(exprs[,which(time==i)]))
  }
  colnames(mean_rpkm)=i
  if (count == 1){
    mean_rpkm_ok=mean_rpkm
  } else {
    mean_rpkm_ok=merge(mean_rpkm_ok,mean_rpkm,by="row.names")
    rownames(mean_rpkm_ok)=mean_rpkm_ok[,1]
    mean_rpkm_ok=mean_rpkm_ok[,-1]
  }
  count=count+1
}
exprs_with_time=as.matrix(mean_rpkm_ok, header=TRUE, sep="\t",row.names=1,as.is=TRUE)

expr<-exprs_with_time

res2<-rcorr(as.matrix(expr),type = "pearson")
get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
}

# Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }

reorder_cormat <- function(cormat){
# Use correlation between variables as distance
dd <- as.dist((1-cormat)/2)
hc <- hclust(dd)
cormat <-cormat[hc$order, hc$order]
}

cormat <- res2$r
min(cormat)
anno4<-read.table("cluster.new.tcga",header=T,sep="\t")

num<-NULL
for ( i in as.character(colnames(cormat))){
    a<-which(anno4$cellline == i)
    num<-c(num,a)
}
anno3<-anno4[num,]

anno2<-data.frame(subtype=anno3$tcga,atac=anno3$new_cluster)
anno2$atac<-as.factor(anno2$atac)
rownames(anno2)<-anno4$cellline

breaks<-seq(0.5,1,by=0.025)
colour<-colorRampPalette(brewer.pal(7,"RdYlBu"))(21)
pdf("sample.correlation2.pdf")
pheatmap(cormat,fontsize_row=4,show_colnames= F,treeheight_col = 0,show_rownames=T,breaks=breaks,color=rev(colour),annotation_row=anno2,annotation_colors = list(subtype=c(CL="#1600ff",MS="#ff0000",PN="#811788"),atac=c("1"="#1f78b4","2"="#976D0B","3"="#FFAE33")),border_color = "NA",annotation_col=anno2)
dev.off()

out<-pheatmap(cormat,fontsize_row=4,show_colnames= F,treeheight_col = 0,show_rownames=T,breaks=breaks,color=rev(colour),annotation_row=anno2,annotation_colors = list(subtype=c(CL="#F8766D",MS="#00BA38",PN="#619CFF")))
write.table(colnames(cormat[,out$tree_col[["order"]]]),file="order.all_peaks.txt",row.names=F,quote=F)

