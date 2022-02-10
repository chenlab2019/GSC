suppressMessages(library(DESeq2))
suppressMessages(library(edgeR))
library("Hmisc")
library(reshape2)
library(ggplot2)
library( "RColorBrewer" )
library(pheatmap)

cell_name<-read.table("/disk1/xilu/glioblastoma/cell_line_140/140.sample.name2",header=F,sep="\t")
data<-read.table("/disk1/xilu/glioblastoma/cell_line_140_new_pipeline_call_peaks/138.sample.peak.remove.chrun.random.chry.xls",header=T,sep="\t")
peak_name<-paste(data[,1],data[,2],data[,3],sep="_")
data2<-data[,4:141]
rownames(data2)<-peak_name
colnames(data2)<-cell_name[,1]

simple_name<-read.table("/disk1/xilu/glioblastoma/cell_line_140/140.sample.simple.name2",header=F,sep="\t")
cell51_type<-read.table("/disk1/xilu/glioblastoma/cell_line_140_new_pipeline_call_peaks/promoter_enhancer_nmf_50/cellline_exclude3175_new.txt2",header=T,sep="\t")
simple_name_140<-paste(rep("U",138),simple_name[,1],sep="")
cell140_type<-read.table("/disk1/xilu/glioblastoma/cell_line_140/140SampleNameOrder_NEW.txt2",header=T,sep="\t")

num<-NULL
for ( i in as.character(cell51_type[,1])){
    a<-which(simple_name_140 == i)
    num<-c(num,a)
}
cell_104_data<-data2[,num]

type_num<-NULL
for( i in as.character(cell51_type[,1])){
    b<-which(cell140_type[,1] == i)
    type_num<-c(type_num,b)
}
cell104_type_mid<-cell140_type[type_num,]
cell104_type<-data.frame(subtype=cell104_type_mid$subtype_latest)
rownames(cell104_type)<-colnames(cell_104_data)
cell104_type$subtype<-factor(cell104_type$subtype)

print ("begin")
b<-NULL
#for(i in 1:dim(data2)[1]){
#    a<-ifelse(data2[i,] > 5,1,0)
#    b<-rbind(b,a)
#}
data3<-(data2 > 5) 
for(i in 1:ncol(data3)){
    if(is.logical(data3[, i]) == TRUE) data3[, i] <- as.numeric(data3[, i])
    }

print ("end")
cell_104_data<-data2[which(rowSums(data3) > 5),]

#all<-read.table("/disk1/xilu/glioblastoma/cell_line_51/pearson/modefied.peak.remove.random.chrUn.bed",header=F,sep="\t")
#all2<-paste(all[,1],all[,2],all[,3],sep="_")
#cell_104_data2<-cell_104_data[as.character(all2),]


dat_edgeR_norm <- cpm(cell_104_data)
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
anno4<-read.table("/disk1/xilu/glioblastoma/corrected_order_50/cluster.new.tcga",header=T,sep="\t")

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

