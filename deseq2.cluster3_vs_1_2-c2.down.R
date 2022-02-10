data2<-readRDS("/disk1/xilu/glioblastoma/cell_line_140_new_pipeline_call_peaks/differential_peaks/deseq2.data.rds")
data<-data2[,-which(colnames(data2) %in% c("U3175"))]
type<-read.table("/disk1/xilu/glioblastoma/corrected_order_50/new.cluster.50",header=F)
classfication<-read.table("/disk1/xilu/glioblastoma/cell_line_51/cellline_51_subtype_new.txt",header=T,sep="\t")
colnames(classfication)[1]<-"cellline"
colnames(classfication)[2]<-"subtype"
colnames(type)[1]<-"cellline"
colnames(type)[2]<-"cluster"
merged<-merge(classfication,type,by="cellline")
data2<-data
num<-NULL
for( i in colnames(data2)){
    order<-which(merged[,1] == i)
    num<-c(num,order)
}
merged_new<-merged[num,]
merged_new2<-data.frame(subtype=factor(merged_new[,2]),cluster=factor(merged_new[,3]))
rownames(merged_new2)<-merged_new[,1]
merged_new2$cluster[grep("2",merged_new2$cluster)]<-1
merged_new2$cluster<- factor(merged_new2$cluster)

library(DESeq2)
#atacDDS <- DESeqDataSetFromMatrix(round(data2,digits=0), merged_new2, ~ cluster + subtype)
atacDDS <- DESeqDataSetFromMatrix(round(data2,digits=0), merged_new2, ~ cluster)
#atacDDS_deseq <- DESeq(atacDDS,fitType="local")
#atacDDS_deseq <- DESeq(atacDDS,minReplicatesForReplace = 20,parallel=T)
atacDDS_deseq <- DESeq(atacDDS,minReplicatesForReplace = Inf,parallel=T)
#atac_Rlog <- rlog(atacDDS_deseq,fitType = "local")
result1 <- results(atacDDS_deseq, contrast=c('cluster','3','1'))

save(atacDDS_deseq,file="atacDDS3_vs1_2.RData")
#save(atac_Rlog,file="atac_Rlog_mean.RData")
write.table(result1,file="cluster3_all.results.xls",quote=F,sep="\t")

#result1
maxCooks<-apply(assays(atacDDS_deseq)[["cooks"]], 1, max)
result11<-cbind(result1,maxCooks)
result2<-result11[which(result11$maxCooks < 100),]
resSig <- result2[ which(result2$padj < 0.01 & result2$pvalue < 0.01), ]
#filter 20
resSig2 <- resSig[which(resSig$log2FoldChange < -1 ),]
resSig2 <- resSig2[which(resSig2$baseMean > 20 & resSig2$lfcSE < 1),]
write.table(resSig2,file="cluster3_vs1_2_mean_20.down.xls",quote=F,sep="\t")
#filter 30
resSig2 <- resSig[which(resSig$log2FoldChange < -1 ),]
resSig2 <- resSig2[which(resSig2$baseMean > 30 & resSig2$lfcSE < 1),]
write.table(resSig2,file="cluster3_vs1_2_mean_30.down.xls",quote=F,sep="\t")
#filter 20
resSig2 <- resSig[which(resSig$log2FoldChange < -1 ),]
resSig2 <- resSig2[which(resSig2$baseMean > 40 & resSig2$lfcSE < 1),]
write.table(resSig2,file="cluster3_vs1_2_mean_40.down.xls",quote=F,sep="\t")
