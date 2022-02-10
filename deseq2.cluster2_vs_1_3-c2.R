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
#rownames(data2)<-paste(data$chr,data$start,data$end,sep="_")
num<-NULL
for( i in colnames(data2)){
    order<-which(merged[,1] == i)
    num<-c(num,order)
}
merged_new<-merged[num,]
merged_new2<-data.frame(subtype=factor(merged_new[,2]),cluster=factor(merged_new[,3]))
rownames(merged_new2)<-merged_new[,1]
merged_new2$cluster[grep("1",merged_new2$cluster)]<-3

library(DESeq2)
atacDDS <- DESeqDataSetFromMatrix(round(data2,digits=0), merged_new2, ~ cluster + subtype)
#atacDDS <- DESeqDataSetFromMatrix(round(data2,digits=0), merged_new2, ~ cluster )
#atacDDS_deseq <- DESeq(atacDDS,fitType="local")
#atacDDS_deseq <- DESeq(atacDDS,minReplicatesForReplace = 11,parallel=T)
atacDDS_deseq <- DESeq(atacDDS,minReplicatesForReplace = Inf,parallel=T)
#atac_Rlog <- rlog(atacDDS_deseq,fitType = "local")
result1 <- results(atacDDS_deseq, contrast=c('cluster','2','3'))
write.table(result1,file="cluster2_vs1_3.xls",quote=F,sep="\t")
save(atacDDS_deseq,file="atacDDS2_vs1_3.RData")
#save(atac_Rlog,file="atac_Rlog_mean.RData")

#result1
maxCooks<-apply(assays(atacDDS_deseq)[["cooks"]], 1, max)
result11<-cbind(result1,maxCooks)
result2<-result11[which(result11$maxCooks < 100),]
resSig <- result2[which(result2$padj < 0.01 & result2$pvalue < 0.01), ]
#filter 20
resSig2 <- resSig[which(resSig$log2FoldChange > 1),]
resSig2 <- resSig2[which(resSig2$baseMean > 20 & resSig2$lfcSE < 1),]
write.table(resSig2,file="cluster2_vs1_3_mean_20.xls",quote=F,sep="\t")
#filter 30
resSig2 <- resSig[which(resSig$log2FoldChange > 1),]
resSig2 <- resSig2[which(resSig2$baseMean > 30 & resSig2$lfcSE < 1),]
write.table(resSig2,file="cluster2_vs1_3_mean_30.xls",quote=F,sep="\t")
#filter 40
resSig2 <- resSig[which(resSig$log2FoldChange > 1),]
resSig2 <- resSig2[which(resSig2$baseMean > 40 & resSig2$lfcSE < 1),]
write.table(resSig2,file="cluster2_vs1_3_mean_40.xls",quote=F,sep="\t")

#filter_sd
#ressig3<-resSig2[which(as.numeric(counts_norm_sd) < 60 & as.numeric(counts_norm_sd2) <60  & as.numeric(counts_norm_sd3) < 60),]
#write.table(ressig3,file="/disk1/xilu/glioblastoma/cell_line_51/compared_1_2_3-c/cluster2_vs1_3.sd.xls",quote=F,sep="\t")
