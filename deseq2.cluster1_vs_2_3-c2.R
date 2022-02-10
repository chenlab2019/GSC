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
merged_new2$cluster[grep("2",merged_new2$cluster)]<-3
merged_new2$cluster<-factor(merged_new2$cluster)

library(DESeq2)
#atacDDS <- DESeqDataSetFromMatrix(round(data2,digits=0), merged_new2, ~ cluster + subtype)
atacDDS <- DESeqDataSetFromMatrix(round(data2,digits=0), merged_new2, ~ cluster)
#atacDDS_deseq <- DESeq(atacDDS,minReplicatesForReplace = 51,parallel=T)
atacDDS_deseq <- DESeq(atacDDS,minReplicatesForReplace = Inf,parallel=T)
#atacDDS_deseq <- DESeq(atacDDS,fitType="local")

#atac_Rlog <- rlog(atacDDS_deseq,fitType = "local")
result1 <- results(atacDDS_deseq, contrast=c('cluster','1','3'))
write.table(result1,file="cluster1_vs2_3-2.xls",quote=F,sep="\t")
save(atacDDS_deseq,file="atacDDS1_vs2_3.RData")
#save(atac_Rlog,file="atac_Rlog_mean.RData")

#counts<-assays(atacDDS_deseq)$counts
#result1
cluster1<-which(merged_new2$cluster == 1)
maxCooks<-apply(assays(atacDDS_deseq)[["cooks"]], 1, max) # qf(0.7,1,50)
#cluster1_means<-rowMeans(normalized[,cluster1])
result11<-cbind(result1,maxCooks)
result2<-result11[which(result11$maxCooks < 100),]
#result2 <- result1[as.numeric(which(maxCooks < 2.4)),]
resSig <- result2[which(result2$padj < 0.01 & result2$pvalue < 0.01), ]
#filter 20
resSig2 <- resSig[which(resSig$log2FoldChange > 1),]
resSig2 <- resSig2[which(resSig2$baseMean > 20 & resSig2$lfcSE < 1),]
write.table(resSig2,file="cluster1_vs2_3_mean_20.xls",quote=F,sep="\t")
#filter 30
resSig2 <- resSig[which(resSig$log2FoldChange > 1),]
resSig2 <- resSig2[which(resSig2$baseMean > 30 & resSig2$lfcSE < 1),]
write.table(resSig2,file="cluster1_vs2_3_mean_30.xls",quote=F,sep="\t")
#filter 40 
resSig2 <- resSig[which(resSig$log2FoldChange > 1),]
resSig2 <- resSig2[which(resSig2$baseMean > 40 & resSig2$lfcSE < 1),] 
write.table(resSig2,file="cluster1_vs2_3_mean_40.xls",quote=F,sep="\t")
#counts filter
#counts_filter<-counts[rownames(resSig2),]
#ressig3<-resSig2[which(as.numeric(counts_norm_sd) < 50 & as.numeric(counts_norm_sd2) < 50 & as.numeric(counts_norm_sd3) < 50),]
#write.table(ressig3,file="compared_1_2_3-c/cluster1_vs2_3.sd.xls",sep="\t",quote=F)

#mergerd_new3<-merged_new2[order(merged_new2$cluster,decreasing=F),]
#a<-counts(atacDDS_deseq,normalized=T)
#a1<-a[,order(merged_new2$cluster,decreasing=F)]
#counts_norm<-a1[rownames(resSig2),]
#counts_norm_sd<-apply(counts_norm[,1:20],1,sd)
#counts_norm_sd2<-apply(counts_norm[,21:31],1,sd)
#counts_norm_sd3<-apply(counts_norm[,32:51],1,sd)
#ressig3<-resSig2[which(as.numeric(counts_norm_sd) < 50 & as.numeric(counts_norm_sd2) < 50 & as.numeric(counts_norm_sd3) < 50),]


#cooks's distance 
#maxCooks<-apply(assays(atacDDS_deseq)[["cooks"]], 1, max)
#pdf("basemean_cooks.pdf")
#plot(mcols(atacDDS_deseq)$baseMean,maxCooks)
#dev.off()
#pdf("lfc_cooks.pdf")
#plot(result1$log2FoldChange,maxCooks)
#dev.off()
#cutoff
#qf(0.99,2,51-2)
#result1$pvalue[maxCooks > 2.4] <-NA
#result1$pvalue[res$baseMean < 16] <- NA
#result1$padj <- p.adjust(result1$pvalue, method="BH")

