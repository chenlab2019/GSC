args<-commandArgs(T)
input<-args[1]
input2<-args[2]

data<-readRDS(input)
samples<-read.table(input2,header=T,sep="\t")
sample_2<-data.frame(duplicates=samples[,2],comparison=samples[,3])
rownames(sample_2)<-samples[,1]
sample_2<-unique(sample_2)
w<-NULL
for ( i in colnames(data)){
    a<-which(sample_2$duplicates == i)
    w<-c(w,a)
}
sample_4<-sample_2[w,]
sample_3<-data.frame(group=sample_4[,2])
rownames(sample_3)<-sample_4[,1]

library(DESeq2)
atacDDS <- DESeqDataSetFromMatrix(round(data,digits=0), sample_3, ~ group)
#atacDDS <- collapseReplicates ( atacDDS2, groupby = atacDDS2$collapse)
#atacDDS <- DESeqDataSetFromMatrix(round(data2,digits=0), merged_new2, ~ cluster )
#atacDDS_deseq <- DESeq(atacDDS,fitType="local")
#atacDDS_deseq <- DESeq(atacDDS,minReplicatesForReplace = 11,parallel=T)
atacDDS_deseq <- DESeq(atacDDS,minReplicatesForReplace = Inf,parallel=T)
vsd <- vst(atacDDS_deseq, blind=FALSE)
matrix<-assay(vsd,blind=FALSE)

#atac_Rlog <- rlog(atacDDS_deseq,fitType = "local")
result1 <- results(atacDDS_deseq, contrast=c('group','cluster1','others'))

#filter 20 
#resSig <- result1[which(result1$padj < 0.05 & result1$pvalue < 0.05), ]
resSig <- result1[which( result1$padj < 0.05), ]
resSig2 <- resSig[which(resSig$log2FoldChange < -1),]
resSig2 <- resSig2[which(resSig2$baseMean > 20 & resSig2$lfcSE < 1),]
write.table(resSig2,file="DE_control_vs_NDST_DE.low.xls",quote=F,sep="\t")
write.table(result1,file="DE_control_vs_NDST_DE.all.xls",quote=F,sep="\t")
#filter 

