args<-commandArgs(TRUE)
suppressMessages(library(rGREAT))
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(clusterProfiler)
library(ChIPseeker)
library(org.Hs.eg.db)
library(annotate)

input<-args[1]
title<-unlist(strsplit(input,".",fixed=TRUE))
output_gene<-paste(title[1],"500000.gene.xls",sep=".")
output_anno<-paste(title[1],"500000.anno.genes.xls",sep=".")
output_kegg<-paste("kegg",title[1],"500000.pathway.xls",sep=".")
output_rea<-paste("reatomePA",title[1],"500000.pathway.xls",sep=".")

cluster1<-read.table(input,header=F,sep="\t")
cluster1_anno<-read.table("/disk1/xilu/glioblastoma/cell_line_140_new_pipeline_call_peaks/differential_peaks/anno.list",header=FALSE,sep="\t")
cluster2_anno<-data.frame(seqnames=cluster1_anno[,2],start=cluster1_anno[,3],end=cluster1_anno[,4],width=cluster1_anno[,5],strand=rep("*",length(cluster1_anno[,2])))
rownames(cluster2_anno)<-cluster1_anno[,1]
cluster1_anno<-cluster2_anno
cluster_c<-cluster1_anno[which(rownames(cluster1_anno) %in% cluster1[,1]),]
cluster_c_g<-makeGRangesFromDataFrame(cluster_c, keep.extra.columns = TRUE)
anno_LiverMinusHindbrain <- annotatePeak(cluster_c_g,TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene)
anno_LiverMinusHindbrain_GRanges <- as.GRanges(anno_LiverMinusHindbrain)
a<-getSYMBOL(anno_LiverMinusHindbrain_GRanges$geneId, data='org.Hs.eg')
anno_LiverMinusHindbrain_GRanges$genes<-a
b<-anno_LiverMinusHindbrain_GRanges[which(anno_LiverMinusHindbrain_GRanges$distanceToTSS < 500000 & anno_LiverMinusHindbrain_GRanges$distanceToTSS > -500000),]
write.table(data.frame(geneid=b$geneId,genes=b$genes),sep="\t",quote=F,file=output_gene,row.names=F)
write.table(as.data.frame(b),sep="\t",quote=F,file=output_anno)


library(clusterProfiler)
kk2<-enrichKEGG(b$geneId,organism='hsa', pvalueCutoff = 0.05)
kk3<-setReadable(kk2,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
write.table(kk2,file=output_kegg,quote=F,sep="\t")

#reactomePA
library(ReactomePA)
x <- enrichPathway(b$geneId,pvalueCutoff=0.05, readable=T)
write.table(x,file=output_rea,quote=F,sep="\t")

