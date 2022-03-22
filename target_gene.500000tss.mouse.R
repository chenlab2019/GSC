args<-commandArgs(TRUE)
suppressMessages(library(rGREAT))
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library(clusterProfiler)
library(ChIPseeker)
library(org.Mm.eg.db)
library(annotate)

input2<-args[1]
input<-read.table(input2,header=F)
title<-unlist(strsplit(input2,"%",fixed=TRUE))
output_gene<-paste(title[1],"500000.gene.xls",sep=".")
output_anno<-paste(title[1],"500000.anno.genes.xls",sep=".")
output_kegg<-paste("kegg",title[1],"500000.pathway.xls",sep=".")
output_kegg_pdf<-paste("kegg",title[1],"500000.pathway.pdf",sep=".")
output_kegg_pdf2<-paste("kegg",title[1],"500000.pathway.heatplot.pdf",sep=".")
output_rea<-paste("reatomePA",title[1],"500000.pathway.xls",sep=".")
output_rea_pdf<-paste("reatomePA",title[1],"500000.pathway.pdf",sep=".")
output_rea_pdf2<-paste("reatomePA",title[1],"500000.pathway.heatplot.pdf",sep=".")
output_propotion<-paste(title[1],"500000.propotion.xls",sep=".")
cluster2_anno<-data.frame(seqnames=input[,1],start=input[,2],end=input[,3],width=input[,3]-input[,2],strand=rep("*",length(input[,2])))
rownames(cluster2_anno)<-paste(input[,1],input[,2],input[,3],sep="_")
cluster_c<-cluster2_anno
cluster_c_g<-makeGRangesFromDataFrame(cluster_c, keep.extra.columns = TRUE)
anno_LiverMinusHindbrain <- annotatePeak(cluster_c_g,TxDb = TxDb.Mmusculus.UCSC.mm9.knownGene)
anno_LiverMinusHindbrain_GRanges <- as.GRanges(anno_LiverMinusHindbrain)
a<-getSYMBOL(anno_LiverMinusHindbrain_GRanges$geneId, data='org.Mm.eg')
anno_LiverMinusHindbrain_GRanges$genes<-a
b<-anno_LiverMinusHindbrain_GRanges[which(anno_LiverMinusHindbrain_GRanges$distanceToTSS < 500000 & anno_LiverMinusHindbrain_GRanges$distanceToTSS > -500000),]
write.table(data.frame(geneid=b$geneId,genes=b$genes),sep="\t",quote=F,file=output_gene,row.names=F)
write.table(as.data.frame(b),sep="\t",quote=F,file=output_anno)


library(clusterProfiler)
kk2<-enrichKEGG(b$geneId,organism='mmu', pvalueCutoff = 0.05)
kk3<-setReadable(kk2,OrgDb = org.Mm.eg.db,keyType = "ENTREZID")
write.table(kk2,file=output_kegg,quote=F,sep="\t")
pdf(output_kegg_pdf)
dotplot(kk2, showCategory=20)

#reactomePA
library(ReactomePA)
x <- enrichPathway(b$geneId,organism='mouse',pvalueCutoff=0.05, readable=T)
write.table(x,file=output_rea,quote=F,sep="\t")

