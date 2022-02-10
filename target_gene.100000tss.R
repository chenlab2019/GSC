args<-commandArgs(TRUE)
suppressMessages(library(rGREAT))
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(ChIPseeker)
library(org.Hs.eg.db)
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
anno_LiverMinusHindbrain <- annotatePeak(cluster_c_g,TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene)
anno_LiverMinusHindbrain_GRanges <- as.GRanges(anno_LiverMinusHindbrain)
a<-getSYMBOL(anno_LiverMinusHindbrain_GRanges$geneId, data='org.Hs.eg')
anno_LiverMinusHindbrain_GRanges$genes<-a
b<-anno_LiverMinusHindbrain_GRanges[which(anno_LiverMinusHindbrain_GRanges$distanceToTSS < 500000 & anno_LiverMinusHindbrain_GRanges$distanceToTSS > -500000),]
write.table(data.frame(geneid=b$geneId,genes=b$genes),sep="\t",quote=F,file=output_gene,row.names=F)
write.table(as.data.frame(b),sep="\t",quote=F,file=output_anno)
write.table(anno_LiverMinusHindbrain,file=paste(input2,"annotation.txt",sep=""),quote=F,sep="\t")
pdf_name=paste(input2,".pdf",sep="")
pdf(pdf_name)
plotAnnoPie(anno_LiverMinusHindbrain)
dev.off()

#pdf("cluster1_2.e.uniq.pdf")
#plotAnnoPie(MacsCalls_chr20_filteredAnno)
#dev.off()

library(clusterProfiler)
kk2<-enrichKEGG(b$geneId,organism='hsa', pvalueCutoff = 0.05)
kk3<-setReadable(kk2,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
kk4 <- setReadable(kk2, 'org.Hs.eg.db', 'ENTREZID')
write.table(kk4,file=output_kegg,quote=F,sep="\t")
pdf(output_kegg_pdf)
dotplot(kk4, showCategory=30)
dev.off()
pdf(output_kegg_pdf2)
#heatplot(kk4)
dev.off()
#library(enrichplot)
#emapplot(kk2, showCategory = 36)
#cnetplot(kk3,showCategory =c("TGF-beta signaling pathway","PI3K-Akt signaling pathway","EGFR tyrosine kinase inhibitor resistance","ErbB signaling pathway","Ras signaling pathway"),categorySize=7)

#GSEA
#anno_LiverMinusHindbrain_GRanges$LFC<-cluster1$V3
#genelist<-b$LFC
#names(genelist)<-b$geneId
#genelist2<-genelist[order(genelist,decreasing=T)]
#y <- gsePathway(genelist2, nPerm=10000,pvalueCutoff=0.2,pAdjustMethod="BH",verbose=FALSE)
#pdf("gsea.pdf")
#gseaplot(y, geneSetID = "R-HSA-162582")
#dev.off()

#reactomePA
library(ReactomePA)
x <- enrichPathway(b$geneId,pvalueCutoff=0.05, readable=T)
edox <- setReadable(x, 'org.Hs.eg.db', 'ENTREZID')
write.table(x,file=output_rea,quote=F,sep="\t")
pdf(output_rea_pdf)
dotplot(x, showCategory=30)
dev.off()
pdf(output_rea_pdf2)
#heatplot(edox)
dev.off()
#pdf("rea.pdf")
#emapplot(x,showCategory = 18)
#dev.off()
#kk3<-setReadable(x,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
#write.table(as.data.frame(kk3),sep="\t",quote=F,file="cluster1.reactome.anno.genes.xls")
#pdf("rea.net.pdf")
#cnetplot(kk3,showCategory =c("Neuronal System","Signalling by NGF","Cell-Cell communication"),categorySize=3)
#dev.off()


######
#library(pathview)
#names(logFC) <- anno_LiverMinusHindbrain_GRanges$geneId
#logFC <- cluster1$V3
#pathview(gene.data = logFC,pathway.id = "hsa04550",species="hsa")

