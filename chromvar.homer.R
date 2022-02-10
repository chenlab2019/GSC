library(TFBSTools)
library(chromVAR)
library(motifmatchr)
library(SummarizedExperiment)
library(Matrix)
library(ggplot2)
library(BiocParallel)
library(pheatmap)
library("gplots")
library("BSgenome.Hsapiens.UCSC.hg19")
library(chromVARmotifs)

#homer_database_from_homer_software(414)
database_homer<-read.table("/disk1/xilu/database/homer/homer.chromvar.motif.R3",header=T,sep="\t")
database_homer<-as.matrix(database_homer)
pwm_homer<-list()
for (i in seq(1,dim(database_homer)[[1]],1)){
    matrix_advanced<-database_homer[i,7]
    num<-as.numeric(unlist(strsplit(as.character(matrix_advanced),",")))
    num_matrix<-matrix(num,byrow=F,nrow=4,dimnames=list(c("A","C","G","T")))
    num_matrix_1<-apply(num_matrix,2,function(x){x/sum(x)})
    matrix<-log(num_matrix_1/0.25)
    tempt<-PWMatrix(ID=as.character(database_homer[i,1]),name=database_homer[i,2],matrixClass=database_homer[i,3],strand=database_homer[i,4],pseudocounts=numeric(),bg=c(A=0.25, C=0.25, G=0.25, T=0.25),tag=list(evidence=as.character(database_homer[i,5]),source=as.character(database_homer[i,6])),profileMatrix=matrix(matrix,byrow=F,nrow=4,dimnames=list(c("A","C","G","T"))))
#    assign(database_homer[i,1],tempt)
    pwm_homer[[database_homer[i,1]]]<-tempt
}
motifs <- do.call(PWMatrixList,pwm_homer)

set.seed(1234)
#chromvar analysis 
df2<-readRDS("/disk1/xilu/glioblastoma/cell_line_140_new_pipeline_call_peaks/differential_peaks/deseq2.data.rds")
df<-df2[,-which(colnames(df2) %in% c("U3175"))]
cluster.s<-read.table("cluster",header=F,sep="\t")
#peak_chr<-paste(df$chr,df$start,df$end,sep="_")
df4<-df[which(rownames(df) %in% cluster.s$V1),]

cluster<-read.table("/disk1/xilu/glioblastoma/corrected_order_50/new.cluster.50",header=F,sep="\t")
cluster<-data.frame(overlap_cellline=cluster$V1,subtype_latest=cluster$V2)
cluster$subtype_latest<-c(rep("cluster1",19),rep("cluster2",14),rep("cluster3",17))
colData<-DataFrame(type=cluster$subtype_latest,row.names=cluster$overlap_cellline)

names<-matrix(unlist(strsplit(rownames(df4),"_")),byrow=T,ncol=3)
names<-as.data.frame(names)
df3<-cbind(names,df4[,1:50])
colnames(df3)[1]<-"chr"
colnames(df3)[2]<-"start"
colnames(df3)[3]<-"end"


doc<-NULL
for (i in colnames(as.matrix(df3[,4:53]))){
    num<-which(rownames(colData) == i )
    doc<-c(doc,num)
    }

colData2<-DataFrame(type=colData$type[doc],depth=as.numeric(colSums(df3[,4:53])),row.names=rownames(colData)[doc])

rowRanges <- GRanges(df3$chr,IRanges(start=as.numeric(as.character(df3$start)),end=as.numeric(as.character(df3$end))),strand=rep("*",length(df3$chr)),seqlengths=c(chr1=249250621,chr2=243199373,chr3=198022430,chr4=191154276,chr5=180915260,chr6=171115067,chr7=159138663,chrX=155270560,chr8=146364022,chr9=141213431,chr10=135534747,chr11=135006516,chr12=133851895,chr13=115169878,chr14=107349540,chr15=102531392,chr16=90354753,chr17=81195210,chr18=78077248,chr20=63025520,chrY=59373566,chr19=59128983,chr22=51304566,chr21=48129895,chr6_ssto_hap7=4928567,chr6_mcf_hap5=4833398,chr6_cox_hap2=4795371,chr6_mann_hap4=4683263,chr6_apd_hap1=4622290,chr6_qbl_hap6=4611984,chr6_dbb_hap3=4610396,chr17_ctg5_hap1=1680828,chr4_ctg9_hap1=590426,chr1_gl000192_random=547496,chrUn_gl000225=211173,chr4_gl000194_random=191469,chr4_gl000193_random=189789,chr9_gl000200_random=187035,chrUn_gl000222=186861,chrUn_gl000212=186858,chr7_gl000195_random=182896,chrUn_gl000223=180455,chrUn_gl000224=179693,chrUn_gl000219=179198,chr17_gl000205_random=174588,chrUn_gl000215=172545,chrUn_gl000216=172294,chrUn_gl000217=172149,chr9_gl000199_random=169874,chrUn_gl000211=166566,chrUn_gl000213=164239,chrUn_gl000220=161802,chrUn_gl000218=161147,chr19_gl000209_random=159169,chrUn_gl000221=155397,chrUn_gl000214=137718,chrUn_gl000228=129120,chrUn_gl000227=128374,chr1_gl000191_random=106433,chr19_gl000208_random=92689,chr9_gl000198_random=90085,chr17_gl000204_random=81310,chrUn_gl000233=45941,chrUn_gl000237=45867,chrUn_gl000230=43691,chrUn_gl000242=43523,chrUn_gl000243=43341,chrUn_gl000241=42152,chrUn_gl000236=41934,chrUn_gl000240=41933,chr17_gl000206_random=41001,chrUn_gl000232=40652,chrUn_gl000234=40531,chr11_gl000202_random=40103,chrUn_gl000238=39939,chrUn_gl000244=39929,chrUn_gl000248=39786,chr8_gl000196_random=38914,chrUn_gl000249=38502,chrUn_gl000246=38154,chr17_gl000203_random=37498,chr8_gl000197_random=37175,chrUn_gl000245=36651,chrUn_gl000247=36422,chr9_gl000201_random=36148,chrUn_gl000235=34474,chrUn_gl000239=33824,chr21_gl000210_random=27682,chrUn_gl000231=27386,chrUn_gl000229=19913,chrM=16571,chrUn_gl000226=15008,chr18_gl000207_random=4262))

rse <- SummarizedExperiment(assays=SimpleList(counts=as.matrix(df3[,4:53])),rowRanges=rowRanges, colData=colData2)
example_counts <- addGCBias(rse,genome =BSgenome.Hsapiens.UCSC.hg19)
counts_filtered <- filterSamples(example_counts, min_depth = 1500,min_in_peaks = 0.15,shiny=FALSE)
counts_filtered <- filterPeaks(counts_filtered)

motif_ix <- matchMotifs(motifs, counts_filtered,genome = BSgenome.Hsapiens.UCSC.hg19)
dev <- computeDeviations(object = counts_filtered, annotations = motif_ix)

#tsne
#tsne_results <- deviationsTsne(dev, threshold = 1.5, perplexity = 11,max_iter = 20000,shiny=FALSE,theta=0.6)
tsne_results <- deviationsTsne(dev, threshold = 1, perplexity = 15,max_iter = 10000,shiny=FALSE,theta=0.6) 
#tsne_results <- deviationsTsne(dev, threshold = 1.5, perplexity = 10,max_iter = 5000,shiny=FALSE) 
#tsne_plots <- plotDeviationsTsne(dev, tsne_results,annotation_name =c("FOSL2","JUNB","JUND"),sample_column ="type",shiny = FALSE)
#tsne_plots <- plotDeviationsTsne(dev, tsne_results,annotation_name =c("NeuroG2(bHLH)","Olig2(bHLH)","GABPA(ETS)","Slug(Zf)"),sample_column ="type",shiny = FALSE)
tsne_plots <- plotDeviationsTsne(dev, tsne_results,annotation_name =c("RUNX2(Runt)","Ascl1(bHLH)","Olig2(bHLH)","Sox3(HMG)","JunB(bZIP)","Fosl2(bZIP)","Tcf3(HMG)","TCFL2(HMG)"),sample_column ="type",shiny = FALSE)
pdf("tsne_samples.pdf")
tsne_plots[[1]]
tsne_plots[[2]]
tsne_plots[[3]]
tsne_plots[[4]]
dev.off()

diff_acc <- differentialDeviations(dev, "type")
variability <- computeVariability(dev)
#filter<-deviations(dev)[intersect(which(diff_acc$p_value_adjusted < 0.05),which(variability$variability >= 1.5)),]
#filter<-assays(dev)$z[intersect(which(diff_acc$p_value_adjusted < 0.05),which(variability$variability >= 1.5)),]
#filter<-deviations(dev)[which(variability$variability >= 1),]
#filter<-assays(dev)$z[intersect(which(variability$p_value_adj < 0.05),which(variability$variability >= 1.5)),]

#top50 variability 
v<-variability[order(variability$variability,decreasing=T),][1:50,]
write.table(v,file="top50.TF.all.rank.xls",sep="\t",quote=F)

#heatmap
filter<-deviations(dev)[order(variability$variability,decreasing=T),][1:57,]
filter<-filter[c(-16,-45,-3,-46,-1,-15),]
merged_new2<-colData2
mergerd_new3<-merged_new2[order(gsub("cluster","",colData2$type),decreasing=F),]
mergerd_new3$type <- factor(mergerd_new3$type)
target2<-filter[,order(gsub("cluster","",colData2$type),decreasing=F)]

pdf("heatmap_deviation.pdf")
col<-as.character(mergerd_new3$type)
col2<-col
col2[grep("cluster1",col)]<-'#FFC0CB'
col2[grep("cluster2",col)]<-'#808080'
col2[grep("cluster3",col)]<-'#0000FF'
heatmap.2(target2,trace="none",scale="none",density="none",col=bluered(20),margins=c(10,7),labCol = FALSE,cexRow=0.5,ColSideColors=col2,Colv=FALSE)
heatmap.2(target2,trace="none",scale="row",density="none",col=bluered(20),margins=c(10,7),labCol = FALSE,cexRow=0.1,ColSideColors=col2,Colv=T)
legend("topright",legend = unique(col),col = unique(col2),lty= 1, lwd = 5,cex=.7)
dev.off()

cluster3.special=rownames(target2[which(rowMeans(target2[,34:50]) > 0 & rowMeans(target2[,20:34]) <0 & rowMeans(target2[,1:19]) <0  ),])
cluster2.special<-rownames(target2[which(rowMeans(target2[,20:33]) > 0 & rowMeans(target2[,1:19]) <0 & rowMeans(target2[,34:50]) <0 ),])
cluster1.special<-rownames(target2[which(rowMeans(target2[,1:19]) > 0 & rowMeans(target2[,20:33]) <0 & rowMeans(target2[,34:50]) <0),])
write.table(cluster3.special,file="cluster3.special.txt",quote=F)
write.table(cluster2.special,file="cluster2.special.txt",quote=F)
write.table(cluster1.special,file="cluster1.special.txt",quote=F)

pdf("heatmap_deviation.new.dfp.pdf")
pdf("heatmap_deviation.new.dfp.pdf",width = 10 , height = 30)
filter2<-deviations(dev)[order(variability$variability,decreasing=T),]
filter<-deviations(dev)[which(variability$variability >= 1),]
#grep("NRF\\(NRF\\)",rownames(deviations(dev)[order(variability$variability,decreasing=T),]))
a<-filter2[-11,]
filter<-a[1:50,]
merged_new2<-data.frame(cluster=factor(colData2[,1]))
rownames(merged_new2)<-rownames(colData2)
mergerd_new3<-data.frame(cluster=merged_new2[order(colData2$type),])
rownames(mergerd_new3)<-rownames(merged_new2)[order(colData2$type)]
mergerd_new3$cluster<-factor(mergerd_new3$cluster)
target2<-filter[,order(colData2$type)]
tf_name<-matrix(unlist(strsplit(rownames(target2),"/")),ncol=3, byrow = T)
tf_final<-matrix(unlist(strsplit(tf_name[,1],"\\(")),ncol=2, byrow = T)
rownames(target2)<-tf_final[,1]
name<-data.frame(rownames(target2))
f <- apply(name,2,toupper)
rownames(target2)<-f[,1]

library("RColorBrewer")
color = colorRampPalette(brewer.pal(3,"RdBu"))(21)
breaks<-seq(-1.5,1.5,by=0.15)
pheatmap(target2,cluster_rows=F,cluster_cols = F,border_color = "grey60",clustering_distance_cols="euclidean",scale="row",show_colnames= F,show_rownames = T,color = rev(color),annotation_col=mergerd_new3,clustering_distance_rows="euclidean", cex=1, clustering_method="complete",fontsize_col=4,breaks=breaks,annotation_colors = list(cluster=c("cluster1"="#ae82b7","cluster2"="#b1b376","cluster3"="#8d96a5")),cellheight=10,cellwidth=10,fontsize = 8,fontsize_row =10 )
pheatmap(target2,cluster_rows=T,cluster_cols = F,border_color = "grey60",clustering_distance_cols="euclidean",scale="row",show_colnames= F,show_rownames = T,color = rev(color),annotation_col=mergerd_new3,clustering_distance_rows="euclidean", cex=1, clustering_method="complete",fontsize_col=4,breaks=breaks,annotation_colors = list(cluster=c("cluster1"="#ae82b7","cluster2"="#b1b376","cluster3"="#8d96a5")),cellheight=10,cellwidth=10,fontsize = 8,fontsize_row =10 )
dev.off()

#self
ix <- intersect(which(variability$p_value_adj <0.05),which(variability$variability >= 1.5))
mat <- deviations(dev)[ix, , drop = FALSE]
tsne_res <- Rtsne::Rtsne(t(mat), perplexity = 16, max_iter = 5000,theta = 0.5)
out <- tsne_res$Y
rownames(out) <- colnames(dev)
tsne_results <-out
tsne_plots <- plotDeviationsTsne(dev, tsne_results,sample_column ="type",shiny = FALSE)
pdf("tsne_samples.test3.pdf")
tsne_plots[[1]]
dev.off()

