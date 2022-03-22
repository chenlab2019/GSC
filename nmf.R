library(NMF)
library(matrixStats)
args<-commandArgs(T)
input<-args[1]

exprs_with_time2<-readRDS(input)
a<-rowVars(exprs_with_time2)
exprs_with_time=exprs_with_time2[rev(order(a))[1:70000],]

############rank = 2
res2_ATAC_allpeak <- nmf(exprs_with_time, 2,'ns',nrun=50,.opt='vP60',seed=122)
pdf(file ="cpm51samles_nrun50_nmf2.pdf",width = 8,height = 8)
consensusmap(res2_ATAC_allpeak)
#consensusmap(estim.r, color = "-RdYlBu", distfun = function(x) as.dist(1 - x), hclustfun = "complete", Rowv = TRUE, Colv = TRUE, info = FALSE,labCol=colnames(estim.r),labRow=colnames(estim.r))
dev.off()
#summary(res2_ATAC_allpeak)
cophcor(consensus(res2_ATAC_allpeak))
dispersion(consensus(res2_ATAC_allpeak))
si <- silhouette(res2_ATAC_allpeak)
summary(si)$avg.width

###############  rank = 3
res3_ATAC_allpeak <- nmf(exprs_with_time, 3,'ns',nrun=50,.opt='vP60',seed=122)
pdf(file ="cpm51samles_nrun50_nmf3.pdf",width = 8,height = 8)
consensusmap(res3_ATAC_allpeak)
dev.off()
#summary(res3_ATAC_allpeak)
cophcor(consensus(res3_ATAC_allpeak))
dispersion(consensus(res3_ATAC_allpeak))
si <- silhouette(res3_ATAC_allpeak)
summary(si)$avg.width

###############  rank = 4
res4_ATAC_allpeak <- nmf(exprs_with_time, 4,'ns',nrun=50,.opt='vP60', seed=122)
pdf(file ="cpm51samles_nrun50_nmf4.pdf",width = 8,height = 8)
consensusmap(res4_ATAC_allpeak)
dev.off()
#summary(res4_ATAC_allpeak)
cophcor(consensus(res4_ATAC_allpeak))
dispersion(consensus(res4_ATAC_allpeak))
si <- silhouette(res4_ATAC_allpeak)
summary(si)$avg.width

###############  rank = 5
res5_ATAC_allpeak <- nmf(exprs_with_time, 5,'ns',nrun=50,.opt='vP60', seed=122)
pdf(file ="cpm51samles_nrun50_nmf5.pdf",width = 8,height = 8)
consensusmap(res5_ATAC_allpeak)
dev.off()
#summary(res5_ATAC_allpeak)
cophcor(consensus(res5_ATAC_allpeak))
dispersion(consensus(res5_ATAC_allpeak))
si <- silhouette(res5_ATAC_allpeak)
summary(si)$avg.width


###############  rank = 6
res5_ATAC_allpeak <- nmf(exprs_with_time, 6,'ns',nrun=50,.opt='vP60', seed=122)
pdf(file ="cpm51samles_nrun50_nmf5.pdf",width = 8,height = 8)
consensusmap(res5_ATAC_allpeak)
dev.off()
#summary(res5_ATAC_allpeak)
cophcor(consensus(res5_ATAC_allpeak))
dispersion(consensus(res5_ATAC_allpeak))
si <- silhouette(res5_ATAC_allpeak)
summary(si)$avg.width

###############  rank = 7
res5_ATAC_allpeak <- nmf(exprs_with_time, 7,'ns',nrun=50,.opt='vP60', seed=122)
pdf(file ="cpm51samles_nrun50_nmf5.pdf",width = 8,height = 8)
consensusmap(res5_ATAC_allpeak)
dev.off()
#summary(res5_ATAC_allpeak)
cophcor(consensus(res5_ATAC_allpeak))
dispersion(consensus(res5_ATAC_allpeak))
si <- silhouette(res5_ATAC_allpeak)
summary(si)$avg.width

###############  rank = 8
res5_ATAC_allpeak <- nmf(exprs_with_time, 8,'ns',nrun=50,.opt='vP60', seed=122)
pdf(file ="cpm51samles_nrun50_nmf5.pdf",width = 8,height = 8)
consensusmap(res5_ATAC_allpeak)
dev.off()
#summary(res5_ATAC_allpeak)
cophcor(consensus(res5_ATAC_allpeak))
dispersion(consensus(res5_ATAC_allpeak))
si <- silhouette(res5_ATAC_allpeak)
summary(si)$avg.width
