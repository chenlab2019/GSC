library(NMF)
library(ConsensusClusterPlus)

exprs_with_time<-readRDS("vst_51.rds")
#estim.r <- nmf(exprs_with_time, 1:8, nrun=30, seed=123456)
cophcor(estim.r)
#consensusmap(estim.r, color = "-RdYlBu", distfun = function(x) as.dist(1 - x), hclustfun = "complete", Rowv = TRUE, Colv = TRUE, info = FALSE,labCol=colnames(estim.r),labRow=colnames(estim.r))
#dev.off()
#save(estim.r,file="estim.r.RData")
#pdf("estim.pdf")
#plot(estim.r)
#dev.off()

###############  rank = 2
res2_ATAC_allpeak <- nmf(exprs_with_time, 2, nrun=40, seed=12)
pdf(file ="./consensusmap(res2_ATAC_65_allpeaks).pdf",width = 8,height = 8)
consensusmap(res2_ATAC_allpeak)
consensusmap(estim.r, color = "-RdYlBu", distfun = function(x) as.dist(1 - x), hclustfun = "complete", Rowv = TRUE, Colv = TRUE, info = FALSE,labCol=colnames(estim.r),labRow=colnames(estim.r))
dev.off()
#summary(res2_ATAC_allpeak)

###############  rank = 3
res3_ATAC_allpeak <- nmf(exprs_with_time, 3, nrun=40, seed=123)
pdf(file ="./consensusmap(res3_ATAC_65_allpeaks).pdf",width = 8,height = 8)
consensusmap(res3_ATAC_allpeak)
dev.off()
#summary(res3_ATAC_allpeak)

###############  rank = 4
res4_ATAC_allpeak <- nmf(exprs_with_time, 4, nrun=40, seed=1234)
pdf(file ="./consensusmap(res4_ATAC_65_allpeaks).pdf",width = 8,height = 8)
consensusmap(res4_ATAC_allpeak)
dev.off()
#summary(res4_ATAC_allpeak)

###############  rank = 5
res5_ATAC_allpeak <- nmf(exprs_with_time, 5, nrun=40, seed=12345)
pdf(file ="./consensusmap(res5_ATAC_65_allpeaks).pdf",width = 8,height = 8)
consensusmap(res5_ATAC_allpeak)
dev.off()
#summary(res5_ATAC_allpeak)

###############  rank = 6
res5_ATAC_allpeak <- nmf(exprs_with_time, 6, nrun=40, seed=12345)
pdf(file ="./consensusmap(res6_ATAC_65_allpeaks).pdf",width = 8,height = 8)
consensusmap(res5_ATAC_allpeak)
dev.off()
#summary(res5_ATAC_allpeak)

###############  rank = 7
res5_ATAC_allpeak <- nmf(exprs_with_time, 7, nrun=40, seed=12345)
pdf(file ="./consensusmap(res7_ATAC_65_allpeaks).pdf",width = 8,height = 8)
consensusmap(res5_ATAC_allpeak)
dev.off()
#summary(res5_ATAC_allpeak)

###############  rank = 8
res5_ATAC_allpeak <- nmf(exprs_with_time, 8, nrun=40, seed=12345)
pdf(file ="./consensusmap(res8_ATAC_65_allpeaks).pdf",width = 8,height = 8)
consensusmap(res5_ATAC_allpeak)
dev.off()
#summary(res5_ATAC_allpeak)