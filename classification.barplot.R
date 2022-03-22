suppressWarnings(library(ggplot2))
suppressWarnings(library(RColorBrewer))

input_file<-read.table("cluster.new.tcga",sep="\t",header=T)
cluster1<-input_file[which(input_file$new_cluster == "1"),]
cluster2<-input_file[which(input_file$new_cluster == "2"),]
cluster3<-input_file[which(input_file$new_cluster == "3"),]

cluster1_n<-data.frame(region=cluster1[,2],num=rep(1,dim(cluster1)[1]))
cluster1_accum<-aggregate(num~region,cluster1_n,sum)
cluster1_final<-data.frame(subtype=cluster1_accum$region,cluster=rep("cluster1(N=19)",dim(cluster1_accum)[1]),num=cluster1_accum$num/sum(cluster1_accum$num))

cluster2_n<-data.frame(region=cluster2[,2],num=rep(1,dim(cluster2)[1]))
cluster2_accum<-aggregate(num~region,cluster2_n,sum)
cluster2_final<-data.frame(subtype=cluster2_accum$region,cluster=rep("cluster2(N=14)",dim(cluster2_accum)[1]),num=cluster2_accum$num/sum(cluster2_accum$num))

cluster3_n<-data.frame(region=cluster3[,2],num=rep(1,dim(cluster3)[1]))
cluster3_accum<-aggregate(num~region,cluster3_n,sum)
cluster3_final<-data.frame(subtype=cluster3_accum$region,cluster=rep("cluster3(N=17)",dim(cluster3_accum)[1]),num=cluster3_accum$num/sum(cluster3_accum$num))

cluster<-rbind(cluster1_final,cluster2_final,cluster3_final)
cluster$num<-round(cluster$num,4)
fill<-c("#639EFC","#E2675E","#7676EE")
pdf("classification.barplot.pdf")
#ggplot(cluster, aes(x = cluster, y = num,fill = subtype))+geom_col()+scale_fill_brewer(palette="RdBu")+geom_text(aes(label=scales::percent(num)),position=position_stack(vjust=0.5))+ylab("percent")+theme_classic()
ggplot(cluster, aes(x = cluster, y = num,fill = subtype))+geom_col()+scale_fill_manual(values=fill)+geom_text(aes(label=scales::percent(num)),position=position_stack(vjust=0.5))+ylab("percent")+theme_classic()
dev.off()

