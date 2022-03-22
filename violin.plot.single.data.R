library(ggplot2)
library("Useful2me")
args<-commandArgs(T)
input<-args[1]
input2<-args[2]

data<-read.csv(input,header=F)

library(data.table)
data$V5 <- as.numeric(data$V3 %like% "Promoter")
data2<-data[,c(1,2,4,5)]
colnames(data2)<-c("position","mean","gene","region")

vst<-read.table(input2,header=F,sep="\t")
colnames(vst)[1]<-"position"
a<-merge(vst,data2,by="position")
a2<-a[,c(1:52,54,55)]
#a2<-a[,c(1:51,53,54)]

#a3<-a2[-c(which(a2$gene == "CADM2-AS2"),which(a2$gene == "DBI")),]
#a2<-a3

list_new<-NULL
for (i  in unique(a2$gene)){
combine<-data.frame(accessibility=c(as.numeric(unlist(a2[which(a2$region == 1 & a2$gene == i),c(2:52)]))),region=rep("1",length(as.numeric(unlist(a2[which(a2$region == 1 & a2$gene == i),c(2:52)])))),gene=rep(i,length(as.numeric(unlist(a2[which(a2$region == 1 & a2$gene == i),c(2:52)])))))
combine2<-data.frame(accessibility=c(as.numeric(unlist(a2[which(a2$region == 0 & a2$gene == i),c(2:52)]))),region=rep("0",length(as.numeric(unlist(a2[which(a2$region == 0 & a2$gene == i),c(2:52)])))),gene=rep(i,length(as.numeric(unlist(a2[which(a2$region == 0 & a2$gene == i),c(2:52)])))))
new<-rbind(combine,combine2)
list_new<-rbind(list_new,new)
}

list_new2<-cbind(list_new,data.frame(name=rep(input,dim(list_new)[[1]])))
name<-paste("violin",input2,"graph.pdf",sep="")
name_table<-paste(input2,".table.out.xls",sep="")
write.table(list_new2,file=name_table,quote=F,sep="\t")

pdf(name)
ggplot(list_new, aes(x=gene, y=accessibility , fill=region)) +geom_split_violin()+theme_classic()+xlab("")+ylab("normalized chromatin accessibility")+scale_fill_manual(values=c("#9ecae1", "#3182bd"))+theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
dev.off()



