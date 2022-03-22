library(ggpubr)
library(Rmisc)
args<-commandArgs(T)
input<-args[1]

data<-read.table(input,header=F)
data$V3<-as.numeric(gsub(",",".",data$V3))
tgc2 <- summarySE(data, measurevar="V3", groupvars=c("V2"))
data$V2<-factor(data$V2,levels=c("PN","MS","CL"))
my_comparisons<-list(c("PN","MS"),c("PN","CL"),c("MS","CL"))


name<-paste(input,"p.format.pdf")
g<-ggplot(data,aes(x=V2,y=V3,shape=as.factor(V2),fill=as.factor(V2),color=as.factor(V2)))+geom_jitter(size=4,aes(colour=V2),width = 0.1, height = 0.2)+scale_color_manual(values=c("CL"="blue", "MS"="red", "PN"="purple"))+stat_compare_means(comparisons = my_comparisons,label="p.format",method="t.test",paired=F)+geom_errorbar(data = tgc2,fun.data="mean_se",aes(x=V2,y=V3,ymin=V3-sd,ymax=V3+sd,mean=V3),size=0.3,width=0.2)+stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,geom = "crossbar",width = 0.5,size=0.2)+theme_classic()
pdf(name)
g
dev.off()

#my_comparisons<-list(c("1","2"),c("1","3"),c("2","3"))
#ata<-read.csv2(input,header=F)
#data$V2<-as.character(data$V2)
#tgc2 <- summarySE(data, measurevar="V3", groupvars=c("V2"))
#g<-ggplot(data,aes(x=V2,y=V3,shape=V2,fill=V2))+geom_jitter(size=4,aes(colour=V2),width = 0.1, height = 0.2)+scale_color_manual(values=c("1"="blue", "2"="red", "3"="purple"))+stat_compare_means(comparisons = my_comparisons,label="p.signif",method="t.test",paired=F)+geom_errorbar(data = tgc2,fun.data="mean_se",aes(x=V2,y=V3,ymin=V3-sd,ymax=V3+sd,mean=V3),size=0.3,width=0.2)+stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,geom = "crossbar",width = 0.5,size=0.2)+theme_classic()
#g<-ggplot(data,aes(x=V2,y=V3,shape=V2,fill=V2))+geom_jitter(size=4,aes(colour=V2),width = 0.3, height = 0.2)+scale_color_manual(values=c("CL"="blue", "MS"="red", "PN"="purple"))+stat_compare_means(comparisons = my_comparisons,label="p.format",method="t.test",paired=F)+geom_errorbar(data = tgc2,fun.data="mean_se",aes(x=V2,y=V3,ymin=V3-1/2*sd,ymax=V3+1/2*sd,mean=V3),size=0.3,width=0.2)+stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,geom = "crossbar",width = 0.5,size=0.2)+scale_shape_manual(values=c(19,15,17))+ylim(-50,200)+theme_classic()
name2<-paste(input,"p.signif.pdf")
g<-ggplot(data,aes(x=V2,y=V3,shape=as.factor(V2),fill=as.factor(V2),color=as.factor(V2)))+geom_jitter(size=4,aes(colour=V2),width = 0.1, height = 0.2)+scale_color_manual(values=c("CL"="blue", "MS"="red", "PN"="purple"))+stat_compare_means(comparisons = my_comparisons,label="p.signif",method="t.test",paired=F)+geom_errorbar(data = tgc2,fun.data="mean_se",aes(x=V2,y=V3,ymin=V3-sd,ymax=V3+sd,mean=V3),size=0.3,width=0.2)+stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,geom = "crossbar",width = 0.5,size=0.2)+theme_classic()
pdf(name2)
g
dev.off()
