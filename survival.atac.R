library("survminer")
require("survival")

args<-commandArgs(T)
input<-args[1]

data<-read.table(input,header=F)
colnames(data)<-c("cellline","cluster","survival_days")
data$survival_days<-as.numeric(as.character(data$survival_days))
fit <- survfit(Surv(survival_days) ~ cluster, data = data)

names<-paste(input,".all.pdf")
a<-surv_pvalue(fit,data)
b<-a$pval
print (c(names,b))

pdf(names)
ggsurvplot(fit, data = data,pval = TRUE,break.y.by = 0.2)
dev.off()

data2<-data[c(which(data$cluster =="1"),which(data$cluster =="2")),]
fit <- survfit(Surv(survival_days) ~ cluster, data = data2)
names2<-paste(input,".1vs2.pdf")
a<-surv_pvalue(fit)
b<-a$pval
print (c(names2,b))

pdf(names2)
ggsurvplot(fit, data = data2,pval = TRUE)
dev.off()

data2<-data[c(which(data$cluster =="2"),which(data$cluster =="3")),]
fit <- survfit(Surv(survival_days) ~ cluster, data = data2)

names3<-paste(input,".2vs3.pdf") 
a<-surv_pvalue(fit)
b<-a$pval
print (c(names3,b))

pdf(names3)
ggsurvplot(fit, data = data2,pval = TRUE)
dev.off()

data2<-data[c(which(data$cluster =="1"),which(data$cluster =="3")),]
fit <- survfit(Surv(survival_days) ~ cluster, data = data2)
names4<-paste(input,".1vs3.pdf")
a<-surv_pvalue(fit)
b<-a$pval
print (c(names4,b))

pdf(names4)
ggsurvplot(fit, data = data2,pval = TRUE)
dev.off()
