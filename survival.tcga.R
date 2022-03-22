library("survminer")
require("survival")

args<-commandArgs(T)
input<-args[1]

data<-read.table(input,header=F)
colnames(data)<-c("cellline","cluster","survival_days")
data$survival_days<-as.numeric(as.character(data$survival_days))
fit <- survfit(Surv(survival_days) ~ cluster, data = data)

a<-surv_pvalue(fit,data)
b<-a$pval

names<-paste(input,".all.pdf")
print (c(names,b))
pdf(names)
ggsurvplot(fit, data = data,pval = TRUE,break.y.by=0.2)
dev.off()

data2<-data[c(which(data$cluster =="PN"),which(data$cluster =="CL")),]
fit <- survfit(Surv(survival_days) ~ cluster, data = data2)
names2<-paste(input,".PNvsCL.pdf")

a<-surv_pvalue(fit,data2)
b<-a$pval
print (c(names2,b))

pdf(names2)
ggsurvplot(fit, data = data2,pval = TRUE)
dev.off()

data2<-data[c(which(data$cluster =="PN"),which(data$cluster =="MS")),]
fit <- survfit(Surv(survival_days) ~ cluster, data = data2)

a<-surv_pvalue(fit)
b<-a$pval

names3<-paste(input,".PNvsMS.pdf") 
print (c(names3,b))
pdf(names3)
ggsurvplot(fit, data = data2,pval = TRUE,break.y.by = 0.2)
dev.off()

data2<-data[c(which(data$cluster =="MS"),which(data$cluster =="CL")),]
fit <- survfit(Surv(survival_days) ~ cluster, data = data2)

a<-surv_pvalue(fit,data2)
b<-a$pval

names4<-paste(input,".MSvsCL.pdf")
print (c(names4,b))
pdf(names4)
ggsurvplot(fit, data = data2,pval = TRUE)
dev.off()
