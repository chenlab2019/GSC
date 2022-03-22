args<-commandArgs(T)
input<-args[1]
output<-args[2]

table<-read.table(input,header=F,sep="\t")
table$V3[grep("1",table$V3)]<-"mp4"
table$V3[grep("2",table$V3)]<-"mp2"
table$V3[grep("3",table$V3)]<-"mp1"
table$V3[grep("mp4",table$V3)]<-"mp3"
table$V2[grep("1",table$V2)]<-"atac1"
table$V2[grep("2",table$V2)]<-"atac2"
table$V2[grep("3",table$V2)]<-"atac3"

a<-table(table[,2],table[,3])
write.table(a,file=output,quote=F,sep="\t")
