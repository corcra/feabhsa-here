args<-commandArgs(TRUE)

genes<-read.table(args[1])
from<-read.table(args[2],header=T)
to<-read.table(args[3],header=T)

# merge the timepoints... need info from both
both<-merge(from,to,by=2)

# add the strand information
strand_info<-genes[,c(4,6)]
names(strand_info)<-c("name","strand")
both<-merge(both,strand_info,by=1)

# look at the transitions...
starts<-ifelse(both$strand=="+",floor(both$gene_start.x+both$transition.x),floor(both$gene_end.x-both$transition.y))
ends<-ifelse(both$strand=="+",ceiling(both$gene_start.x+both$transition.y),ceiling(both$gene_end.x-both$transition.x))
chr<-both$chr.x
name<-both$name
filler<-rep("NA",length(name))

# combinate
region<-data.frame(chr,starts,ends,name,filler)
region<-subset(region,ends>starts)

# save
write.table(region,file="regions.bed.temp",row.names=F,col.names=F,quote=F)
