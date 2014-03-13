# Use Andre's bigWig to make a metaplot!
# bedfile: regions we want to get the metaplot around
# bw files: plus and minus signal files (e.g. chipseq, groseq, etc.)

library(bigWig)
library(ggplot2)
library(grid)
mytheme<-theme(plot.margin=unit(c(2,2,2,2),"cm"),axis.title.x=element_text(vjust=-0.2),axis.title.y=element_text(angle=90,vjust=0.2),plot.title=element_text(size=15,vjust=1.4))+theme_bw()
bw_path<-'/Users/stephanie/ll/data/'

# Parse the input!
args<-commandArgs(TRUE)
bedfilepath<-args[1]

datatype<-args[2]
datatime<-args[3]
filename<-args[4]

bw_plus<-paste0(bw_path,datatype,'/V6.5_',datatime,ifelse(datatype=="controls","",datatype),'_Plus.bw')
bw_minus<-paste0(bw_path,datatype,'/V6.5_',datatime,ifelse(datatype=="controls","",datatype),'_Minus.bw')

# Constants

hw<-2000
stepp<-1

# Load the data!
plus<-load.bigWig(bw_plus)
minus<-load.bigWig(bw_minus)
bedfile<-read.table(bedfilepath)

# in the event that bedfile has no strand information, we should be ignoring the second bw argument in either of these... 
cat("make_metaplot.r: plus strand!\n")
result_plus<-meta.subsample(bedfile,plus,minus,halfWindow=hw,step=stepp)
cat("make_metaplot.r: minus strand!\n")
result_minus<-meta.subsample(bedfile,minus,plus,halfWindow=hw,step=stepp)

# prepare for plotting
p<-colMeans(result_plus[[1]])
m<-(-1)*colMeans(result_minus[[1]])
col<-rep(c(TRUE,FALSE),each=length(p))
lv<-c(p,m)
uq_p<-rep(result_plus[[2]],2)
uq_m<-rep((-1)*result_minus[[2]],2)
lq_p<-rep(result_plus[[3]],2)
lq_m<-rep((-1)*result_minus[[3]],2)
d<-rep(seq(-length(p)/2,length(p)/2-1))
dl<-data.frame(d,lv,col,uq_p,uq_m,lq_p,lq_m)
write.table(dl,row.names=F,col.names=T,quote=F,file=paste0(filename,".txt"))

cat("make_metaplot.r: data constructed! Ready to plot.\n")
ggplot(dl,aes(x=d,y=lv,colour=col))+geom_point(cex=0.5,show_guide=FALSE)+mytheme+xlab("Distance to centre (bp)")+ylab("Average GRO-seq signal")+ggtitle(paste0("Metagene plot (",datatype,", ",datatime,")"))+scale_colour_manual(values=c("blue","red"))+geom_ribbon(aes(ymin=lq_p,ymax=uq_p),alpha=0.5,fill="red",colour=NA)+geom_ribbon(aes(ymin=lq_m,ymax=uq_m),alpha=0.5,fill="blue",colour=NA)+geom_line(y=0,linetype="dashed",show_guide=FALSE,colour="black")
ggsave(paste0(filename,".pdf"))
