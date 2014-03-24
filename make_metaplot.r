# Use Andre's bigWig to make a metaplot!
# bedfile: regions we want to get the metaplot around
# bw files: plus and minus signal files (e.g. chipseq, groseq, etc.)
#
# note: here we are doing appropriate normalisation using a spikein (the denominators)
get_rd<-function(datatype,datatime){
    if (datatype=="FP"){
        rds<-list("2min"=35290786/3528924,"5min"=37460495/3740931,"12.5min"=38534246/4854250,"25min"=33348710/5599230,"50min"=29422671/4845825)
    } else if(datatype=="TRP"){
        rds<-list("12.5min"=14299837/696930,"25min"=8798933/805210,"50min"=4463951/558215)
    } else if(datatype=="controls"){
        rds<-list("untreated"=43409559/4341312,"DMSO"=18091023/613884)
    } else if(datatype=="ChIPseq"){
        rds<-list("H3K4me1"=1,"H3K4me3"=1,"H3K27ac"=1)
    }else{
        print("what crazy datatype did you give me?")
        print(datatype)
        quit("no")
    }
    cat(paste("Read depth for",datatype,"and",datatime,"is",rds[datatime],"\n"))
    return(rds[datatime])
}

library(bigWig)
library(ggplot2)
library(grid)
mytheme<-theme(plot.margin=unit(c(2,2,2,2),"cm"),axis.title.x=element_text(vjust=-0.2),axis.title.y=element_text(angle=90,vjust=0.2),plot.title=element_text(size=15,vjust=1.4))+theme_bw()
bw_path<-'/Users/stephanie/ll/data/'

# Parse the input!
args<-commandArgs(TRUE)
bedfilepath<-args[1]

datatype<-args[2]
if(datatype=="ChIPseq"){
    bothstrands<-FALSE
} else{
    bothstrands<-TRUE
}
datatime<-args[3]
read_depth<-as.numeric(get_rd(datatype,datatime))

titleinfo<-args[4]
filename<-paste0(datatype,"_",datatime,"_",titleinfo)

if (file.exists(paste0(filename,".txt"))){
    cat("make_metaplot.r: Metadata already exists! Loading from there. :)\n")
    dl<-read.table(paste0(filename,".txt"),header=T) 
} else{
    cat("make_metaplot.r: No metadata found! Creating...\n")
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
    r_p<-meta.normalize(result_plus,1/read_depth)
    cat("make_metaplot.r: minus strand!\n")
    result_minus<-meta.subsample(bedfile,minus,plus,halfWindow=hw,step=stepp)
    r_m<-meta.normalize(result_minus,1/read_depth)

    # prepare for plotting
    p<-colMeans(r_p[[1]])
    m<-(-1)*colMeans(r_m[[1]])
    col<-rep(c(TRUE,FALSE),each=length(p))
    lv<-c(p,m)
    uq_p<-rep(r_p[[2]],2)
    uq_m<-rep((-1)*r_m[[2]],2)
    lq_p<-rep(r_p[[3]],2)
    lq_m<-rep((-1)*r_m[[3]],2)
    d<-rep(seq(-length(p)/2,length(p)/2-1))
    dl<-data.frame(d,lv,col,uq_p,uq_m,lq_p,lq_m)
    write.table(dl,row.names=F,col.names=T,quote=F,file=paste0(filename,".txt"))
}

cat("make_metaplot.r: data constructed! Ready to plot.\n")
ggplot(dl,aes(x=d,y=lv,colour=col))+geom_point(cex=0.5,show_guide=FALSE)+mytheme+xlab("Distance to centre (bp)")+ylab("Average normalised GRO-seq signal")+ggtitle(paste0("Metagene plot (",datatype,", ",datatime,", ",titleinfo,")"))+scale_colour_manual(values=c("blue","red"))+geom_ribbon(aes(ymin=lq_p,ymax=uq_p),alpha=0.5,fill="red",colour=NA)+geom_ribbon(aes(ymin=lq_m,ymax=uq_m),alpha=0.5,fill="blue",colour=NA)+geom_line(y=0,linetype="dashed",show_guide=FALSE,colour="black")+ylim(-0.01,0.01)
cat(paste0("Saving to ",filename,".pdf"),"\n")
ggsave(paste0(filename,".pdf"))
