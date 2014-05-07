# This takes the results of run_HMM.sh, extracts rates, and (will) perform analysis.
library(ggplot2)

#args<-commandArgs(TRUE)
suffix<-"40m_SVM"
#time_point<-args[1]
time_point<-"mid"

# ---- rate data ---- #
five<-read.table(paste0('/Users/stephanie/ll/from_hojoong/result/V6.5_5min_bin500.',suffix,'.txt'),header=T,as.is=TRUE)
twelve<-read.table(paste0('/Users/stephanie/ll/from_hojoong/result/V6.5_12.5min_bin1000.',suffix,'.txt'),header=T,as.is=TRUE)
twentyfive<-read.table(paste0('/Users/stephanie/ll/from_hojoong/result/V6.5_25min_bin2000.',suffix,'.txt'),header=T,as.is=TRUE)
fifty<-read.table(paste0('/Users/stephanie/ll/from_hojoong/result/V6.5_50min_bin5000.',suffix,'.txt'),header=T,as.is=TRUE)
time<-c(rep("5min",nrow(five)),rep("12.5min",nrow(twelve)),rep("25min",nrow(twentyfive)),rep("50min",nrow(fifty)))
time<-factor(time,c("5min","12.5min","25min","50min"))
all<-rbind(five,twelve,twentyfive,fifty)
all<-data.frame(all,time)
ggplot(all,aes(x=transition/1000,fill=time))+geom_histogram(binwidth=2)+facet_grid(time~.,scales="free")+xlab("transition point (kb)")+theme_bw()+scale_fill_manual(values=rainbow(4,v=0.8,s=0.5))
ggsave("transitions.pdf")

# --- other covariates --- #
covar_folder<-paste0('/Users/stephanie/ll/results/rate_analysis/',time_point,'/')
n_gb<-read.table(paste0(covar_folder,'n_gb.txt'),header=T)
int1len<-read.table(paste0(covar_folder,'int1len.txt'),header=T)
n_exon<-read.table(paste0(covar_folder,'n_exon.txt'),header=T)
#
# --- functions --- #
get_rate<-function(from,to,timediff){
    # note, the name field should be the second one...
    both<-merge(from,to,by=2)
    name<-both$name
    chr<-both$chr.x
    start<-both$gene_start.x
    end<-both$gene_end.x
    rate<-(both$transition.y-both$transition.x)/(timediff*1000)
    name<-both$name
    df<-data.frame(chr,start,end,name,rate)
    df_pos<-subset(df,rate>0)
    return(df_pos)
}

# taken from hojoong's script
roz<-function(d)
{
    r<-qnorm(rank(d+runif(length(d),0,0.1))/(length(d)+1))
    return(r)
}

add_covars<-function(data,n_gb,int1len,n_exon){
    print("Adding covariates!")
    data<-merge(data,n_gb,by.x=4,by.y=1)
    data<-merge(data,int1len,by=1)
    data<-merge(data,n_exon,by=1)

    #data$rate<-roz(data$rate)
    #data$n_gb<-roz(data$n_gb)
    #data$int1len<-roz(data$int1len)
    #data$n_exon<-roz(data$n_exon)
    return(data)
}

rateplot<-function(data,title){
    ggplot(data,aes(x=n_gb,y=rate,group=n_gb))+geom_boxplot()+theme_bw()+xlab("# dREG hits")+ylab("rate (kb/min)")+ggtitle(title)
}
# --- get rates --- #
if (time_point=="early"){
    data<-get_rate(five,twelve,7.5)
} else if(time_point=="mid"){
    data<-get_rate(twelve,twentyfive,12.5)
} else if(time_point=="late"){
    data<-get_rate(twentyfive,fifty,25)
} else{
    print("what timepoint did you give me? :(")
    quit("no")
}

# --- combine the covariates --- #
data_covar<-add_covars(data,n_gb,int1len,n_exon)

# --- do comparisons --- #
