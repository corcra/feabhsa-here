# This takes the results of run_HMM.sh, extracts rates, and performs analysis.
library(ggplot2)
library(ggm)
library(grid)
library(scales)

# --- load functions --- #
source('fns_rate_analysis.r')
mytheme<-theme(plot.margin=unit(c(2,2,2,2),"cm"),axis.title.x=element_text(vjust=-0.2),axis.title.y=element_text(angle=90,vjust=0.2),plot.title=element_text(size=15,vjust=1.4))+theme_bw()

# --- options --- #
#args<-commandArgs(TRUE)
suffix<-"40m_SVM"
#time_point<-args[1]
#time_point<-"early"
time_point<-"mid"
#time_point<-"late"
#time_point<-"to12"
#time_point<-"to25"
replicates<-c("B","R")
#replicates<-c("B","B")
#replicates<-c("R","R")
#replicates<-c("C","C")
#get_covariates<-TRUE
get_covariates<-FALSE
#get_region<-TRUE
get_region<-FALSE

# --- gene data --- #
genes<-read.table(paste0("/Users/stephanie/ll/results/",suffix,"/active_genes_Rable.bed"),as.is=TRUE)

# ---- rate data (replicates) ---- #
five_1<-read.table(paste0('/Users/stephanie/ll/from_hojoong/result/',replicates[1],'_5min_bin500.',suffix,'.txt'),header=T,as.is=TRUE)
five_2<-read.table(paste0('/Users/stephanie/ll/from_hojoong/result/',replicates[2],'_5min_bin500.',suffix,'.txt'),header=T,as.is=TRUE)
twelve_1<-read.table(paste0('/Users/stephanie/ll/from_hojoong/result/',replicates[1],'_1230min_bin1000.',suffix,'.txt'),header=T,as.is=TRUE)
twelve_2<-read.table(paste0('/Users/stephanie/ll/from_hojoong/result/',replicates[2],'_1230min_bin1000.',suffix,'.txt'),header=T,as.is=TRUE)
twentyfive_1<-read.table(paste0('/Users/stephanie/ll/from_hojoong/result/',replicates[1],'_25min_bin2000.',suffix,'.txt'),header=T,as.is=TRUE)
twentyfive_2<-read.table(paste0('/Users/stephanie/ll/from_hojoong/result/',replicates[2],'_25min_bin2000.',suffix,'.txt'),header=T,as.is=TRUE)
fifty_1<-read.table(paste0('/Users/stephanie/ll/from_hojoong/result/',replicates[1],'_50min_bin5000.',suffix,'.txt'),header=T,as.is=TRUE)
fifty_2<-read.table(paste0('/Users/stephanie/ll/from_hojoong/result/',replicates[2],'_50min_bin5000.',suffix,'.txt'),header=T,as.is=TRUE)

# --- restrict to genes consistent across replicates --- #
five<-consistent_over_replicates(five_1,five_2,"five")
twelve<-consistent_over_replicates(twelve_1,twelve_2,"twelve")
twentyfive<-consistent_over_replicates(twentyfive_1,twentyfive_2,"twentyfive")
fifty<-consistent_over_replicates(fifty_1,fifty_2,"fifty")

# --- save the region_of_relevance (formerly get_FP_affected_region.r) --- #
if (get_region){
    get_region_of_relevance(five,twelve,genes,"early",1000,30000)
    get_region_of_relevance(twelve,twentyfive,genes,"mid",7500,30000)
    get_region_of_relevance(twentyfive,fifty,genes,"late",NA,NA)
    get_region_of_relevance(NA,twelve,genes,"to12",500,7500)
    get_region_of_relevance(NA,twentyfive,genes,"to25",500,30000)
}

# --- plot the transitions--- #
#plot_transitions(five,twelve,twentyfive,fifty)

# --- comparison data --- #
comp_folder<-'/Users/stephanie/ll/from_hojoong/regression/data/'
comp_rates<-read.table(paste0(comp_folder,'level-rate.txt'),header=T)
comp_rates<-subset(comp_rates,rate>0)

# --- other covariates --- #
if(get_covariates){
    cat("Getting covariates!\n")
    system(paste0("./get_covar_minimal.sh ",time_point))
}

covar_folder<-paste0('/Users/stephanie/ll/results/rate_analysis/',time_point,'/')
n_gb<-read.table(paste0(covar_folder,'n_gb.txt'),header=T)
int1len<-read.table(paste0(covar_folder,'int1len.txt'),header=T)
n_exon<-read.table(paste0(covar_folder,'n_exon.txt'),header=T)
CpG_gb<-read.table(paste0(covar_folder,'CpG_gb.txt'),header=T)
H3K79me2_gb<-read.table(paste0(covar_folder,'H3K79me2_gb.txt'),header=T)

# --- get rates --- #
if (time_point=="early"){
    data<-get_rate(five,twelve,7.5)
} else if(time_point=="mid"){
    data<-get_rate(twelve,twentyfive,12.5)
} else if(time_point=="late"){
    data<-get_rate(twentyfive,fifty,25)
} else if(time_point=="to12"){
    data<-get_rate(NA,twelve,12.5)
} else if(time_point=="to25"){
    data<-get_rate(NA,twentyfive,25)
} else{
    print("what timepoint did you give me? :(")
    quit("no")
}
cat(nrow(data),"genes at",time_point,"time point.\n")

# --- delete that one weird point --- #
data<-subset(data,rate<4)

# --- combine the covariates --- #
data_covar<-add_covars(data,n_gb,int1len,n_exon,CpG_gb,H3K79me2_gb)
cat(nrow(data_covar),"genes with covar data.\n")
cat("Correlation between rate and n_exon:",cor(data_covar$rate,data_covar$n_exon,method="spearman"),"\n")
cat("Correlation between rate and n_gb:",cor(data_covar$rate,data_covar$n_gb,method="spearman"),"\n")

# --- only dREG hits --- #
data_withhit<-subset(data_covar,n_gb>0)
cat("There are",nrow(data_withhit),"genes with dREG hits.\n")
#cat("Correlation between rate and n_exon:",cor(data_withhit$rate,data_withhit$n_exon,method="spearman"),"\n")
#cat("Correlation between rate and n_gb:",cor(data_withhit$rate,data_withhit$n_gb,method="spearman"),"\n")

# --- only with exons --- #
data_withexon<-subset(data_covar,n_exon>0)
#cat("There are",nrow(data_withexon),"genes with exons.\n")
#cat("Correlation between rate and n_exon:",cor(data_withexon$rate,data_withexon$n_exon,method="spearman"),"\n")
#cat("Correlation between rate and n_gb:",cor(data_withexon$rate,data_withexon$n_gb,method="spearman"),"\n")

# --- do comparisons --- #
in_theirs<-merge(data_covar,comp_rates,by=1)
cat(nrow(in_theirs),"genes overlap with previous results.\n")
cat("Correlation between rates:",cor(in_theirs$rate.x,in_theirs$rate.y,method="spearman"),"\n")
cat("Correlation between rate and n_exon (overlapped):",cor(in_theirs$rate.x,in_theirs$n_exon,method="spearman"),"\n")

new<-!data_covar$name%in%comp_rates$geneID
data_new<-subset(data_covar,new)

# --- make model --- #
m_full<-lm(rate ~ n_gb + int1len + n_exon + CpG + H3K79me2,data=data_covar)
m_withhit<-lm(rate ~ n_gb + int1len + n_exon + CpG + H3K79me2,data=data_withhit)
m_withexon<-lm(rate ~ n_gb + int1len + n_exon + CpG + H3K79me2,data=data_withexon)
roz.m_full<-lm(roz(rate) ~ roz(n_gb) + roz(int1len) + roz(n_exon) + roz(CpG) + roz(H3K79me2),data=data_covar)
roz.m_withhit<-lm(roz(rate) ~ roz(n_gb) + roz(int1len) + roz(n_exon) + roz(CpG) + roz(H3K79me2),data=data_withhit)

# --- analyse model --- #
qplot(sample=m_full$residuals/sqrt(var(m_full$residuals)))+geom_abline(intercept=0,slope=1,linetype="dashed",col="red")+mytheme+xlab("Normal quantiles")+ylab("Observed quantiles")+ggtitle(paste0("QQ Plot of Residuals (",time_point,")"))
ggsave(paste0("residual_qqplot_",time_point,".pdf"))

# --- plot some covars --- #
ggplot(data_covar,aes(x=n_exon,y=rate,group=n_exon))+geom_boxplot()+mytheme+xlab("Number of exons")+ylab("Rate (kb/min)")+ggtitle(paste0("Rate as a function of exon number (",time_point,")"))
ggsave(paste0("boxplot_exon_",time_point,".pdf"))
ggplot(data_covar,aes(x=n_gb,y=rate,group=n_gb))+geom_boxplot()+mytheme+xlab("Number of dREG hits")+ylab("Rate (kb/min)")+ggtitle(paste0("Rate as a function of dREG hit number (",time_point,")"))
ggsave(paste0("boxplot_gb_",time_point,".pdf"))

# --- compare covariates --- #
#ggplot(data_covar,aes(x=int1len,fill=n_gb>0))+geom_histogram(position="identity",aes(y=..density..),alpha=0.7,binwidth=5000)+theme_bw()+ggtitle("Length of intron 1 (with/without dREG hit)")
#ggsave("int1len_hist.pdf")
#ggplot(data_covar,aes(x=n_exon,fill=n_gb>0))+geom_histogram(position="identity",aes(y=..density..),alpha=0.7,binwidth=0.02)+theme_bw()+ggtitle("Exon density (per kb) (with/without dREG hit)")
#ggsave("n_exon_hist.pdf")
#ggplot(data_covar,aes(x=CpG,fill=n_gb>0))+geom_histogram(position="identity",aes(y=..density..),alpha=0.7,binwidth=1)+theme_bw()+ggtitle("CpG density (per kb) (with/without dREG hit)")
#ggsave("CpG_hist.pdf")
#ggplot(data_covar,aes(x=H3K79me2,fill=n_gb>0))+geom_histogram(position="identity",aes(y=..density..),alpha=0.7,binwidth=2)+theme_bw()+ggtitle("H3K79me2 density (per kb) (with/without dREG hit)")
#ggsave("H3K79me2_hist.pdf")
#ggplot(data_covar,aes(x=rate,fill=n_gb>0))+geom_histogram(position="identity",aes(y=..density..),alpha=0.7,binwidth=0.1)+theme_bw()+ggtitle("rate (kb/min) (with/without dREG hit)")
#ggsave("rate_hist.pdf")

# --- gaussian graphical model --- #
roz.covar_array<-get_covar_array(data_covar)
ggm<-get_ggm(data_covar$rate,roz.covar_array,0.05)
g<-ggm$"g"
plot(g,vertex.size=8,vertex.color="gray90",vertex.label.family="Helvetica",vertex.label.cex=0.5,vertex.label.color="green3",vertex.frame.color=NA)

# --- nonlinear model (from Andre) --- #
s=list("b_gb"=-7,"b_exon"=-3,"baserate"=2.5,"b_int1len"=0.01,"b_CpG"=-1,"b_H3K79me2"=0.3)
nl<-nls(distance~(12.5+(b_gb*n_gb+b_exon*n_exon))*(baserate+b_int1len*int1len+b_CpG*CpG+b_H3K79me2*H3K79me2),data=data_covar,start=s)
# (fuller)
sp=list("b_n_gb"=-7,"b_d_gb"=-0.3,"b_n_exon"=-3,"b_d_exon"=-0.1,"baserate"=2.5,"b_int1len"=0.01,"b_CpG"=-1,"b_H3K79me2"=0.3)
nlp<-nls(distance~(12.5+(b_n_gb*n_gb+b_n_exon*n_exon))*(baserate+b_int1len*int1len+b_CpG*CpG+b_H3K79me2*H3K79me2+b_d_gb*d_gb+b_d_exon*d_exon),data=data_covar,start=sp)

# (using only exons and no dREG)
e=list("b_n_exon"=-3,"baserate"=2.5,"b_int1len"=0.01,"b_CpG"=-1,"b_H3K79me2"=0.3)
ne<-nls(distance~(12.5+(b_n_exon*n_exon))*(baserate+b_int1len*int1len+b_CpG*CpG+b_H3K79me2*H3K79me2),data=data_covar,start=e)
ep=list("b_n_exon"=-3,"b_d_exon"=-0.3,"baserate"=2.5,"b_int1len"=0.01,"b_CpG"=-1,"b_H3K79me2"=0.3)
nep<-nls(distance~(12.5+(b_n_exon*n_exon))*(baserate+b_int1len*int1len+b_CpG*CpG+b_H3K79me2*H3K79me2+b_d_exon*d_exon),data=data_covar,start=ep)
