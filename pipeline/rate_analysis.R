# This takes the results of run_HMM.sh, extracts rates, and (will) perform analysis.
library(ggplot2)

# --- load functions --- #
source('fns_rate_analysis.r')

# --- options --- #
#args<-commandArgs(TRUE)
suffix<-"40m_SVM"
#time_point<-args[1]
#time_point<-"early"
time_point<-"mid"
#time_point<-"late"
replicates<-c("B","R")
#replicates<-c("B","B")
#replicates<-c("R","R")
#replicates<-c("C","C")
get_covariates<-FALSE
#get_covariates<-TRUE

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
get_region_of_relevance(five,twelve,genes,"early")
get_region_of_relevance(twelve,twentyfive,genes,"mid")
get_region_of_relevance(twentyfive,fifty,genes,"late")

# --- plot the transitions--- #
plot_transitions(five,twelve,twentyfive,fifty)

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
} else{
    print("what timepoint did you give me? :(")
    quit("no")
}
cat(nrow(data),"genes at",time_point,"time point.\n")

# --- combine the covariates --- #
data_covar<-add_covars(data,n_gb,int1len,n_exon,CpG_gb,H3K79me2_gb)
cat(nrow(data_covar),"genes with covar data.\n")
cat("Correlation between rate and n_exon:",cor(data_covar$rate,data_covar$n_exon),"\n")
cat("Correlation between rate and n_gb:",cor(data_covar$rate,data_covar$n_gb),"\n")

# --- only dREG hits --- #
data_withhit<-subset(data_covar,n_gb>0)
cat("There are",nrow(data_withhit),"genes with dREG hits.\n")
cat("Correlation between rate and n_exon:",cor(data_withhit$rate,data_withhit$n_exon),"\n")
cat("Correlation between rate and n_gb:",cor(data_withhit$rate,data_withhit$n_gb),"\n")

# --- do comparisons --- #
in_theirs<-merge(data_covar,comp_rates,by=1)
cat(nrow(in_theirs),"genes overlap with previous results.\n")
cat("Correlation between rates:",cor(in_theirs$rate.x,in_theirs$rate.y),"\n")
cat("Correlation between rate and n_exon (overlapped):",cor(in_theirs$rate.x,in_theirs$n_exon),"\n")

new<-!data_covar$name%in%comp_rates$geneID
data_new<-subset(data_covar,new)

# --- make model --- #
m_full<-lm(rate ~ n_gb + int1len + n_exon + CpG + H3K79me2,data=data_covar)
m_withhit<-lm(rate ~ n_gb + int1len + n_exon + CpG + H3K79me2,data=data_withhit)
