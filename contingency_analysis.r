# Test for enrichment of hits in gene bodies for later FP timepoints...
#
time_thresh<-list("X0"=0,"X2"=0,"X5"=0,"X12.5"=7000,"X25"=30000,"X50"=100000)

# Different ways of getting contingency tables
upto_and_after<-function(data,time_points,nt,t,region){
    yes_upto<-sum((data[,region]!=0)&(!rowSums(cbind(rep(0,nrow(data)),data[,time_points[1:t]]))==0))
    no_upto<-sum((data[,region]==0)&(!rowSums(cbind(rep(0,nrow(data)),data[,time_points[1:t]]))==0))
    #yes_upto<-sum((data[,region]!=0)&(rowSums(cbind(rep(0,nrow(data)),data[,time_points[(t+1):nt]]))==0))
    #no_upto<-sum((data[,region]==0)&(rowSums(cbind(rep(0,nrow(data)),data[,time_points[(t+1):nt]]))==0))
    yes_after<-sum((data[,region]!=0)&(rowSums(cbind(rep(0,nrow(data)),data[,time_points[1:t]]))==0))
    no_after<-sum((data[,region]==0)&(rowSums(cbind(rep(0,nrow(data)),data[,time_points[1:t]]))==0))
    contingency.table<-matrix(c(yes_upto,yes_after,no_upto,no_after),2,2,byrow=TRUE)
    return(contingency.table)
}

compare_to_baseline<-function(data,time_points,nt,t,region){
    yes_base<-sum((data[,region]!=0)&(data[,time_points[1]]==1)&(data[,time_points[t+1]]==0))
    no_base<-sum((data[,region]==0)&(data[,time_points[1]]==1)&(data[,time_points[t+1]]==0))
    yes_when<-sum((data[,region]!=0)&(data[,time_points[1]]==0)&(data[,time_points[t+1]]==1))
    no_when<-sum((data[,region]==0)&(data[,time_points[1]]==0)&(data[,time_points[t+1]]==1))
    contingency.table<-matrix(c(yes_base,yes_when,no_base,no_when),2,2,byrow=TRUE)
    return(contingency.table)
}

# this one is really dumb tbh
presence<-function(data,time_points,nt,t,region){
    yes_no<-sum((data[,region]!=0)&(data[,time_points[t]]==0))
    no_no<-sum((data[,region]==0)&(data[,time_points[t]]==0))
    yes_yes<-sum((data[,region]!=0)&(data[,time_points[t]]==1))
    no_yes<-sum((data[,region]==0)&(data[,time_points[t]]==1))
    contingency.table<-matrix(c(yes_no,yes_yes,no_no,no_yes),2,2,byrow=TRUE)
    return(contingency.table)
}

# influenced
influenced_by_FP<-function(data,time_points,nt,t){
    # exclude the hits which have a distance of 0, because what even is that
    d<-data[((data$"gene_body"!=0)&(data[,time_points[t+1]]==1)&(data$"distance"!=0)),]
    #d<-data[((data$"gene_body"!=0)&(data[,time_points[t+1]]==1)),]

#    wt_inf<-sum((d[,time_points[1]]==1)&(d$"distance"<=time_thresh[time_points[t+1]]))
#    wt_ninf<-sum((d[,time_points[1]]==1)&(d$"distance">time_thresh[time_points[t+1]]))
#    nwt_inf<-sum((d[,time_points[1]]==0)&(d$"distance"<=time_thresh[time_points[t+1]]))
#    nwt_ninf<-sum((d[,time_points[1]]==0)&(d$"distance">time_thresh[time_points[t+1]]))
    
    wt_inf<-sum((rowSums(cbind(rep(0,nrow(d)),d[,time_points[1:t]]))!=0)&(d$"distance"<=time_thresh[time_points[t+1]]))
    wt_ninf<-sum((rowSums(cbind(rep(0,nrow(d)),d[,time_points[1:t]]))!=0)&(d$"distance">time_thresh[time_points[t+1]]))
    nwt_inf<-sum((rowSums(cbind(rep(0,nrow(d)),d[,time_points[1:t]]))==0)&(d$"distance"<=time_thresh[time_points[t+1]]))
    nwt_ninf<-sum((rowSums(cbind(rep(0,nrow(d)),d[,time_points[1:t]]))==0)&(d$"distance">time_thresh[time_points[t+1]]))

    contingency.table<-matrix(c(wt_inf,nwt_inf,wt_ninf,nwt_ninf),2,2,byrow=TRUE)
    return(contingency.table)
}
datapath<-'/Users/stephanie/ll/results/40m_SVM/dREG_regions_confident_0.8.bed.gz'
#datapath<-'/Users/stephanie/ll/results/40m_SVM/dREG_regions_maybe_0.5.bed.gz'
#datapath<-'/Users/stephanie/ll/results/40m_SVM/dREG_regions_preQC_0.8.bed.gz'

data<-read.table(datapath,header=T)

time_points<-c("X0","X2","X5","X12.5","X25","X50")
nt<-length(time_points)
gs.pvals<-list()
gb.pvals<-list()
ng.pvals<-list()
fp.pvals<-list()
gs.tables<-list()
gb.tables<-list()
ng.tables<-list()
fp.tables<-list()
for (t in seq(nt-1)){
    # test gene_start enrichment
    gs.contingency.table<-compare_to_baseline(data,time_points,nt,t,"gene_start")
    chisq.result<-chisq.test(gs.contingency.table)
    stat<-chisq.result$statistic
    pval<-chisq.result$p.value
    cat("Gene start?\n")
    print(gs.contingency.table)
    gs.pvals[time_points[t+1]]<-pval
    gs.tables[[time_points[t+1]]]<-gs.contingency.table

    # test gene_body enrichment
    gb.contingency.table<-compare_to_baseline(data,time_points,nt,t,"gene_body")
    chisq.result<-chisq.test(gb.contingency.table)
    stat<-chisq.result$statistic
    pval<-chisq.result$p.value
    cat("gene body?\n")
    print(gb.contingency.table)
    cat("Chi-sq stat:",stat,"\n")
    cat("p-value of association:",pval,"\n")
    gb.pvals[time_points[t+1]]<-pval
    gb.tables[[time_points[t+1]]]<-gb.contingency.table

    # test non_gene enrichment
    ng.contingency.table<-compare_to_baseline(data,time_points,nt,t,"non_gene")
    chisq.result<-chisq.test(ng.contingency.table)
    stat<-chisq.result$statistic
    pval<-chisq.result$p.value
    cat("non_gene enrichment?\n")
    print(ng.contingency.table)
    cat("Chi-sq stat:",stat,"\n")
    cat("p-value of association:",pval,"\n")
    ng.pvals[time_points[t+1]]<-pval
    ng.tables[[time_points[t+1]]]<-ng.contingency.table

    fp.contingency.table<-influenced_by_FP(subset(data,gene_body!=0),time_points,nt,t)
    chisq.result<-chisq.test(fp.contingency.table)
    stat<-chisq.result$statistic
    pval<-chisq.result$p.value
    cat("effects of FP?\n")
    print(fp.contingency.table)
    cat("Chi-sq stat:",stat,"\n")
    cat("p-value of association:",pval,"\n")
    fp.pvals[time_points[t+1]]<-pval
    fp.tables[[time_points[t+1]]]<-fp.contingency.table
}
