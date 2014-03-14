# Test for enrichment of hits in gene bodies for later FP timepoints...
#

# Different ways of getting contingency tables
upto_and_after<-function(data,time_points,nt,t,region){
    yes_upto<-sum((data[,region]!=0)&(rowSums(cbind(rep(0,nrow(data)),data[,time_points[(t+1):nt]]))==0))
    no_upto<-sum((data[,region]==0)&(rowSums(cbind(rep(0,nrow(data)),data[,time_points[(t+1):nt]]))==0))
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
    print(contingency.table)
    return(contingency.table)
}

datapath<-'/Users/stephanie/ll/results/FP/dREG_regions_marked.bed.gz'
data<-read.table(datapath,header=T)

time_points<-c("X0","X2","X5","X12.5","X25","X50")
nt<-length(time_points)
gs.pvals<-list()
gb.pvals<-list()
ng.pvals<-list()
gs.tables<-list()
gb.tables<-list()
ng.tables<-list()
for (t in seq(nt-1)){
    # test gene_start enrichment
    gs.contingency.table<-compare_to_baseline(data,time_points,nt,t,"gene_start")
    chisq.result<-chisq.test(gs.contingency.table)
    stat<-chisq.result$statistic
    pval<-chisq.result$p.value
    cat("Gene start?\n")
    print(gs.contingency.table)
    print(stat)
    print(pval)
    gs.pvals[time_points[t]]<-pval
    gs.tables[[time_points[t]]]<-gs.contingency.table

    # test gene_body enrichment
    gb.contingency.table<-compare_to_baseline(data,time_points,nt,t,"gene_body")
    chisq.result<-chisq.test(gb.contingency.table)
    stat<-chisq.result$statistic
    pval<-chisq.result$p.value
    cat("gene body?\n")
    print(gb.contingency.table)
    print(stat)
    print(pval)
    gb.pvals[time_points[t]]<-pval
    gb.tables[[time_points[t]]]<-gb.contingency.table

    # test non_gene enrichment
    ng.contingency.table<-compare_to_baseline(data,time_points,nt,t,"non_gene")
    chisq.result<-chisq.test(ng.contingency.table)
    stat<-chisq.result$statistic
    pval<-chisq.result$p.value
    cat("non_gene enrichment?\n")
    print(ng.contingency.table)
    print(stat)
    print(pval)
    ng.pvals[time_points[t]]<-pval
    ng.tables[[time_points[t]]]<-ng.contingency.table

}