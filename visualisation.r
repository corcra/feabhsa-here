# Makes histograms: region sizes, distance to TSS, dREG scores
library(ggplot2)
library(grid)
source('vis_fns.r')

#args<-commandArgs(TRUE)
#datapath<-args[1]
datapath<-'/Users/stephanie/ll/results/FP/dREG_regions_marked.bed.gz'
data<-read.table(datapath,header=T)
n<-nrow(data)

# Timepoints (optionally from argument)
time_points<-c("X0","X2","X5","X12.5","X25","X50")
time_readable<-c("WT","2min","5min","12.5min","25min","50min")
time_thresh<-list("X0"=0,"X2"=0,"X5"=0,"X12.5"=7000,"X25"=30000,"X50"=100000)

# ChIP-seq mods (eg on histones, but could be TFs)
mods<-c("H3K4me1","H3K4me3","H3K27ac")
norms<-list("H3K4me1"=4935147/1e06,"H3K4me3"=8046120/1e06,"H3K27ac"=21839696/1e06)

# General style thing
mytheme<-theme(plot.margin=unit(c(2,2,2,2),"cm"),axis.title.x=element_text(vjust=-0.2),axis.title.y=element_text(angle=90,vjust=0.2),plot.title=element_text(size=15,vjust=1.4))+theme_bw()

#Quick check: how many multiply-assigned regions do we have?
double_regions<-sum((data$TSS+1*(data$gene_body!=0)+data$non_gene)==2)
triple_regions<-sum((data$TSS+1*(data$gene_body!=0)+data$non_gene)==3)
cat("visualisation.r: There are",double_regions,"doubly-assigned and",triple_regions,"triply-assigned regions.\n")

# Condense regions
regions<-ifelse(data$TSS==1,"TSS",ifelse(data$gene_body!=0,"gene_body","non_gene"))
regions<-factor(regions,c("TSS","gene_body","non_gene"))
data<-cbind(data,regions)
region_cols<-c("firebrick1","orange1","turquoise3","gray")

# First appearances
first_appearance<-time_points[apply(data[,time_points],1,function(x) which(x==1)[1])]
first_appearance<-factor(first_appearance,time_points)
data<-cbind(data,first_appearance)

# How often are the regions seen? (aross timepoints)
cat("visualisation.r: How often is each region seen?\n")
howmany_prefac<-rowSums(data[,time_points])
howmany<-factor(howmany_prefac,c(1,2,3,4,5,6))
data_howfreq<-data.frame(regions,howmany,first_appearance)
ggplot(data_howfreq,aes(x=howmany,fill=regions))+geom_histogram(position="dodge")+scale_fill_manual(values=region_cols)+mytheme+ggtitle("How many timepoints do hits occur in?")+xlab("Number of timepoints")+facet_grid(~regions,margins=TRUE)
ggsave("repetition_regions.pdf",width=10)

# For regions that are only seen once... when are they seen?
cat("visualisation.r: Loner regions\n")
data_loners<-data_howfreq[howmany_prefac==1,]
ggplot(data_loners,aes(x=first_appearance,fill=regions))+geom_histogram(position="dodge")+scale_fill_manual(values=region_cols)+mytheme+xlab("Timepoint")+ggtitle("For regions appearing only once, when was it?")
ggsave("loner_regions.pdf",width=10)

# How many violations? (e.g. disappearance after appearance)
#cat("visualisation.r: violations\n")
#violations<-factor(viol_prefac,c(0,1,2,3,4))
#data_viol<-data.frame(regions,violations)
#ggplot(data_viol,aes(x=violations,fill=regions))+geom_histogram(position="dodge")+scale_fill_manual(values=region_cols)+mytheme+facet_grid(~regions,margins=TRUE)
#ggsave("violations_regions.pdf",width=10)

#How consistent are they?
cat("visualisation.r: When do consistent regions appear?\n")
viol_prefac<-get_violations(data[,time_points])
consistent<-viol_prefac==0
consistent<-factor(consistent,c(TRUE,FALSE))
data_consistent<-data.frame(regions,consistent,first_appearance)
ggplot(data_consistent,aes(x=first_appearance,fill=consistent))+facet_grid(consistent~regions)+geom_histogram(position="dodge")+mytheme+scale_fill_manual(values=c("darkorchid1","grey60"))+ggtitle("When do consistent regions appear? (consistent meaning that once they appear, they do not disappear)")
ggsave("consistent_regions.pdf",width=10)

# Now, are they earlier (<12.5) or later(>12.5)? (do this better!) (STRICTLY early v. STRICTLY late...)
cat("visualisation.r: regions only appearing in early/late timepoints?\n")
early<-rowMeans(data[,c("X12.5","X25","X50")])==0
late<-rowMeans(data[,c("X0","X2","X5")])==0
all<-rowMeans(data[,time_points])==1
when_restricted<-ifelse(all,"all",ifelse(early,"early only",ifelse(late,"late only","mixed")))
when_restricted<-factor(when_restricted,c("all","early only","late only","mixed"))
data_whenrestricted<-data.frame(regions,when_restricted)
ggplot(data_whenrestricted,aes(x=when_restricted,fill=regions))+geom_histogram(position="dodge")+mytheme+scale_fill_manual(values=c(region_cols))+xlab("Only appears when...")+ggtitle("Many regions are only seen at early timepoints (WT,2min,5min)")
ggsave("earlylate_regions.pdf",width=10)

# Expand regions - this is a giant dataframe with duplicated entries whenever regions appear in multiple timepoints - the purpose is to be able to simultaneously plot not just new points but all points at a time, acces times via 'when'
data_exp<-data.frame()
for (time in time_points){
    time_data<-data[data[,time]==1,-which(names(data)%in%time_points)]
    en<-nrow(time_data)
    when<-rep(time,en)
    t_d<-data.frame(time_data,when)
    data_exp<-rbind(data_exp,t_d)
}

# When do the regions first appear? (note, they may disappear and reappear - this takes very first appearance)
cat("visualisation.r: First appearance of regions!\n")
ggplot(data,aes(x=first_appearance,fill=regions))+geom_bar(position="dodge")+scale_fill_manual(values=region_cols)+mytheme+facet_grid(~regions)
ggsave("first_appearance.pdf",width=10)
# What fraction of the regions appear for the first time at this timepoint?
new_bool<-data_exp$when==data_exp$first_appearance
new<-ifelse(new_bool,ifelse(data_exp$regions=="TSS", "new (TSS)", ifelse(data_exp$regions=="gene_body", "new (gene body)", "new (non gene)")),"old")
new<-factor(new,c("new (TSS)","new (gene body)","new (non gene)","old"))
new_cols<-c(region_cols[1:3],"grey90")
data_new<-cbind(data_exp,new)
ggplot(data_new,aes(x=when,fill=new))+geom_bar(position="fill")+scale_fill_manual(values=new_cols)+mytheme+facet_grid(~regions)+ylab("Fraction of regions which are new at this timepoint")+xlab("Timepoint")+ggtitle("Does time after flavopiridol treatment aid discovery of intragenic regions?")
ggsave("fraction_new.pdf",width=10)

new<-ifelse(new_bool,ifelse(data_exp$regions=="TSS", "new (TSS)", ifelse(data_exp$regions=="gene_body", "new (gene body)", "new (non gene)")),ifelse(data_exp$regions=="TSS", "old (TSS)", ifelse(data_exp$regions=="gene_body", "old (gene body)","old (non gene)")))
new<-factor(new,c("new (TSS)","old (TSS)","new (gene body)","old (gene body)", "new (non gene)", "old (non gene)"))
data_new<-cbind(data_exp,new)
new_cols<-c(region_cols[1],"firebrick4",region_cols[2],"orange3",region_cols[3],"turquoise4")
ggplot(data_new,aes(x=when,fill=new))+geom_bar(position="dodge")+scale_fill_manual(values=new_cols)+mytheme+facet_grid(~regions)+ylab("Number of regions  at this timepoint")+xlab("Timepoint")
ggsave("total_regions.pdf",width=10)

# Region sizes
cat("visualisation.r: Region sizes!\n")
sizes<-data$end-data$start
size_dat<-data.frame(sizes,regions)
ggplot(size_dat,aes(x=sizes,fill=regions))+geom_histogram(binwidth=50)+mytheme+xlim(0,2000)+scale_fill_manual(values=region_cols)+ggtitle("dREG Region sizes")+xlab("Size of dREG region (bp)")+ylab("Counts")+facet_grid(~regions,margins=FALSE,scales="free")
ggsave("region_sizes.pdf",width=10)
ggplot(size_dat,aes(x=sizes))+geom_histogram(binwidth=50,fill="gray")+mytheme+xlim(0,2000)+ggtitle("dREG Region sizes (all)")+xlab("Size of dREG region (bp)")+ylab("Counts")
ggsave("all_sizes.pdf")

# Inflation factor
cat("visualisation.r: Inflation!\n")
inflation<-data$inflation
inf_dat<-data.frame(inflation,regions,first_appearance)
ggplot(inf_dat,aes(x=inflation,fill=regions))+geom_histogram(binwidth=0.1)+mytheme+scale_fill_manual(values=region_cols)+ggtitle("Ratio by which smallest contributing region was embiggened")+xlab("Inflation factor")+ylab("Counts")+facet_grid(first_appearance~regions,margins=FALSE,scales="free")+xlim(0,7.5)
ggsave("region_inflation.pdf",width=10)
ggplot(inf_dat,aes(x=inflation))+geom_histogram(binwidth=0.5,fill="gray")+mytheme+ggtitle("Ratio by which smallest contributing region was embiggened (all)")+xlab("Inflation factor")+ylab("Counts")+xlim(0,20)
ggsave("all_inflation.pdf")

# Distance to closest TSS
cat("visualisation.r: Distance to closest TSS!\n")
distance<-data_exp$distance
regions_exp<-data_exp$regions
when<-data_exp$when
dist_dat_pre<-data.frame(distance,regions_exp,when)
dist_dat<-dist_dat_pre[dist_dat_pre$regions_exp!="TSS",]
affected<-dist_dat$distance<time_thresh[dist_dat$when]
thresh<-ifelse(affected,"area cleared by FP treatment",ifelse(dist_dat$regions_exp=="gene_body","unaffected (gene body)","unaffected (non gene)"))
thresh_cols<-c("darkorchid1",region_cols[2],region_cols[3])
dist_dat<-cbind(dist_dat,thresh)
names(dist_dat)<-c("distance","regions","when","thresh")
ggplot(dist_dat,aes(x=distance,fill=thresh))+geom_histogram(binwidth=250)+mytheme+scale_fill_manual(values=thresh_cols)+ggtitle("Distance to closest TSS (either direction)")+xlab("Distance (bp)")+ylab("Counts")+facet_grid(when~regions,margins=FALSE,scale="free",space="free")+xlim(0,25000)
ggsave("region_distance.pdf",width=10)
ggplot(dist_dat,aes(x=distance))+geom_histogram(fill="gray",binwidth=500)+mytheme+ggtitle("Distance to closest TSS (either direction) (gene body and non-gene)")+xlab("Distance")+ylab("Counts")+xlim(0,25000)
ggsave("all_distance.pdf",width=10)

# dREG scores
#cat("make_histograms.r: dREG scores!\n")
#dREG_preds<-read.table(paste(folder,name,".scores.txt",sep=""))
#scores<-dREG_preds[,1]
#thresh<-ifelse(scores>0.8,"over","under")

#score_dat<- data.frame(scores,thresh)
#qplot(scores,data=score_dat,geom="histogram",fill=thresh)+theme_bw()
#ggsave(paste(folder,name,".dREG_scores.pdf",sep=""))

#  Histone modifications! (have to make long form)
counts_per_region<-function(start,stop,counts,norm){
    size<-stop-start
    norm_counts<-counts/(size*norm)
    return(norm_counts)
}
#  Question: how to normalise?!
cat("visualisation.r: Histone modifications!\n")
long_counts<-vector()
for (mod in mods){
    mod_norm<-as.numeric(norms[mod])
    long_counts<-c(long_counts,counts_per_region(data$start,data$end,data[,mod],mod_norm))
}
n_mods<-length(mods)
mod_names<-rep(mods,each=n)
long_regions<-rep(regions,n_mods)
long_first_appearance<-rep(first_appearance,n_mods)
long_IDs<-rep(data$dREG_id,n_mods)
mod_dat<-data.frame(long_IDs,mod_names,long_counts,long_regions,long_first_appearance)
# Remove NaNs (no mapped regions)
mod_dat<-mod_dat[is.finite(mod_dat$long_counts),]
names(mod_dat)<-c("dREG_id","mods","counts","regions","first_appearance")
print(mod_dat[which(mod_dat$counts==max(mod_dat[mod_dat$mods=="H3K4me1",]$counts)),]$dREG_id)
ggplot(mod_dat,aes(x=mods,y=counts,colour=mods,group=mods))+geom_boxplot(notch=TRUE)+facet_grid(~.regions)+scale_y_log10()+mytheme+ylab("Histone marks per base (per million read) : log10")+xlab("Marks and genomic regions")
ggsave("histone_regions.pdf")
