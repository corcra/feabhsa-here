# Makes histograms: region sizes, distance to TSS, dREG scores
library(ggplot2)
library(grid)
args<-commandArgs(TRUE)
datapath<-args[1]
data<-read.table(datapath,header=T)
n<-nrow(data)

# Timepoints (optionally from argument)
time_points<-c("X0","X2","X5","X12.5","X25","X50")
time_readable<-c("WT","2min","5min","12.5min","25min","50min")
time_thresh<-list("X0"=0,"X2"=0,"X5"=0,"X12.5"=7000,"X25"=30000,"X50"=100000)

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

# Expand regions
data_exp<-data.frame()
for (time in time_points){
    time_data<-data[data[,time]==1,-which(names(data)%in%time_points)]
    en<-nrow(time_data)
    when<-rep(time,en)
    t_d<-data.frame(time_data,when)
    data_exp<-rbind(data_exp,t_d)
}
print(names(data_exp))


# First appearances
first_appearance<-time_points[apply(data[,time_points],1,function(x) which(x==1)[1])]
first_appearance<-factor(first_appearance,time_points)
data<-cbind(data,first_appearance)
cat("visualisation.r: First appearance of regions!\n")
ggplot(data,aes(x=first_appearance,fill=regions))+geom_bar(position="dodge")+scale_fill_manual(values=region_cols)+mytheme
ggsave("first_appearance.pdf",width=10)

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
distance<-data$distance
dist_dat_pre<-data.frame(distance,regions,first_appearance)
dist_dat<-dist_dat_pre[dist_dat_pre$regions!="TSS",]
affected<-dist_dat$distance<time_thresh[dist_dat$first_appearance]
thresh<-ifelse(affected,"area cleared by FP treatment",ifelse(dist_dat$regions=="gene_body","unaffected (gene body)","unaffected (non gene)"))
thresh_cols<-c("darkorchid1",region_cols[2],region_cols[3])
dist_dat<-cbind(dist_dat,thresh)
ggplot(dist_dat,aes(x=distance,fill=thresh))+geom_histogram(binwidth=250)+mytheme+scale_fill_manual(values=thresh_cols)+ggtitle("Distance to closest TSS (either direction)")+xlab("Distance (bp)")+ylab("Counts")+facet_grid(first_appearance~regions,margins=FALSE,scale="free")+xlim(0,25000)
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
