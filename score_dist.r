library(ggplot2)

# this one weird value will appear in your data for no reason
weirdval=0.179973493896357
#weirdval=-0.229395386677444

timedata<-function(time){
     gs<-read.table(paste0("FP_",time,"min.gs.predictions.bedGraph.gz"))
     gb<-read.table(paste0("FP_",time,"min.gb.predictions.bedGraph.gz"))
     ng<-read.table(paste0("FP_",time,"min.ng.predictions.bedGraph.gz"))
     region<-c(rep("gene_start",nrow(gs)),rep("gene_body",nrow(gb)),rep("non_gene",nrow(ng)))
     region<-factor(region,c("gene_start","gene_body","non_gene"))
     score<-c(gs[,4],gb[,4],ng[,4])
     time<-rep(paste0(time,"min"),length(score))
     td<-data.frame(score,region,time)
     return(td)
}

alldata<-timedata(0)
cat("Loaded 0min\n")
alldata<-rbind(alldata,timedata(2))
cat("Loaded 2min\n")
alldata<-rbind(alldata,timedata(5))
cat("Loaded 5min\n")
alldata<-rbind(alldata,timedata(12.5))
cat("Loaded 12.5min\n")
alldata<-rbind(alldata,timedata(25))
cat("Loaded 25min\n")
alldata<-rbind(alldata,timedata(50))
cat("Loaded 50min\n")

# Check for that weird value, delete it...
print(mean(alldata$score==weirdval))
alldata<-alldata[alldata$score!=weirdval,]

region_cols<-c("firebrick1","orange1","turquoise3","gray")
ggplot(alldata,aes(x=score,fill=region))+geom_histogram(binwidth=0.01)+facet_wrap(time~region,ncol=3,scales="free")+scale_fill_manual(values=region_cols)+theme_bw()
ggsave("dREG_scores_hist.pdf")
