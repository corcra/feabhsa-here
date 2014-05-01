# Perform quality control (of arbitrary stlye)

args<-commandArgs(TRUE)
# Options
consistent<-TRUE

# Data?
suffix<-"40m_SVM"           # 40m_SVM or full_SVM
#suffix<-"full_SVM"           # 40m_SVM or full_SVM

# Paths
inpath<-args[1]
outpath<-args[2]
cat("QC.r : performing QC on ",inpath,", outputting ",outpath,"\n")

#inpath<-paste0('/Users/stephanie/ll/results/',suffix,'/FP_dREG_regions_marked.bed.gz')
#outpath<-paste0('/Users/stephanie/ll/results/',suffix,'/FP_dREG_regions_filtered.bed')

# Set up some vars
time_points<-c("X0","X2","X5","X12.5","X25","X50")

# Functions and stuff
get_violations<-function(df){
    n<-nrow(df)
    violations<-rep(0,n)
    one_seen<-(df[,1]==1)
    for (col in seq(2,ncol(df))){
        one_seen<-one_seen|df[,col]==1
        violations<-violations+1*((one_seen)*(df[,col]==0))
    }
    return(violations)
}

# Get the data!
data<-read.table(inpath,header=T,na.strings="NAN")

# Apply filters, etc.
if (consistent){
    print("Restricting to consistent hits!")
    violations<-get_violations(data[,time_points])
    consistent<-violations<=1
    data_cons<-subset(data,consistent==TRUE)
}

write.table(data_cons,file=outpath,quote=F,row.names=F,sep="\t")
