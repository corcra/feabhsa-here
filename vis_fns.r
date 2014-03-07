get_consistent<-function(df,n_viol){
    n<-nrow(df)
#    consistent<-rep(TRUE,n)
    violations<-rep(0,n)
    one_seen<-rep(FALSE,n)
    for (col in seq(2,ncol(df))){
#        consistent<-consistent&(df[,col]>=df[,col-1]) # strictest condition
        one_seen<-one_seen|df[,col]==1
        violations<-violations+1*((one_seen)*(df[,col]==0))
#        violations<-violations+(1*(df[,col]<df[,col-1]))
    }
    consistent<-violations<n_viol
    return(consistent)
}
