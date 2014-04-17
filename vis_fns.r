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
