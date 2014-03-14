#! Generate some fake data...

dat<-read.table('~/Tools/UCSC/mm9.chromSizes',as.is=TRUE)
intdat<-dat[c(1:20,22),]
genome_size<-sum(as.numeric(intdat[,2]))
locs<-sample.int(genome_size,size=25000)
chrs<-rep(NA,25000)

# break up the chrs
chros<-c(0)
cum<-0
for (i in seq(nrow(intdat))){
    cum<-cum+intdat[i,2]
    chros<-c(chros,cum)
}

get_chro<-function(chros,loc){
    for (i in seq(length(chros)-1)){
        if((chros[i]<loc)&(loc<chros[i+1])){
            return(i)
        }
    }
}

fake.chrs<-sapply(locs,function(x) intdat[get_chro(chros,x),1])
fake.dat<-data.frame(fake.chrs,locs,locs+1000*rexp(250000,rate=1.5))
write.table(fake.dat,file="fake.bed",row.names=F,col.names=F)
