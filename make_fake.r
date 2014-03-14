#! Generate some fake data...

N.fake<-25000
dat<-read.table('~/Tools/UCSC/mm9.chromSizes',as.is=TRUE)
intdat<-dat[c(1:20,22),]
genome_size<-sum(as.numeric(intdat[,2]))
locs<-sample.int(genome_size,size=N.fake)
chrs<-rep(NA,N.fake)

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
fake.locs<-sapply(locs,function(x) x-chros[get_chro(chros,x)])
fake.ids<-paste("fake_id_",seq(N.fake),sep="")
header<-c("chr","start","end","dREG_id","0","2","5","12.5","25","50","inflation","when_smallest")
fake.X0<-rbinom(N.fake,1,0.5)
fake.X2<-rbinom(N.fake,1,0.5)
fake.X5<-rbinom(N.fake,1,0.5)
fake.X12.5<-rbinom(N.fake,1,0.5)
fake.X25<-rbinom(N.fake,1,0.5)
fake.inflation<-5*rexp(N.fake,rate=4)
fake.X50<-ifelse((fake.X0+fake.X2+fake.X5+fake.X12.5+fake.X25)==0,1,rbinom(N.fake,1,0.5))
fake.dat<-data.frame(fake.chrs,fake.locs,fake.locs+1+ceiling(1000*rexp(N.fake,rate=1.5)),fake.ids,fake.X0,fake.X2,fake.X5,fake.X12.5,fake.X25,fake.X50,fake.inflation,rep("X0",N.fake))
names(fake.dat)<-header
write.table(fake.dat,file="fake_uniq.bed",row.names=F,col.names=T,quote=F,sep="\t")
