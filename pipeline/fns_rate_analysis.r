cat("rate_analysis_fns.r: Loading functions!\n")

get_region_of_relevance<-function(from,to,genes,timepoint){
# merge timepoints
    both<-merge(from,to,by=1)
        strand_info<-genes[,c(4,6)]
        names(strand_info)<-c("name","strand")
        both<-merge(both,strand_info,by=1)

# look at transitions
        starts<-ifelse(both$strand=="+",floor(both$gene_start.x+both$transition.x),floor(both$gene_end.x-both$transition.y))
        ends<-ifelse(both$strand=="+",ceiling(both$gene_start.x+both$transition.y),ceiling(both$gene_end.x-both$transition.x))
        chr<-both$chr.x
        name<-both$name
        filler<-rep("NA",length(name))

# combinate
        region<-data.frame(chr,starts,ends,name,filler)
        region<-subset(region,ends>starts)

        filename<-paste0(timepoint,"_region_of_relevance.bed.temp")
        cat("Saving region of relevance to",filename,"\n")
        write.table(region,file=filename,row.names=F,col.names=F,quote=F)
        cat("Sorting to",timepoint,"_region_of_relevance.bed!\n")
        system(paste0("sort-bed ",filename," > ",timepoint,"_region_of_relevance.bed"))
}

plot_transitions<-function(five,twelve,twentyfive,fifty){
    time<-c(rep("5min",nrow(five)),rep("12.5min",nrow(twelve)),rep("25min",nrow(twentyfive)),rep("50min",nrow(fifty)))
        time<-factor(time,c("5min","12.5min","25min","50min"))
        all<-rbind(five,twelve,twentyfive,fifty)
        all<-data.frame(all,time)
        ggplot(all,aes(x=transition/1000,fill=time))+geom_histogram(binwidth=2)+facet_grid(time~.,scales="free")+xlab("transition point (kb)")+theme_bw()+scale_fill_manual(values=rainbow(4,v=0.8,s=0.5))
        cat("Plotting transitions.\n")
        ggsave("transitions.pdf",width=6,height=4)
}

consistent_over_replicates<-function(data_1,data_2,time){
    both<-merge(data_1,data_2,by=2)
        sd<-sqrt(var(both$transition.x-both$transition.y))
#sd<-10000000
        if (sd==0){
            cat("Transitions are identical!\n")
                name<-data_1$name
                chr<-data_1$chr
                gene_start<-data_1$gene_start
                gene_end<-data_1$gene_end
                transition<-data_1$transition
                cons<-data.frame(name,chr,gene_start,gene_end,transition)
                return(cons)
        } else{
#        pdf(paste0("transition_reproducibility_",time,".pdf"))
#        hist(both$transition.x-both$transition.y,breaks=100,xlab="Difference in transition points between replicates (bp)",main=paste0("Reproducibility of transition points (",time," min)"))
#        abline(v=sd,col="red")
#        abline(v=-sd,col="red")
#        graphics.off()
            cons<-subset(both,abs(transition.x-transition.y)<2*sd)
                lost<-subset(both,abs(transition.x-transition.y)>=2*sd)
                name<-cons$name
                transition<-0.5*(cons$transition.x+cons$transition.y)
                gene_start<-cons$gene_start.x
                gene_end<-cons$gene_end.x
                chr<-cons$chr.x
                cat(nrow(lost),"genes lost due to irreproducibility. (",time,")\n")
                cons<-data.frame(name,chr,gene_start,gene_end,transition)
                return(cons)
        }
}

get_rate<-function(from,to,timediff){
# note, the name field should be the second one...
    both<-merge(from,to,by=1)
        name<-both$name
        chr<-both$chr.x
        start<-both$gene_start.x
        end<-both$gene_end.x
        rate<-(both$transition.y-both$transition.x)/(timediff*1000)
        name<-both$name
        df<-data.frame(chr,start,end,name,rate)
        df_pos<-subset(df,rate>0)
        return(df_pos)
}

# taken from hojoong's script
roz<-function(d)
{
    r<-qnorm(rank(d+runif(length(d),0,0.1))/(length(d)+1))
        return(r)
}

# also taken from hojoong's script
myplot<-function(x,y,xlabel='',ylabel='',arrow=TRUE,labcol='black')
{
    plot(c(0,1),c(0,1),type='n',xlim=c(-4,4),ylim=c(-4,4),ylab=ylabel,xlab=xlabel,col.lab=labcol,xaxs='i',yaxs='i',main='',yaxt='n',xaxt='n',cex.lab=0.8)
    points(x,y,pch=16,cex=0.3)
    axis(side = 1, labels = NA, tck = -0.02, at=seq(-4,4,2),cex.axis=0.8)
    axis(side = 2, labels = NA, tck = -0.02, at=seq(-4,4,2),cex.axis=0.8)
    pcaresult<-princomp(cbind(x,y))
    if(pcaresult$loadings[1,1]>0){
        arrowdir<--1
    } else{
        arrowdir<-1
    }
    r<-(pcaresult$sdev[1]^2/(pcaresult$sdev[1]^2+pcaresult$sdev[2]^2)-0.5)*2
    if(arrow){
        arrows(pcaresult$center[1],pcaresult$center[2],pcaresult$center[1]-arrowdir*pcaresult$loadings[1,1]*pcaresult$sdev[1]*r*4,pcaresult$center[2]-arrowdir*pcaresult$loadings[2,1]*pcaresult$sdev[1]*r*4,col='red',length=0.05,lwd=1.5)
    }
    text(2,3.5,sprintf("r = %1.2f",r),cex=0.8)
#p<-coef(summary(lm(y~x)))[2,4]
#if(p<0.001) text(2,-3.5,sprintf("p < 0.001"),cex=0.7)
#else text(2,-3.5,sprintf("p = %1.3f",p),cex=0.7)
}

add_covars<-function(data,n_gb,int1len,n_exon,CpG_gb,H3K79me2){
    data<-merge(data,n_gb,by.x=4,by.y=1)
        data<-merge(data,int1len,by=1)
        data<-merge(data,n_exon,by=1)
        data<-merge(data,CpG_gb,by=1)
        data<-merge(data,H3K79me2_gb,by=1)

        roz.rate<-roz(data$rate)
        roz.n_gb<-roz(data$n_gb)
        roz.int1len<-roz(data$int1len)
        roz.n_exon<-roz(data$n_exon)
        roz.CpG<-roz(data$CpG)
        roz.H3K79me2<-roz(data$H3K79me2)

        data<-data.frame(data,roz.rate,roz.n_gb,roz.int1len,roz.n_exon,roz.CpG,roz.H3K79me2)
        return(data)
}

rateplot<-function(data,title){
    ggplot(data,aes(x=n_gb,y=rate,group=n_gb))+geom_boxplot()+theme_bw()+xlab("# dREG hits")+ylab("rate (kb/min)")+ggtitle(title)
}
