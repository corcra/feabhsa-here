library(ggplot2)
library(grid)
library(bigWig)

preprocess<-TRUE
hw<-3000
stepp<-1

get_metaplot<-function(bedfile,bw,hw,stepp,region_name){
    cat("Metaplotting! (",region_name,")\n")
    res<-meta.subsample(bedfile,bw,NULL,halfWindow=hw,step=stepp)
    lv<-colMeans(res[[1]])
    uq<-res[[2]]
    lq<-res[[3]]
    d<-seq(-length(lv)/2,length(lv)/2-1)
    region<-rep(region_name,length(lv))
    region<-factor(region_name,c("gene_start","gene_body","non_gene"))
    dl<-data.frame(d,lv,uq,lq,region_name)
    return(dl)
}

plot_metaplot<-function(dl,dataname){
    cat("Plotting!(",dataname,")!\n")
    mytheme<-theme(plot.margin=unit(c(2,2,2,2),"cm"),axis.title.x=element_text(vjust=-0.2),axis.title.y=element_text(angle=90,vjust=0.2),plot.title=element_text(size=15,vjust=1.4))+theme_bw()
    region_cols<-c("firebrick1","orange1","turquoise3","gray")

    browser()
    ggplot(dl,aes(x=d,y=lv,colour=region_name,group=region_name))+geom_point(cex=0.5,show_guide=FALSE)+facet_grid(~region_name)+mytheme+xlab("Distance to centre (bp)")+ylab("Average ChIP-seq signal")+scale_colour_manual(values=region_cols)+ggtitle(paste0("Metagene plot (",dataname,")"))+geom_ribbon(aes(ymin=lq,ymax=uq,alpha=0.5),show_guide=FALSE)
    ggsave(paste0(dataname,".pdf"))
}

if(!preprocess){
    
    H3K4me1_data<-read.table(file="H3K4me1_metaplot.txt",header=TRUE,as.is=TRUE)
    H3K4me1_data$region_name<-factor(H3K4me1_data$region_name,c("gene_start","gene_body","non_gene"))
    plot_metaplot(H3K4me1_data,"H3K4me1")

    H3K4me3_data<-read.table(file="H3K4me3_metaplot.txt",header=TRUE,as.is=TRUE)
    H3K4me3_data$region_name<-factor(H3K4me3_data$region_name,c("gene_start","gene_body","non_gene"))
    plot_metaplot(H3K4me3_data,"H3K4me3")

    H3K27ac_data<-read.table(file="H3K27ac_metaplot.txt",header=TRUE,as.is=TRUE)
    H3K27ac_data$region_name<-factor(H3K27ac_data$region_name,c("gene_start","gene_body","non_gene"))
    plot_metaplot(H3K27ac_data,"H3K27ac")

    H3_data<-read.table(file="H3_metaplot.txt",header=TRUE,as.is=TRUE)
    H3_data$region_name<-factor(H3_data$region_name,c("gene_start","gene_body","non_gene"))
    plot_metaplot(H3_data,"H3")
} else{
    gene_start<-read.table('/Users/stephanie/ll/results/FP/dREG_gs.bed')
    gene_body<-read.table('/Users/stephanie/ll/results/FP/dREG_gb.bed')
    non_gene<-read.table('/Users/stephanie/ll/results/FP/dREG_ng.bed')

    cat("H3K4me1!\n")
    H3K4me1<-load.bigWig('/Users/stephanie/ll/data/ChIPseq/bw/H3K4me1_Meiss.sorted.bw')
    gs_dat<-get_metaplot(gene_start,H3K4me1,hw,stepp,"gene_start")
    gb_dat<-get_metaplot(gene_start,H3K4me1,hw,stepp,"gene_body")
    ng_dat<-get_metaplot(gene_start,H3K4me1,hw,stepp,"non_gene")
    H3K4me1_data<-rbind(gs_dat,gb_dat,ng_dat)
    write.table(H3K4me1_data,row.names=F,col.names=T,quote=F,file="H3K4me1_metaplot.txt")
    plot_metaplot(H3K4me1_data,"H3K4me1")

    cat("H3K4me3!\n")
    H3K4me3<-load.bigWig('/Users/stephanie/ll/data/ChIPseq/bw/H3K4me3_Marson.sorted.bw')
    gs_dat<-get_metaplot(gene_start,H3K4me3,hw,stepp,"gene_start")
    gb_dat<-get_metaplot(gene_start,H3K4me3,hw,stepp,"gene_body")
    ng_dat<-get_metaplot(gene_start,H3K4me3,hw,stepp,"non_gene")
    H3K4me3_data<-rbind(gs_dat,gb_dat,ng_dat)
    write.table(H3K4me3_data,row.names=F,col.names=T,quote=F,file="H3K4me3_metaplot.txt")

    cat("H3K27ac!\n")
    H3K427ac<-load.bigWig('/Users/stephanie/ll/data/ChIPseq/bw/H3K27Ac.sorted.bw')
    gs_dat<-get_metaplot(gene_start,H3K27ac,hw,stepp,"gene_start")
    gb_dat<-get_metaplot(gene_start,H3K27ac,hw,stepp,"gene_body")
    ng_dat<-get_metaplot(gene_start,H3K27ac,hw,stepp,"non_gene")
    H3K27ac_data<-rbind(gs_dat,gb_dat,ng_dat)
    write.table(H3K27ac_data,row.names=F,col.names=T,quote=F,file="H3K27ac_metaplot.txt")

    cat("H3!\n")
    H3<-load.bigWig('/Users/stephanie/ll/data/ChIPseq/bw/H3_align.sorted.bw')
    gs_dat<-get_metaplot(gene_start,H3,hw,stepp,"gene_start")
    gb_dat<-get_metaplot(gene_start,H3,hw,stepp,"gene_body")
    ng_dat<-get_metaplot(gene_start,H3,hw,stepp,"non_gene")
    H3_data<-rbind(gs_dat,gb_dat,ng_dat)
    write.table(H3_data,row.names=F,col.names=T,quote=F,file="H3_metaplot.txt")
}

