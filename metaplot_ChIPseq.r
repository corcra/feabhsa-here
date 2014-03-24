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
    lq<-res[[2]]
    d<-seq(-length(lv)/2,length(lv)/2-1)
    region<-rep(region_name,length(lv))
    dl<-data.frame(d,lv,uq,lq,region_name)
    return(dl)
}

plot_metaplot<-function(dl,dataname){
    cat("Plotting!(",dataname,")!\n")
    mytheme<-theme(plot.margin=unit(c(2,2,2,2),"cm"),axis.title.x=element_text(vjust=-0.2),axis.title.y=element_text(angle=90,vjust=0.2),plot.title=element_text(size=15,vjust=1.4))+theme_bw()
    region_cols<-c("firebrick1","orange1","turquoise3","gray")

    browser()
    ggplot(dl,aes(x=d,y=lv,colour=region_name,group=region_name))+geom_point(cex=0.5,show_guide=FALSE)+facet_grid(~region_name)+mytheme+xlab("Distance to centre (bp)")+ylab("Average ChIP-seq signal")+scale_colour_manual(values=region_name_cols)+geom_line(y=0,linetype="dashed",show_guide=FALSE,colour="black")+ggtitle(paste0("Metagene plot (",dataname,")"))+geom_ribbon(aes(ymin=lq,ymax=uq,alpha=0.5))
    ggsave(paste0(dataname,".pdf"))
}

if(!preprocess){
    # this needs more tweaking
    basepath<-'/Users/stephanie/ll/data/ChIPseq/bw/'
    gs<-read.table(paste0(basepath,name,'_genestart.txt'),header=TRUE)
    gb<-read.table(paste0(basepath,name,'_genebody.txt'),header=TRUE)
    ng<-read.table(paste0(basepath,name,'_nongene.txt'),header=TRUE)

    region<-rep("gene_start",nrow(gs))
    gs<-cbind(gs,region)
    region<-rep("gene_body",nrow(gb))
    gb<-cbind(gb,region)
    region<-rep("non_gene",nrow(ng))
    ng<-cbind(ng,region)
    dl<-rbind(gs,gb,ng)

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
    plot_metaplot(H3K4me1_data,"H3K4me1")
    browser()

    cat("H3K4me3!\n")
    H3K4me3<-load.bigWig('/Users/stephanie/ll/data/ChIPseq/bw/H3K4me3_Marson.sorted.bw')
    gs_dat<-get_metaplot(gene_start,H3K4me3,hw,stepp,"gene_start")
    gb_dat<-get_metaplot(gene_start,H3K4me3,hw,stepp,"gene_body")
    ng_dat<-get_metaplot(gene_start,H3K4me3,hw,stepp,"non_gene")
    H3K4me3_data<-rbind(gs_dat,gb_dat,ng_dat)

    cat("H3K27ac!\n")
    H3K427ac<-load.bigWig('/Users/stephanie/ll/data/ChIPseq/bw/H3K27Ac.sorted.bw')
    gs_dat<-get_metaplot(gene_start,H3K27ac,hw,stepp,"gene_start")
    gb_dat<-get_metaplot(gene_start,H3K27ac,hw,stepp,"gene_body")
    ng_dat<-get_metaplot(gene_start,H3K27ac,hw,stepp,"non_gene")
    H3K27ac_data<-rbind(gs_dat,gb_dat,ng_dat)

    cat("H3!\n")
    H3<-load.bigWig('/Users/stephanie/ll/data/ChIPseq/bw/H3_align.sorted.bw')
    gs_dat<-get_metaplot(gene_start,H3,hw,stepp,"gene_start")
    gb_dat<-get_metaplot(gene_start,H3,hw,stepp,"gene_body")
    ng_dat<-get_metaplot(gene_start,H3,hw,stepp,"non_gene")
    H3_data<-rbind(gs_dat,gb_dat,ng_dat)
}

