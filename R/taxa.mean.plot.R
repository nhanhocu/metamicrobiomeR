#' Plot mean taxa abundance
#'
#' This function visualize mean relative abundance by group as stacked plots.
#' @param tabmean table of mean abundance generated from taxa.meansdn.
#' @param sumvar variable to be plotted. Options are c("taxa","path"). Default is "taxa"
#' @param tax.select list of selected taxa/pathways to be plotted. Default is "none" or plot all taxa/pathways.
#' @param tax.lev taxa level to be visualized. Options are from "l2" (phylum) to "l7" (species). Default is "l2". If sumvar="path", all pathways will be visualized.
#' @param comvar main variable for comparison.
#' @param groupvar variable for stratifying.
#' @param mean.filter mean abundance filtering threshold (only plot those with mean abundance>threshold and plot all those with mean abundance <threshold as "other").
#' @param pallete.by.phylum whether each pallete of color for each phylum. Default is FALSE (plot distinc colors).
#' @param show.taxname whether show "full" taxa name or "short" name. Default is "full".
#' @param legend.position position of legend. Options are c("right", "left","bottom","top","none") as in ggplot2. Default is "right".
#' @return a list of ggplot2 object and list of taxa/pathways plotted (those with mean abundance >mean.filter).
#' @keywords mean taxa abundance plot.
#' @export
#' @examples
#' #Load summary tables of bacterial taxa relative abundance from Bangladesh data
#' data(taxtab.rm7)
#' taxa.meansdn.rm<-taxa.meansdn(taxtab=taxtab.rm[[5]],sumvar="bf",groupvar="age.sample")
#' p.bf.l2<-taxa.mean.plot(tabmean=taxa.meansdn.rm,tax.lev="l2", comvar="bf", groupvar="age.sample",mean.filter=0.005)
#' p.bf.l2$p

taxa.mean.plot<-function(tabmean,sumvar="taxa",tax.select="none",tax.lev="l2", comvar, groupvar,mean.filter=0.005, pallete.by.phylum=FALSE, show.taxname="full",legend.position="right",xlab="Chronological age (month)",ylab="Relative abundance"){
  #sapply(c("ggplot2","reshape2","RColorBrewer","plyr"), require, character.only = TRUE)
  taxname1<-c(colnames(tabmean)[grep("mean",colnames(tabmean))])
  tabm<-apply(tabmean[,taxname1],2,mean)
  tabm<-tabm[tabm>=mean.filter]
  if (sumvar=="taxa"){
    l2<-taxname1[-grep("c__",taxname1)]
    l3<-taxname1[-grep("o__",taxname1)]
    l4<-taxname1[-grep("f__",taxname1)]
    l5<-taxname1[-grep("g__",taxname1)]
    l6<-taxname1
    taxnamel<-list(l2=l2,l3=l3,l4=l4,l5=l5,l6=l6)
    tabmeanl<-list(l2=tabmean[,colnames(tabmean) %in% c(comvar, groupvar, l2)],l3=tabmean[,colnames(tabmean) %in% c(comvar, groupvar,l3)],l4=tabmean[,colnames(tabmean) %in% c(comvar, groupvar,l4)],l5=tabmean[,colnames(tabmean) %in% c(comvar, groupvar,l5)],l6=tabmean[,colnames(tabmean) %in% c(comvar, groupvar,l6)])
    tab<-tabmeanl[[tax.lev]]
    taxuse<-taxnamel[[tax.lev]]
    if (tax.lev=="l2"){
      taxuse<-taxuse
    }
    if (tax.lev=="l7"){
      stringt<-"s__"
      taxuse<-taxuse[grep(".s__",taxuse)]
    }
    if (tax.lev=="l6"){
      stringt<-"g__"
      taxuse<-taxuse[grep(".g__",taxuse)]
    }
    if (tax.lev=="l5"){
      stringt<-"f__"
      taxuse<-taxuse[grep(".f__",taxuse)]
      taxs<-sub(".*f__", "",taxuse)
      taxuse<-taxuse[!taxs %in% "_mean"]
    }
    if (tax.lev=="l4"){
      stringt<-"o__"
      taxuse<-taxuse[grep(".o__",taxuse)]
      taxs<-sub(".*o__", "",taxuse)
      taxuse<-taxuse[!taxs %in% "_mean"]
    }
  }
  if (sumvar=="path"){
    taxuse=taxname1
  }
  if (tax.select=="none"){
    taxuse<-taxuse[taxuse %in% names(tabm)]
  }
  if (tax.select!="none"){
    taxuse<-taxuse[taxuse %in% paste(tax.select,"_mean",sep="")]
  }
  taxuse.rm<-gsub("_mean","",taxuse)
  tab[,"others_mean"]<-1-apply(tab[,taxuse],1,sum)
  tab.l <- reshape(tab,
                   varying = c("others_mean",taxuse),
                   v.names = "rel_abund",
                   timevar = "taxa",
                   times = c("others_mean",taxuse),
                   new.row.names = 1:(nrow(tab)*length(c("others_mean",taxuse))),
                   direction = "long")
  tab.l<-as.data.frame(tab.l)
  tab.l$taxa<-gsub("_mean","",tab.l$taxa)
  tab.l<-tab.l[,c(comvar,groupvar,"taxa","rel_abund")]
  tab.l$taxa<-gsub("k__bacteria.p__","",tab.l$taxa)
  tab.l$taxa<-as.factor(as.character(tab.l$taxa))
  tab.l$taxa<-factor(tab.l$taxa,levels=c(levels(tab.l$taxa)[!levels(tab.l$taxa) %in% "others"],"others"))
  tab.l$phylum<-sub(".c__.*","",tab.l$taxa)
  tab.l<-tab.l[order(tab.l$taxa),]
  if (tax.lev!="l2"){
    tab.l$taxas<-sub(paste(".*",stringt,sep=""), "",tab.l$taxa)
  }
  if (tax.lev=="l2"){
    tab.l$taxas<-tab.l$taxa
  }
  tab.l$taxas<-paste0(toupper(substr(as.character(tab.l$taxas), 1, 1)), substr(as.character(tab.l$taxas), 2, nchar(as.character(tab.l$taxas))))
  tab.l$taxas<-factor(tab.l$taxas,levels=unique(tab.l$taxas))
  #get color vector
  n <- nlevels(tab.l$taxa)
  qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  tab.l$phylum<-factor(tab.l$phylum,levels=c(unique(tab.l$phylum)[!unique(tab.l$phylum) %in% "others"],"others"))
  if (pallete.by.phylum==TRUE){
    tab.l$colpal<-plyr::mapvalues(tab.l$phylum, from=c("actinobacteria","bacteroidetes","firmicutes","proteobacteria","others" ),to=c("Greens","Purples","Reds","Blues","Greys"))
    ngroup<-length(unique(tab.l[,comvar]))*length(unique(tab.l[,groupvar]))
    colset<-list()
    for (i in 1: nlevels(tab.l$colpal)){
      #get number of unique taxa in a phylum (number of colors in a pallette)
      coltab<-length(unique(tab.l$taxa[tab.l$colpal %in% levels(tab.l$colpal)[i]]))
      colset[[i]]<-rev(RColorBrewer::brewer.pal(9,levels(tab.l$colpal)[i]))[1:coltab]
    }
    tab.l$col<-plyr::mapvalues(tab.l$taxa,from=levels(tab.l$taxa),to=c(unlist(colset)))
    col_vector<-c(unlist(colset))
  }
  # display taxa names
  if(show.taxname=="short"){
    tab.l$taxa<-tab.l$taxas
  }
  if (is.numeric(tab.l[,groupvar])){
    #stacked plot by numeric variable (e.g.age)
    p<-ggplot2::ggplot(tab.l, ggplot2::aes(x=get(as.character(groupvar)),y=rel_abund))+
      ggplot2::geom_area(ggplot2::aes(fill=taxa))+
      ggplot2::scale_fill_manual(values=col_vector)+ ggplot2::xlab(xlab)+ggplot2::ylab(ylab)+
      ggplot2::labs(fill='')+
      ggplot2::theme(legend.position = legend.position)+
      ggplot2::theme(legend.text = ggplot2::element_text(colour="black", size = 8))+
      ggplot2::scale_x_continuous(breaks=seq(from=0,to=24,by=3),
                         labels=seq(from=0,to=24,by=3))+
      ggplot2::theme(legend.key.size = ggplot2::unit(0.5, "cm"),
            axis.line = ggplot2::element_line(colour = "black"),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            panel.border = ggplot2::element_blank(),
            panel.background = ggplot2::element_blank(),
            strip.background =ggplot2::element_rect(fill="white"))+
      ggplot2::guides(fill=ggplot2::guide_legend(ncol=1)) +
      ggplot2::facet_wrap(~get(as.character(comvar)), ncol = 1)
  }
  else{
    #barplot by group variable
    p<-ggplot2::ggplot(tab.l, ggplot2::aes(x = get(as.character(comvar)), y = rel_abund, fill = taxa)) +
      ggplot2::geom_bar(stat = "identity")+
      ggplot2::scale_fill_manual(values=col_vector)+
      ggplot2::xlab(comvar)+
      ggplot2::ylab(ylab)+
      ggplot2::labs(fill='')+
      ggplot2::guides(fill=ggplot2::guide_legend(ncol=1))+
      ggplot2::theme(legend.text = ggplot2::element_text(colour="black", size = 8))+
      ggplot2::theme(legend.position = legend.position)+
      ggplot2::theme(legend.key.size = ggplot2::unit(0.5, "cm"),
            axis.line = ggplot2::element_line(colour = "black"),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            panel.border = ggplot2::element_blank(),
            panel.background = ggplot2::element_blank(),
            strip.background =ggplot2::element_rect(fill="white"))+
      ggplot2::guides(fill=ggplot2::guide_legend(ncol=1)) +
      ggplot2::facet_wrap(~get(as.character(groupvar)), ncol = 1)
  }
  return(list(p=p,taxuse.rm=taxuse.rm))
}
