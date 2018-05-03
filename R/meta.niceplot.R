#' Nice meta-analysis plots.
#'
#' This function displays meta-analysis results of relative abundance as a nice combined heatmap and forest plot.
#' @param metadat output data from metatab.show.
#' @param sumtype Either "taxa" for taxa and "path" for pathway.
#' @param level "main" for main level such as phylum or "sub" for higher level such as species. Default is "main".
#' @param p name of variable for p-values
#' @param p.adjust name of variable for multiple testing adjusted p-values
#' @param phyla.col type of color for main level (phylum). Options are "select" (default) or "rainbow".
#' @param leg.key.size legdend key size for heatmap.
#' @param leg.text.size legend text size for heatmap.
#' @param heat.text.x.size heatmap x label text size.
#' @param heat.text.x.angle heatmap x label text angle.
#' @param forest.axis.text.y forest plot y label text size.
#' @param forest.axis.text.x forest plot x label text size.
#' @return combined heatmap forest plot.
#' @keywords meta-analysis heatmap forest plot.
#' @export
#' @examples
#' #Load saved results of four studies for the comparison of bacterial taxa relative abundance between genders adjusted for breastfeeding and infant age at sample collection
#' data(taxacom.rm.sex.adjustbfage)
#' data(taxacom.ha.sex.adjustbfage)
#' data(taxacom6.zi.usbmk.sex.adjustbfage)
#' data(taxacom6.unc.sex.adjustedbfage)
#' taxacom6.zi.rm.sex.adjustbfage$study<-"Subramanian et al 2014 (Bangladesh)"
#' taxacom6.zi.rm.sex.adjustbfage$pop<-"Bangladesh"
#' taxacom.zi.ha.sex.adjustbfage$study<-"Bender et al 2016 (Haiti)"
#' taxacom.zi.ha.sex.adjustbfage$pop<-"Haiti"
#' taxacom6.zi.usbmk.sex.adjustbfage$study<-"Pannaraj et al 2017 (USA(CA_FL))"
#' taxacom6.zi.usbmk.sex.adjustbfage$pop<-"USA(CA_FL)"
#' taxacom6.zi.unc.sex.adjustedbfage$study<-"Thompson et al 2015 (USA(NC))"
#' taxacom6.zi.unc.sex.adjustedbfage$pop<-"USA(NC)"
#' tabsex4<-plyr::rbind.fill(taxacom6.zi.rm.sex.adjustbfage,taxacom.zi.ha.sex.adjustbfage,taxacom6.zi.usbmk.sex.adjustbfage,taxacom6.zi.unc.sex.adjustedbfage)
#' #Meta-analysis (take time to run)
#' metab.sex<-meta.taxa(taxcomdat=tabsex4,summary.measure="RR",pool.var="id",studylab="study",backtransform=FALSE,percent.meta=0.5,p.adjust.method="fdr")
#' #nice plot phylum level
#' metadat<-metatab.show(metatab=metab.sex$random,com.pooled.tab=tabsex4,tax.lev="l2",showvar="genderMale",p.cutoff.type="p", p.cutoff=1,display="data")
#' meta.niceplot(metadat=metadat,sumtype="taxa",level="main",p="p",p.adjust="p.adjust",phyla.col="rainbow",leg.key.size=1,leg.text.size=8,heat.text.x.size=7,heat.text.x.angle=0,forest.axis.text.y=8,forest.axis.text.x=7)
#' #nice plot family level
#' metadat<-metatab.show(metatab=metab.sex$random,com.pooled.tab=tabsex4,tax.lev="l5",showvar="genderMale",p.cutoff.type="p", p.cutoff=1,display="data")
#' meta.niceplot(metadat=metadat,sumtype="taxa",level="sub",p="p",p.adjust="p.adjust",phyla.col="rainbow",leg.key.size=1,leg.text.size=8,heat.text.x.size=7,forest.axis.text.y=8,forest.axis.text.x=7)


meta.niceplot<-function(metadat,sumtype="taxa",level="main",p,p.adjust,phyla.col=c("select","rainbow"),leg.key.size=1,leg.text.size=8,heat.text.x.size=8,heat.text.x.angle=0,forest.axis.text.y=8,forest.axis.text.x=8){
  #require(ggplot2);require(gridExtra);require("gplots");require(reshape2); require(gdata)
  #heatmap
  test<-metadat$taxsig.all
  test$taxa<-test$id
  test$taxa<-gsub("k__bacteria.p__","",test$taxa)
  if (sumtype=="taxa"){
    test$esticat<-cut(test$estimate, breaks=c(-Inf, -1,-0.5,-0.1,0,0.1,0.5,1, Inf),
                      labels=c("<-1", "[-1,-0.5)","[-0.5,-0.1)","[-0.1,0)", "[0,0.1)", "[01,0.5)", "[0.5,1)", ">=1"))
    test$esticol<-plyr::mapvalues(test$esticat,from=c("<-1", "[-1,-0.5)","[-0.5,-0.1)","[-0.1,0)", "[0,0.1)", "[01,0.5)", "[0.5,1)", ">=1"),
                            to=c("#006d2c", "#2ca25f", "#66c2a4","#b2e2e2", "#fecc5c","#fd8d3c", "#f03b20","#bd0026"))
  }
  if (sumtype=="path"){
    test$esticat<-cut(test$estimate, breaks=c(-Inf, -0.5,-0.1,-0.05,0,0.05,0.1,0.5, Inf),
                      labels=c("<-0.5)", "[-0.5,-0.1)","[-0.1,-0.05)","[-0.05,0)","[0,0.05)","[0.05,0.1)", "[0.1,0.5)", ">=0.5"))
    test$esticol<-plyr::mapvalues(test$esticat,from=c("<-0.5)", "[-0.5,-0.1)","[-0.1,-0.05)","[-0.05,0)","[0,0.05)","[0.05,0.1)", "[0.1,0.5)", ">=0.5"),
                            to=c("#006d2c", "#2ca25f", "#66c2a4","#b2e2e2", "#fecc5c","#fd8d3c", "#f03b20","#bd0026"))
  }
  test$esticat<-gdata::drop.levels(test$esticat, reorder=FALSE)
  test$esticol<-gdata::drop.levels(test$esticol, reorder=FALSE)
  poplev<-levels(factor(test$pop))
  test$pop[is.na(test$pop)]<-"Pooled"
  test$pop<-factor(test$pop,levels=c(poplev[poplev!="Pooled"],"Pooled"))
  if (sumtype=="taxa"){
    if (level=="main"){
      test<-test[order(test$taxa,decreasing = FALSE),]
      test$taxas<-factor(test$taxa,levels=unique(test$taxa))
      test$plotvar<-test$taxas
    }
    if (level=="sub"){
      test$taxas<-sub(".c__.*f__", " ",as.character(test$taxa))
      test<-test[order(test$taxa,decreasing = FALSE),]
      test$taxas<-factor(test$taxas,levels=unique(test$taxas))
      test$taxa<-factor(test$taxa,levels=unique(test$taxa))
      test$plotvar<-test$taxas
    }
  }
  if (sumtype=="path"){
    test<-test[order(test$taxa,decreasing = FALSE),]
    test$taxas<-factor(test$taxa,levels=unique(test$taxa))
    test$plotvar<-test$taxas
  }
  test$study[is.na(test$study)]<-"Meta_analysis"
  test$pdot<-cut(test$p, breaks=c(0, 0.0001,0.05, 1),labels=c("**", "*",""),include.lowest = TRUE)
  nstudy<-length(unique(test$study[!is.na(test$study)]))
  my.lines<-data.frame(x=(nstudy-0.5), y=0.5, xend=(nstudy-0.5), yend=(length(unique(test$taxa))+0.5))
  h<-ggplot2::ggplot(test, ggplot2::aes(pop, plotvar)) +
    ggplot2::geom_tile(ggplot2::aes(fill=esticat)) +
    ggplot2::geom_text(ggplot2::aes(label = pdot)) +
    ggplot2::scale_fill_manual(breaks=levels(test$esticat),
                      values = levels(test$esticol),
                      labels = levels(test$esticat),
                      name = "log(OR)")+
    #ggplot2::scale_y_discrete(limits = rev(levels(test$taxas)))+
    ggplot2::geom_segment(data=my.lines, ggplot2::aes(x,y,xend=xend, yend=yend), size=2, inherit.aes=F)+
    ggplot2::ylab("") +ggplot2::xlab("")+
    ggplot2::theme(legend.title = ggplot2::element_text(size = 10),
          legend.text = ggplot2::element_text(size = leg.text.size),
          plot.title = ggplot2::element_text(size=16),
          axis.title=ggplot2::element_text(size=14,face="bold"),
          legend.position="left",
          plot.background = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.grid.major = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank(),
          panel.border = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank(),
          axis.text.x =ggplot2::element_text(size=heat.text.x.size, angle=heat.text.x.angle, hjust = 1),
          legend.key.size = ggplot2::unit(leg.key.size, "cm"))+
    ggplot2::guides(fill=ggplot2::guide_legend(ncol=1))
  #forest plot
  testf<-metadat$taxsig
  testf$taxa<-testf$id
  testf$taxa<-gsub("k__bacteria.p__","",testf$taxa)
  testf<-testf[testf$taxa %in% unique(test$taxa),]
  testf$taxa<-gdata::drop.levels(testf$taxa, reorder=FALSE)
  if (sumtype=="taxa"){
    if (level=="main"){
      testf$taxa2<-paste0(toupper(substr(as.character(testf$taxa), 1, 1)), substr(as.character(testf$taxa), 2, nchar(as.character(testf$taxa))))
      testf<-testf[order(testf$taxa,decreasing = FALSE),]
      testf$taxa2<-factor(testf$taxa2,unique(testf$taxa2))
      testf$taxa<-factor(testf$taxa,unique(testf$taxa))
      if (phyla.col=="select"){
        testf$colp<-plyr::mapvalues(testf$taxa,from=c("actinobacteria","bacteroidetes","cyanobacteria","firmicutes","fusobacteria","proteobacteria","verrucomicrobia",".thermi."),
                              to=c("#dd1c77","#31a354","#91003f","#d95f0e","#636363","#2ef0e7","#862ef0","#000"))
        testf$colp<-as.character(testf$colp)
      }
      if (phyla.col=="rainbow"){
        testf$colp<-plyr::mapvalues(testf$taxa,from=levels(testf$taxa),to=rainbow(nlevels(testf$taxa)))
        testf$colp<-as.character(testf$colp)
      }
      testf$plotvar<-testf$taxa2
    }
    if (level=="sub"){
      testf$taxas1<-sub(".c__.*f__", " ",as.character(testf$taxa))
      testf$taxas2<-sub(".*f__", "",as.character(testf$taxa))
      testf$taxas2<-paste0(toupper(substr(as.character(testf$taxas2), 1, 1)), substr(as.character(testf$taxas2), 2, nchar(as.character(testf$taxas2))))
      #replace empty truncated names by original names
      testf$taxas2[testf$taxas2 %in% c("",".g__")]<-testf$taxa[testf$taxas2 %in% c("",".g__")]
      testf<-testf[order(testf$taxa,decreasing = FALSE),]
      testf$taxas<-factor(testf$taxas2,levels=testf$taxas2)
      testf$phylum<-sub(".c__.*", "",as.character(testf$taxa))
      testf$phylum<-as.factor(testf$phylum)
      if (phyla.col=="select"){
        testf$colp<-plyr::mapvalues(testf$phylum,from=c("actinobacteria","bacteroidetes","cyanobacteria","firmicutes","fusobacteria","proteobacteria","verrucomicrobia",".thermi."),
                              to=c("#dd1c77","#31a354","#91003f","#d95f0e","#636363","#2ef0e7","#862ef0","#000"))
        testf$colp<-as.character(testf$colp)
      }
      if (phyla.col=="rainbow"){
        testf$colp<-plyr::mapvalues(testf$phylum,from=levels(testf$phylum),to=rainbow(nlevels(testf$phylum)))
        testf$colp<-as.character(testf$colp)
      }
      testf$plotvar<-testf$taxas
    }
  }
  if (sumtype=="path"){
    testf<-testf[order(testf$taxa,decreasing = FALSE),]
    testf$taxas<-factor(testf$taxa,levels=testf$taxa)
    testf$colp=1
    testf$plotvar<-testf$taxas
  }
  testf$pcut<-testf[,p]
  testf$p.adjustcut<-testf[,p.adjust]
  testf$psig<-cut(testf$pcut, breaks=c(0,0.05,1),include.lowest = TRUE, right = FALSE)
  testf$psigcol<-plyr::mapvalues(testf$psig,from=c("[0,0.05)","[0.05,1]"),to=c("red", "black"))
  testf$psigcol<-gdata::drop.levels(testf$psigcol,reorder=FALSE)
  testf$padjustsig<-cut(testf$p.adjustcut, breaks=c(0,0.1,1),include.lowest = TRUE, right = FALSE)
  testf$estimate<-as.numeric(as.character(testf$estimate))
  testf$padjustsign<-plyr::mapvalues(testf$padjustsig,from=c("[0,0.1)","[0.1,1]"),to=c("17","16"))
  testf$padjustsize<-plyr::mapvalues(testf$padjustsig,from=c("[0,0.1)","[0.1,1]"),to=c("2","1"))
  # dirty truncate large estimate, LL and UL for better plot view
  testf[,c("estimate","ll","ul")]<-apply(testf[,c("estimate","ll","ul")],2,function(x){x[x>=5]=5;x[x<=-5]=-5;x})
  f<-ggplot2::ggplot(data=testf,ggplot2::aes(x=estimate,y=plotvar,colour=psigcol))+
    ggplot2::geom_point(shape=as.numeric(as.character(testf$padjustsign)),size=as.numeric(as.character(testf$padjustsize)))+
    ggplot2::scale_y_discrete(position = "right")+ #,limits = rev(levels(testf$taxa2))
    ggplot2::geom_errorbarh(ggplot2::aes(xmin=ll,xmax=ul,colour=psigcol),height=0.0)+
    ggplot2::geom_vline(xintercept=0,linetype="dashed")+
    ggplot2::scale_colour_manual(breaks=testf$psigcol,values = levels(testf$psigcol))+
    ggplot2::theme(legend.position="none",
          plot.background = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank(),
          axis.ticks.y= ggplot2::element_blank(),
          axis.title = ggplot2::element_blank(),
          axis.text.y =ggplot2::element_text(size=forest.axis.text.y, colour = testf$colp),
          axis.text.x =ggplot2::element_text(size=forest.axis.text.x))
  return(gridExtra::grid.arrange(h,f,nrow=1))
}
