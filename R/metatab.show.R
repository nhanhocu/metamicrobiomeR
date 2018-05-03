#' Display meta-analysis results.
#'
#' This function displays meta-analysis results of relative abundance as heatmap, forest plot, table or data.
#' @param metatab matrice of taxa/pathway abundance comparison meta-analysis results generated from meta.taxa.
#' @param com.pooled.tab matrice of taxa/pathway abundance comparison generated from taxa.compare or pathway.compare combined from all included studies.
#' @param sumvar Either "taxa" for taxa and "path" for pathway.
#' @param tax.lev taxa level to be displayed. Options are from "l2" (phylum) to "l7" (species). Default is "l2".
#' @param showvar variable (string pattern) in the model to be displayed.
#' @param estimate.pattern string pattern for estimates. Default is "Estimate.".
#' @param se.pattern string pattern for standard error. Default is "Std. Error.".
#' @param readjust.p multiple testing re-adjustment for only the level to be displayed (TRUE) or keep original multiple testing adjustment for all taxa of all levels (FALSE).Default is FALSE.
#' @param p.cutoff.type type of p-value for cutoff. Options are "p" for p-value or "p.adjust" for multiple testing adjusted p-value. Default is "p".
#' @param p.cutoff cutoff p-value to be displayed. Default is 0.05.
#' @param display type of display. Options are display=c("plot","table","data")
#' @param plot type of plot. Options are plot=c("heatmap","forest").
#' @param fill.value name of legend.
#' @param grid whether multiple plots will be displayed alongside. Default is FALSE.
#' @param digit digit for estimates and 95 CI. Default is 2.
#' @param p.digit digit for p-values. Default is 4.
#' @return plot table or data.
#' @keywords abundance meta-analysis display.
#' @export
#' @examples
#' #Load saved results of four studies for the comparison of bacterial taxa relative abundance between genders adjusted for breastfeeding and infant age at sample collection
#' data(taxacom.rm.sex.adjustbfage)
#' data(taxacom.ha.sex.adjustbfage)
#' data(taxacom6.zi.usbmk.sex.adjustbfage.)
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
#' metatab.show(metatab=metab.sex$random,com.pooled.tab=tabsex4,tax.lev="l2",showvar="genderMale",p.cutoff.type="p", p.cutoff=0.2,display="table")
#' metatab.show(metatab=metab.sex$random,com.pooled.tab=tabsex4,tax.lev="l2",showvar="genderMale",p.cutoff.type="p", p.cutoff=1,display="plot",plot="heatmap")


metatab.show<-function(metatab,com.pooled.tab,sumvar="taxa",tax.lev="l2",showvar,estimate.pattern="Estimate.",se.pattern="Std. Error.",p.pattern="Pr(>|t|)",readjust.p=FALSE,p.cutoff.type="p", p.cutoff=0.05,display="plot",plot="heatmap",fill.value="log(OR)",grid=FALSE,digit=2,p.digit=4,...){
  #require(ggplot2);require(plyr)
  mtaba1<-as.data.frame(metatab[[showvar]])
  #remove row with NA values (no meta-analysis)
  mtaba1<-mtaba1[!is.na(mtaba1[,"estimate"]),]
  if (sumvar=="taxa"){
    al2<-mtaba1$id[-grep("c__",mtaba1$id)]
    al3<-mtaba1$id[-grep("o__",mtaba1$id)]
    al4<-mtaba1$id[-grep("f__",mtaba1$id)]
    al5<-mtaba1$id[-grep("g__",mtaba1$id)]
    #al6<-mtaba1$id[-grep("s__",mtaba1$id)]
    mtabal<-list(l2=mtaba1[mtaba1$id %in% al2,],l3=mtaba1[mtaba1$id %in% al3,],l4=mtaba1[mtaba1$id %in% al4,],l5=mtaba1[mtaba1$id %in% al5,],l6=mtaba1)
    mtaba<-mtabal[[tax.lev]]
    if (tax.lev=="l2"){
      mtab<-mtaba[(mtaba$id %in% mtaba$id[grep("p__",mtaba$id)]),]
      if (readjust.p==TRUE){
        mtab[,"p.adjust"]<-p.adjust(mtab[,"p"],method="fdr")
      }
      taxsig<-mtab[mtab[,p.cutoff.type]<=p.cutoff& !is.na(mtab[,p.cutoff.type]),c("estimate",'se','ll','ul','z','p','p.adjust','id')]
    }
    if (tax.lev=="l3"){
      mtab<-mtaba[(mtaba$id %in% mtaba$id[grep("c__",mtaba$id)]),]
      if (readjust.p==TRUE){
        mtab[,"p.adjust"]<-p.adjust(mtab[,"p"],method="fdr")
      }
      taxsig<-mtab[(mtab$id %in% mtab$id[grep("c__",mtab$id)]) &mtab[,p.cutoff.type]<=p.cutoff& !is.na(mtab[,p.cutoff.type]),c("estimate",'se','ll','ul','z','p','p.adjust','id')]
    }
    if (tax.lev=="l4"){
      mtab<-mtaba[(mtaba$id %in% mtaba$id[grep("o__",mtaba$id)]),]
      if (readjust.p==TRUE){
        mtab[,"p.adjust"]<-p.adjust(mtab[,"p"],method="fdr")
      }
      taxsig<-mtab[(mtab$id %in% mtab$id[grep("o__",mtab$id)]) &mtab[,p.cutoff.type]<=p.cutoff& !is.na(mtab[,p.cutoff.type]),c("estimate",'se','ll','ul','z','p','p.adjust','id')]
    }
    if (tax.lev=="l5"){
      mtab<-mtaba[(mtaba$id %in% mtaba$id[grep("f__",mtaba$id)]),]
      if (readjust.p==TRUE){
        mtab[,"p.adjust"]<-p.adjust(mtab[,"p"],method="fdr")
      }
      taxsig<-mtab[(mtab$id %in% mtab$id[grep("f__",mtab$id)]) &mtab[,p.cutoff.type]<=p.cutoff& !is.na(mtab[,p.cutoff.type]),c("estimate",'se','ll','ul','z','p','p.adjust','id')]
    }
    if (tax.lev=="l6"){
      mtab<-mtaba[(mtaba$id %in% mtaba$id[grep("g__",mtaba$id)]),]
      if (readjust.p==TRUE){
        mtab[,"p.adjust"]<-p.adjust(mtab[,"p"],method="fdr")
      }
      taxsig<-mtab[(mtab$id %in% mtab$id[grep("g__",mtab$id)]) &mtab[,p.cutoff.type]<=p.cutoff& !is.na(mtab[,p.cutoff.type]),c("estimate",'se','ll','ul','z','p','p.adjust','id')]
    }
    if (tax.lev=="l7"){
      mtab<-mtaba[(mtaba$id %in% mtaba$id[grep("s__",mtaba$id)]),]
      if (readjust.p==TRUE){
        mtab[,"p.adjust"]<-p.adjust(mtab[,"p"],method="fdr")
      }
      taxsig<-mtab[(mtab$id %in% mtab$id[grep("s__",mtab$id)]) &mtab[,p.cutoff.type]<=p.cutoff& !is.na(mtab[,p.cutoff.type]),c("estimate",'se','ll','ul','z','p','p.adjust','id')]
    }
  }
  if (sumvar=="path"){
    mtab<-mtaba1[!is.na(mtaba1$estimate),]
    if (readjust.p==TRUE){
      mtab[,"p.adjust"]<-p.adjust(mtab[,"p"],method="fdr")
    }
    taxsig<-mtab[mtab[,p.cutoff.type]<=p.cutoff& !is.na(mtab[,p.cutoff.type]),c("estimate",'se','ll','ul','z','p','p.adjust','id')]
  }
  if (nrow(taxsig)>0){
    taxsig$study<-"Meta_analysis"
    taxsig$pop<-"Pooled"
    taxsig<-taxsig[order(taxsig$id,decreasing = TRUE),]
    taxsig$id<-as.factor(taxsig$id)
    #each study
    pooltab<-com.pooled.tab[,c("id","study","pop",colnames(com.pooled.tab)[grep(showvar,colnames(com.pooled.tab))])]
    colnames(pooltab)[grep(estimate.pattern,colnames(pooltab))]<-"estimate"
    colnames(pooltab)[grep(se.pattern,colnames(pooltab))]<-"se"
    colnames(pooltab)[grep(p.pattern,colnames(pooltab))]<-"p"
    taxsig.stud<-pooltab[pooltab$id %in% taxsig$id,c("estimate","se","p","id","study","pop")]
    taxsig.all<-plyr::rbind.fill(taxsig.stud,taxsig)
    taxsig.all$id<-as.factor(taxsig.all$id)
    taxsig.all$study<-factor(taxsig.all$study, levels=c(unique(taxsig.stud$study),"Meta_analysis"))
    taxsig.all$pop<-factor(taxsig.all$pop, levels=c(unique(taxsig.stud$pop),"Pooled"))
    rownames(taxsig.all)<-NULL
  }
  if (nrow(taxsig)==0){
    display="table"
  }
  if (display=="plot"){
    if (plot=="heatmap"){
      p<-ggplot2::ggplot(taxsig.all, ggplot2::aes(study, id)) +
        ggplot2::geom_tile(ggplot2::aes(fill = estimate)) +
        ggplot2::scale_fill_gradient2(low = "blue", high = "red")+
        ggplot2::ylab("") +
        ggplot2::xlab("") +
        ggplot2::theme(legend.title = ggplot2::element_text(size = 10),
              legend.text = ggplot2::element_text(size = 8),
              plot.title = ggplot2::element_text(size=16),
              axis.title=ggplot2::element_text(size=14,face="bold"),
              axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
              legend.position="bottom") +
        ggplot2::labs(fill = fill.value)
      if (grid==TRUE){
        p<-ggplot2::ggplot(taxsig.all, ggplot2::aes(study, id)) +
          ggplot2::geom_tile(ggplot2::aes(fill = estimate)) +
          ggplot2::scale_fill_gradient2(low = "blue", high = "red")+
          ggplot2::ylab("") +
          ggplot2::xlab("") +
          ggplot2::theme(legend.title = ggplot2::element_blank(),
                legend.text = ggplot2::element_blank(),
                plot.title = ggplot2::element_text(size=16),
                axis.title=ggplot2::element_text(size=14,face="bold"),
                axis.text.x = ggplot2::element_blank(),
                axis.ticks.x = ggplot2::element_blank(),
                legend.position="none") +
          ggplot2::labs(fill = fill.value)
      }
    }
    if (plot=="forest"){
      p<-ggplot2::ggplot(data=taxsig,ggplot2::aes(x=estimate,y=id))+
        ggplot2::geom_point(shape=16)+
        ggplot2::geom_errorbarh(ggplot2::aes(xmin=ll,xmax=ul),height=0.0, colour="blue")+
        ggplot2::geom_vline(xintercept=0,linetype="dashed")+
        ggplot2::xlab(fill.value)
      if (grid==TRUE){
        p<-ggplot2::ggplot(data=taxsig,ggplot2::aes(x=estimate,y=id))+
          ggplot2::geom_point(shape=16)+
          ggplot2::geom_errorbarh(ggplot2::aes(xmin=ll,xmax=ul),height=0.0, colour="blue")+
          ggplot2::geom_vline(xintercept=0,linetype="dashed")+
          ggplot2::xlab(fill.value)+ggplot2::ylab("") +
          ggplot2::theme(axis.title.y=ggplot2::element_blank(),
                axis.text.y=ggplot2::element_blank())
      }
    }
    return(p)
  }
  if (display=="table"){
    taxsigo<-taxsig[,c('id',"estimate",'ll','ul','p','p.adjust')]
    taxsigo[,c("estimate",'ll','ul')]<-round(taxsigo[,c("estimate",'ll','ul')],digit)
    taxsigo[,c('p','p.adjust')]<-round(taxsigo[,c('p','p.adjust')],p.digit)
    taxsigo<-taxsigo[order(taxsigo[,p.cutoff.type]),]
    return(taxsigo)
  }
  if (display=="data"){
    return(data=list(taxsig.all=taxsig.all,taxsig=taxsig))
  }
}
