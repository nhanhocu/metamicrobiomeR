#' Display abundance comparison results.
#'
#' This function displays taxa/pathway abundance comparison results as table.
#' @param taxcomtab table of taxa abundance comparison generated from taxa.compare.
#' @param sumvar Options are "taxa" for bacterial taxa and "path" for pathway. Default is "taxa"
#' @param tax.lev taxa level to be displayed. Options are from "l2" (phylum) to "l7" (species). Default is "l2".
#' @param tax.select selected list of taxa to be displayed. Default is "none" or display all available taxa.
#' @param showvar variable (pattern) in the model to be displayed.
#' @param readjust.p multiple testing re-adjustment for only the level to be displayed (TRUE) or keep original multiple testing adjustment for all taxa of all levels (FALSE).
#' @param p.adjust.method method for multiple testing adjustment. Available options are those of the p.adjust function. Default is "fdr".
#' @param p.cutoff cutoff p-value to be displayed. Default is 0.05.
#' @param digit digit for estimates and 95 CI. Default is 2.
#' @param p.digit digit for p-values. Default is 4.
#' @return a table of results.
#' @keywords abundance comparison display.
#' @export
#' @examples
#' #Load summary tables of bacterial taxa relative abundance from Bangladesh data
#' data(taxtab.rm7)
#' #Comparison of bacterial taxa relative abundance using LMEM or GAMLSS (take time to run)
#' taxacom6.rmg<-taxa.compare(taxtab=taxtab6.rm[[5]],taxsum="rel",propmed.rel="lm",comvar="bf",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
#' taxcomtab.show(taxcomtab=taxacom6.rmg,tax.select="none", showvar="bfNon_exclusiveBF", tax.lev="l2",readjust.p=TRUE,p.adjust.method="fdr")


taxcomtab.show<-function(taxcomtab,sumvar="taxa", tax.lev="l2",tax.select="none",showvar,readjust.p=FALSE,p.adjust.method="fdr",p.cutoff=0.05,digit=2,p.digit=4,...){
  mtaba1<-taxcomtab[,c("id",colnames(taxcomtab)[grep(showvar,colnames(taxcomtab))])]
  if (sumvar=="taxa"){
    al2<-mtaba1$id[-grep("c__",mtaba1$id)]
    al3<-mtaba1$id[-grep("o__",mtaba1$id)]
    al4<-mtaba1$id[-grep("f__",mtaba1$id)]
    al5<-mtaba1$id[-grep("g__",mtaba1$id)]
    mtabal<-list(l2=mtaba1[mtaba1$id %in% al2,],l3=mtaba1[mtaba1$id %in% al3,],l4=mtaba1[mtaba1$id %in% al4,],l5=mtaba1[mtaba1$id %in% al5,],l6=mtaba1)
    mtaba<-mtabal[[tax.lev]]
    if (tax.lev=="l7"){
      taxuse<-mtaba$id[grep(".s__",mtaba$id)]
    }
    if (tax.lev=="l6"){
      taxuse<-mtaba$id[grep(".g__",mtaba$id)]
    }
    if (tax.lev=="l5"){
      taxuse<-mtaba$id[grep(".f__",mtaba$id)]
    }
    if (tax.lev=="l4"){
      taxuse<-mtaba$id[grep(".o__",mtaba$id)]
    }
    if (tax.lev=="l3"){
      taxuse<-mtaba$id[grep(".c__",mtaba$id)]
    }
    if (tax.lev=="l2"){
      taxuse<-mtaba$id
    }
  }
  if (sumvar=="path"){
    mtaba<-mtaba1
    taxuse<-mtaba$id
  }
  if (tax.select!="none"){
    taxuse<-taxuse[taxuse %in% tax.select]
  }
  datuse<-as.data.frame(mtaba[mtaba$id %in% taxuse,])
  if (readjust.p==TRUE){
    datuse[,sub('.*\\.', 'pval.adjust.',colnames(datuse)[grep("Pr(>|t|)",colnames(datuse))])]<-p.adjust(datuse[,colnames(datuse)[grep("Pr(>|t|)",colnames(datuse))]],method = p.adjust.method)
  }
  datuse<-as.data.frame(datuse)
  datsig<-datuse[datuse[,colnames(datuse)[grep("Pr(>|t|)",colnames(datuse))]]<p.cutoff,]
  datsig<-datsig[order(datsig[,colnames(datsig)[grep("Pr(>|t|)",colnames(datsig))]]),]
  datsig[,"ll"]<-datsig[,colnames(datsig)[grep("Estimate",colnames(datsig))]]-1.96*datsig[,colnames(datsig)[grep("Std. Error",colnames(datsig))]]
  datsig[,"ul"]<-datsig[,colnames(datsig)[grep("Estimate",colnames(datsig))]]+1.96*datsig[,colnames(datsig)[grep("Std. Error",colnames(datsig))]]
  datsig[,c(colnames(datsig)[grep("Estimate",colnames(datsig))],"ll","ul")]<-round(datsig[,c(colnames(datsig)[grep("Estimate",colnames(datsig))],"ll","ul")],digit)
  datsig[,c(colnames(datsig)[grep("Pr(>|t|)",colnames(datsig))],colnames(datsig)[grep("pval.adjust.",colnames(datsig))])]<-round(datsig[,c(colnames(datsig)[grep("Pr(>|t|)",colnames(datsig))],colnames(datsig)[grep("pval.adjust.",colnames(datsig))])],p.digit)
  datsig<-datsig[,c("id",colnames(datsig)[grep("Estimate",colnames(datsig))],"ll","ul",colnames(datsig)[grep("Pr(>|t|)",colnames(datsig))],colnames(datsig)[grep("pval.adjust.",colnames(datsig))])]
  return(datsig)
}
