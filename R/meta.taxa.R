#' Meta-analysis of taxa/pathway abundance comparison.
#'
#' This function does meta-analysis based on estimates and standard errors from taxa/pathway abundance comparison using random effect and fixed effect meta-analysis models.
#' @param taxcomdat matrice of estimates and SE of all taxa/pathways combined from all included studies.
#' @param estimate.pattern string pattern for estimates. Default is "Estimate.".
#' @param se.pattern string pattern for standard error. Default is "Std. Error.".
#' @param summary.measure "RR" for estimates from GAMLSS with BEZI family and "RD" for estimates from Linear/linear mixed effect model. Default is "RR"
#' @param pool.var name of id variable for meta-analysis. Default is "id".
#' @param studylab name of variable characterizing included studies. Default is "study".
#' @param backtransform whether or not to perform backtransformation of the estimates. Default is FALSE.
#' @param percent.meta the threshold percentage of number of studies that a taxa is available to do meta-analysis. Default is 0.5
#' @param p.adjust.method method for multiple testing adjustment (available methods of the function p.adjust). Default is "fdr".
#' @return a list of matrices of results for all variables in the comparison models.
#' @keywords abundance meta-analysis.
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
#' metatab.show(metatab=metab.sex$random,com.pooled.tab=tabsex4,tax.lev="l2",showvar="genderMale",p.cutoff.type="p", p.cutoff=0.05,display="table")


meta.taxa<-function(taxcomdat,estimate.pattern="Estimate.",se.pattern="Std. Error.", summary.measure="RR",pool.var="id",studylab="study",backtransform=FALSE,percent.meta=0.5,p.adjust.method="fdr",...){
  #require(meta)
  taxcomdat$studylab<-taxcomdat[,studylab]
  es<-colnames(taxcomdat)[grep(estimate.pattern,colnames(taxcomdat))]
  ses<-colnames(taxcomdat)[grep(se.pattern,colnames(taxcomdat))]
  taxname<-unique(taxcomdat[,pool.var])
  meta.f<-list()
  meta.r<-list()
  for (i in 1:length(es)){
    print(i)
    metatab.f<-NULL
    metatab.r<-NULL
    taxcomdat$estimate<-taxcomdat[,es[i]]
    taxcomdat$se<-taxcomdat[,ses[i]]
    for (j in 1:length(taxname)){
      print(j)
      testdat<-taxcomdat[taxcomdat[,pool.var] %in% taxname[j],]
      # do metanalysis only for taxa exist in >=percent.meta of studies
      if (nrow(testdat)<percent.meta*(length(unique(taxcomdat$studylab)))){
        fitsum.f<-NULL
        fitsum.r<-NULL
        metatab.f<- plyr::rbind.fill(metatab.f,fitsum.f)
        metatab.r<- plyr::rbind.fill(metatab.r,fitsum.r)
      }
      if (nrow(testdat)>=percent.meta*(length(unique(taxcomdat$studylab)))){
        fit.meta<-meta::metagen(estimate, se, studlab=studylab,data=testdat,sm=summary.measure, backtransf=backtransform)
        #fixed effect
        fitsum.f<-as.data.frame(matrix(c(summary(fit.meta)$fixed$TE,summary(fit.meta)$fixed$seTE,summary(fit.meta)$fixed$lower,summary(fit.meta)$fixed$upper,summary(fit.meta)$fixed$z,summary(fit.meta)$fixed$p),nrow=1))
        colnames(fitsum.f)<-c("estimate",'se','ll','ul','z','p')
        fitsum.f[,pool.var]<-taxname[j]
        metatab.f<-plyr::rbind.fill(metatab.f,fitsum.f)
        #random effect
        fitsum.r<-as.data.frame(matrix(c(summary(fit.meta)$random$TE,summary(fit.meta)$random$seTE,summary(fit.meta)$random$lower,summary(fit.meta)$random$upper,summary(fit.meta)$random$z,summary(fit.meta)$random$p),nrow=1))
        colnames(fitsum.r)<-c("estimate",'se','ll','ul','z','p')
        fitsum.r[,pool.var]<-taxname[j]
        metatab.r<-plyr::rbind.fill(metatab.r,fitsum.r)
      }
    }
    metatab.f[,"p.adjust"]<-p.adjust(metatab.f[,"p"],method=p.adjust.method)
    metatab.f<-metatab.f[order(metatab.f[,"p"]),]
    meta.f[[i]]<-metatab.f
    metatab.r[,"p.adjust"]<-p.adjust(metatab.r[,"p"],method=p.adjust.method)
    metatab.r<-metatab.r[order(metatab.r[,"p"]),]
    meta.r[[i]]<-metatab.r
  }
  names(meta.f)<-names(meta.r)<-gsub(estimate.pattern,"",es)
  meta.fr<-list(meta.f,meta.r)
  names(meta.fr)<-c("fixed","random")
  return(meta.fr)
}
