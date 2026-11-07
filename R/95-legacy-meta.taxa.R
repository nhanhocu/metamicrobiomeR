#' Meta-analysis of taxa/pathway abundance comparison.
#'
#' This function does meta-analysis based on estimates and standard errors from taxa/pathway abundance comparison using random effect and fixed effect meta-analysis models.
#' @param taxcomdat matrix of estimates and SE of all taxa/pathways combined from all included studies.
#' @param estimate.pattern string pattern for estimates. Default is "Estimate.".
#' @param se.pattern string pattern for standard error. Default is "Std. Error.".
#' @param summary.measure "RR" for estimates from GAMLSS with BEZI family and "RD" for estimates from Linear/linear mixed effect model. Default is "RR"
#' @param pool.var name of id variable for meta-analysis. Default is "id".
#' @param studylab name of variable characterizing included studies. Default is "study".
#' @param backtransform whether or not to perform backtransformation of the estimates. Default is FALSE.
#' @param percent.meta the threshold percentage of number of studies that a taxa is available to do meta-analysis. Default is 0.5
#' @param p.adjust.method method for multiple testing adjustment (available methods of the function p.adjust). Default is "fdr".
#' @param method.tau estimator of the between-study variance passed to meta::metagen. Default is "DL" (DerSimonian-Laird), the estimator of the 2019 paper. Set explicitly because the meta package now defaults to REML.
#' @return a list of matrices of results for all variables in the comparison models. A taxon is pooled for a term only if at least two studies, and at least `percent.meta` of the studies, have an estimate for that term. `attr(result, "dominance")` lists taxa where one study holds over 90% of the random-effects weight, and the function warns about them, because the pooled result then rests on that study.
#' @keywords abundance meta-analysis.
#' @export
#' @examples
#' # load saved GAMLSS-BEZI results of four studies
#' # for the comparison of bacterial taxa relative abundance between
#' # genders adjusted for breastfeeding and infant age at sample collection
#' data(tabsex4)
#' #select only taxonomies of a small phylum for meta-analysis example
#' # (to save running time)
#' tlm<-tabsex4$id[grep("k__bacteria.p__fusobacteria",tabsex4$id)]
#' # meta-analysis
#' metab.sex<-meta.taxa(taxcomdat=tabsex4[tabsex4$id %in% tlm,],
#' summary.measure="RR", pool.var="id", studylab="study",
#' backtransform=FALSE, percent.meta=0.5, p.adjust.method="fdr")
#' #show results by table and plot
#' #phylum
#' #table
#' metatab.show(metatab=metab.sex$random,com.pooled.tab=tabsex4[tabsex4$id %in% tlm,],
#' tax.lev="l2",showvar="genderMale",p.cutoff.type="p", p.cutoff=1,display="table")


meta.taxa<-function(taxcomdat,estimate.pattern="Estimate.",se.pattern="Std. Error.",
                    summary.measure="RR",pool.var="id",studylab="study",
                    backtransform=FALSE,percent.meta=0.5,p.adjust.method="fdr",method.tau="DL"){
  taxcomdat$studylab<-taxcomdat[,studylab]
  es<-colnames(taxcomdat)[grep(estimate.pattern,colnames(taxcomdat))]
  ses<-colnames(taxcomdat)[grep(se.pattern,colnames(taxcomdat))]
  taxname<-unique(taxcomdat[,pool.var])
  nstudy<-length(unique(taxcomdat$studylab))
  meta.f<-list()
  meta.r<-list()
  dominance<-NULL
  for (i in 1:length(es)){
    metatab.f<-NULL
    metatab.r<-NULL
    n.dropped<-0
    n.failed<-0
    term<-gsub(estimate.pattern,"",es[i])
    taxcomdat$estimate<-taxcomdat[,es[i]]
    taxcomdat$se<-taxcomdat[,ses[i]]
    for (j in 1:length(taxname)){
      testdat<-taxcomdat[taxcomdat[,pool.var] %in% taxname[j],]
      # count only the studies with an estimate for this term
      testdat<-testdat[!is.na(testdat$estimate)&!is.na(testdat$se),]
      if (nrow(testdat)<max(2,percent.meta*nstudy)){
        n.dropped<-n.dropped+1
        next
      }
      fit.meta<-tryCatch(withCallingHandlers(
        meta::metagen(estimate, se, studlab=studylab,data=testdat,sm=summary.measure,
                      backtransf=backtransform,method.tau=method.tau,method.tau.ci=""),
        warning=function(w){
          # the warning below on one study holding most of the weight says this more usefully
          if (grepl("sampling variance extremely large",conditionMessage(w))) invokeRestart("muffleWarning")
        }),error=function(e) NULL)
      if (is.null(fit.meta)){
        n.failed<-n.failed+1
        next
      }
      #fixed effect
      fitsum.f<-data.frame(estimate=fit.meta$TE.common,se=fit.meta$seTE.common,ll=fit.meta$lower.common,
                           ul=fit.meta$upper.common,statistic=fit.meta$statistic.common,p=fit.meta$pval.common)
      fitsum.f[,pool.var]<-taxname[j]
      metatab.f<-plyr::rbind.fill(metatab.f,fitsum.f)
      #random effect
      fitsum.r<-data.frame(estimate=fit.meta$TE.random,se=fit.meta$seTE.random,ll=fit.meta$lower.random,
                           ul=fit.meta$upper.random,statistic=fit.meta$statistic.random,p=fit.meta$pval.random)
      fitsum.r[,pool.var]<-taxname[j]
      metatab.r<-plyr::rbind.fill(metatab.r,fitsum.r)
      #one study dominating the random-effects result
      wt<-fit.meta$w.random/sum(fit.meta$w.random)
      if (max(wt,na.rm=TRUE)>0.9){
        dominance<-rbind(dominance,data.frame(term=term,id=as.character(taxname[j]),
                                              top_study=as.character(testdat$studylab[which.max(wt)]),
                                              max_weight=max(wt,na.rm=TRUE),stringsAsFactors=FALSE))
      }
    }
    if (n.dropped>0){
      warning(term,": ",n.dropped," of ",length(taxname)," taxa were not pooled, because too few studies have an estimate for this term.",call.=FALSE)
    }
    if (n.failed>0){
      warning(term,": ",n.failed," taxa could not be pooled, because the meta-analysis model failed.",call.=FALSE)
    }
    if (is.null(metatab.r)){
      meta.f[[i]]<-metatab.f
      meta.r[[i]]<-metatab.r
      next
    }
    metatab.f[,"p.adjust"]<-stats::p.adjust(metatab.f[,"p"],method=p.adjust.method)
    metatab.f<-metatab.f[order(metatab.f[,"p"]),]
    meta.f[[i]]<-metatab.f
    metatab.r[,"p.adjust"]<-stats::p.adjust(metatab.r[,"p"],method=p.adjust.method)
    metatab.r<-metatab.r[order(metatab.r[,"p"]),]
    meta.r[[i]]<-metatab.r
  }
  names(meta.f)<-names(meta.r)<-gsub(estimate.pattern,"",es)
  meta.fr<-list(meta.f,meta.r)
  names(meta.fr)<-c("fixed","random")
  if (!is.null(dominance)){
    attr(meta.fr,"dominance")<-dominance
    warning(nrow(dominance)," pooled results have one study holding over 90% of the random-effects weight, so they rest on that study. See attr(result, \"dominance\").",call.=FALSE)
  }
  return(meta.fr)
}
