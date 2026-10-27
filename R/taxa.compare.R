#' Compare taxa relative abundance
#'
#' This function compares taxa relative abundance summary tables at all levels between groups using GAMLSS with BEZI or Linear/Linear Mixed Effect models (LM/LMEM) after filtering (using prevalence and relative abundance thresholds).
#' @param taxtab taxa relative abundance table (already merged to mapping file) from phylum to species or any preferred highest taxa level.
#' @param propmed.rel statistical method for comparing relative abundance. Options are "lm" for LM/LMEM or "gamlss" for GAMLSS with BEZI family.
#' @param transform transformation of relative abundance data. Options are "none" for no transformation, "gmpr" for Geometric Mean of Pairwise Ratios (GMPR) normalization, "asin.sqrt" for arcsine transformation, "logit" for logit transformation, "clr" for centered log ratio transformation. Default is "none".
#' @param zeroreplace.method Method for zero replacement implemented in R package *zCompositions*. Options are "none" for no replacement, "multKM" for Multiplicative Kaplan-Meier smoothing spline (KMSS) replacement, "multLN" for Multiplicative lognormal replacement, "multRepl" for Multiplicative simple replacement, "lrEM" for Log-ratio EM algorithm, "lrDA" for Log-ratio DA algorithm. Default is "none".
#' @param comvar main variable for comparison
#' @param adjustvar variables to be adjusted.
#' @param personid name of variable for person id (applicable for longitudinal data)
#' @param longitudinal whether data is longitudinal? Options are "yes" or "no". Default is "yes".
#' @param p.adjust.method method for multiple testing adjustment. Options are those of the p.adjust.methods of stats:: p.adjust function. Default for this function is "fdr".
#' @param percent.filter prevalence threshold (the percentage of number of samples the taxa/pathway available). Default is 0.05.
#' @param relabund.filter relative abundance threshold (the minimum of the average relative abundance for a taxa/pathway to be retained). Default is 0.00005.
#' @return matrice of coefficients, standard errors, p-values and multiple testing adjusted p-values of all variables in the models.
#' @keywords taxa abundance comparison
#' @import gamlss
#' @importFrom matrixStats colMedians
#' @rawNamespace import(lme4, except = refit)
#' @export
#' @examples
#' #Load summary tables of bacterial taxa relative abundance from Bangladesh data
#' data(taxtab6)
#' tab6<-as.data.frame(taxtab6)
#' tl<-colnames(taxtab6)[grep("k__bacteria.p__fusobacteria",colnames(taxtab6))]
#' taxacom.ex<-taxa.compare(taxtab=tab6[,c("personid","x.sampleid","bf","age.sample",tl)],
#' propmed.rel="gamlss",comvar="bf",adjustvar="age.sample",
#' longitudinal="yes",p.adjust.method="fdr")


taxa.compare<-function(taxtab,propmed.rel="gamlss",transform="none",zeroreplace.method="none",
                       comvar,adjustvar,personid="personid",longitudinal="yes",
                       percent.filter=0.05,relabund.filter=0.00005,p.adjust.method="fdr"){
  taxdat<-as.data.frame(taxtab)
  taxdat[,comvar]<-gdata::drop.levels(taxdat[,comvar],reorder=FALSE) #drop missing/unused level and keep level order
  if (longitudinal=="yes"){
    taxdat$personid<-as.factor(taxdat[,personid])
  }
  # get assigned taxa only
  taxlist<-colnames(taxdat)[grep("k__",colnames(taxdat))]
  #filter using percent.filter
  taxtest<-apply(taxdat[,taxlist],2,function(x){length(x[!is.na(x)&x>0])})
  taxget<-taxtest[taxtest>=percent.filter*(nrow(taxdat))]
  #filter using relabund.filter
  taxtestm<-apply(taxdat[,taxlist],2,mean,na.rm=T)
  taxgetm<-taxtestm[taxtestm>relabund.filter]
  taxname<-names(taxget)[names(taxget) %in% names(taxgetm)]
  #transformation of relative abundance
  if (propmed.rel=="gamlss" &transform!="none"){
    stop("gamlss with beta zero-inflated family should only be used for relative abundance without transformation")
  }
  if (transform!="clr" &zeroreplace.method!="none"){
    stop("Zero replacement is only implemented for use with CLR transformation")
  }
  if (transform=="clr" &zeroreplace.method=="none"){
    stop("Zero replacement needs to be done before CLR transformation")
  }
  if (propmed.rel=="lm" &transform=="asin.sqrt"){
    asintransform <- function(p) { asin(sqrt(p)) }
    taxdat[,taxname]<-apply(taxdat[,taxname],2,asintransform)
  }
  if (propmed.rel=="lm" &transform=="logit"){
    logittransform <- function(p) { log(p/(1-p)) }
    taxdat[,taxname]<-apply(taxdat[,taxname],2,logittransform )
  }
  # Centered log ratio transformation
  if (propmed.rel=="lm" &transform=="clr"){
    #zero replacement using package zCompositions
    if (zeroreplace.method=="multLN"){
      test0<-zCompositions::multLN(taxdat[,taxname],label=0,dl=rep(1,length(taxname)))
    }
    if (zeroreplace.method=="multKM"){
      test0<-zCompositions::multKM(taxdat[,taxname],label=0,dl=rep(1,length(taxname)))
    }
    if (zeroreplace.method=="multRepl"){
      test0<-zCompositions::multRepl(taxdat[,taxname],label=0,dl=rep(1,length(taxname)))
    }
    if (zeroreplace.method=="lrEM"){
      test0<-zCompositions::lrEM(taxdat[,taxname],label=0,dl=rep(1,length(taxname)))
    }
    if (zeroreplace.method=="lrDA"){
      test0<-zCompositions::lrDA(taxdat[,taxname],label=0,dl=rep(1,length(taxname)))
    }
    #CLR transformation of imputed data using package compositions
    clrdat<-as.data.frame(compositions::clr(test0))
    taxdat[,taxname]<-clrdat
  }

  #GMPR normalization
  #GMPR script is copied from github repo: "https://github.com/jchen1981/GMPR/blob/master/GMPR.R"
  GMPR <- function (comm, intersect.no = 10, ct.min = 1, trace = TRUE) {
    # mask counts < ct.min
    comm[comm < ct.min] <- 0

    if (is.null(colnames(comm))) {
      colnames(comm) <- paste0('S', 1:ncol(comm))
    }

    if (trace) message('Begin GMPR size factor calculation ...\n')

    comm.no <- numeric(ncol(comm))
    gmpr <- sapply(1:ncol(comm),  function(i) {
      if (i %% 50 == 0) {
        message(i, '\n')
      }
      x <- comm[, i]
      # Compute the pairwise ratio
      pr <- x / comm
      # Handling of the NA, NaN, Inf
      pr[is.nan(pr) | !is.finite(pr) | pr == 0] <- NA
      # Counting the number of non-NA, NaN, Inf
      incl.no <- colSums(!is.na(pr))
      # Calculate the median of PR
      pr.median <- matrixStats::colMedians(pr, na.rm=TRUE)
      # Record the number of samples used for calculating the GMPR
      #comm.no[i] <<- sum(incl.no >= intersect.no)
      # Geometric mean of PR median
      if (comm.no[i] > 1) {
        return(exp(mean(log(pr.median[incl.no >= intersect.no]))))
      } else {
        return(NA)
      }
    }
    )
    if (sum(is.na(gmpr))) {
      warning(paste0('The following samples\n ', paste(colnames(comm)[is.na(gmpr)], collapse='\n'),
                     '\ndo not share at least ', intersect.no, ' common taxa with the rest samples! ',
                     'For these samples, their size factors are set to be NA! \n',
                     'You may consider removing these samples since they are potentially outliers or negative controls!\n',
                     'You may also consider decreasing the minimum number of intersecting taxa and rerun the procedure!\n'))
    }
    if (trace) message('Completed!\n')
    if (trace) message('Please watch for the samples with limited sharing with other samples based on NSS! They may be outliers! \n')
    names(gmpr) <- names(comm.no) <- colnames(comm)
    attr(gmpr, 'NSS') <- comm.no
    return(gmpr)
  }

  #Begin GMPR normalization
  if (transform=="gmpr"){
    reltab<-t(taxdat[,taxname])
    #calculate size factor
    gmpr.sf<-GMPR(comm=reltab,intersect.no = 2, ct.min = 0, trace = TRUE)
    #normalized by size factors
    rel.gmpr<-t(t(reltab)/gmpr.sf)
    taxdat[,taxname]<-rel.gmpr
  }

  # begin models for differential abundance testing
  estisum<-NULL
  for (i in 1: length(taxname)){
    #print(i)
    #linear/linear mixed effect model
    if (propmed.rel=="lm"){
      if (longitudinal=="yes"){
        fitsum<-try(summary(lme4::lmer(stats::as.formula(paste(taxname[i],paste(c(comvar,adjustvar,"(1|personid)"),collapse="+"),sep="~")), data=taxdat)))
        #fitsum<-try(summary(lme4::glmer(stats::as.formula(paste(taxname[i],paste(c(comvar,adjustvar,"(1|personid)"),collapse="+"),sep="~")), data=taxdat,family="gaussian")))
      }
      if (longitudinal=="no"){
        fitsum<-try(summary(stats::glm(stats::as.formula(paste(taxname[i],paste(c(comvar,adjustvar),collapse="+"),sep="~")), data=taxdat,family="gaussian")))
      }
      if (class(fitsum) == "try-error") {
        warning("Error in model fit, NA introduced.\n")
        fitcoefw<-NULL
        estisum<-plyr::rbind.fill(estisum,fitcoefw)
      }
      if (class(fitsum) != "try-error") {
        if (length(which(rownames(fitsum$coefficients)!="(Intercept)"))>1){
          fitcoef<-as.data.frame(fitsum$coefficients[rownames(fitsum$coefficients)!="(Intercept)",]) #remove intercept
          if (longitudinal=="yes"){
            fitcoef[,"Pr(>|t|)"]<-1.96*stats::pnorm(-abs(fitcoef[,"Estimate"]/fitcoef[,"Std. Error"]))
          }
          fitcoef[,"varname"]<-rownames(fitcoef)
          fitcoef[,"id"]<-taxname[i]
          fitcoefw<-stats::reshape(fitcoef, idvar="id", timevar="varname", direction="wide")
        }
        #handling issue when there is one row
        if (length(which(rownames(fitsum$coefficients)!="(Intercept)"))==1){
          fitcoef<-as.data.frame(matrix(fitsum$coefficients[rownames(fitsum$coefficients)!="(Intercept)",],ncol=ncol(fitsum$coefficients)))
          rownames(fitcoef)<-rownames(fitsum$coefficients)[rownames(fitsum$coefficients)!="(Intercept)"]
          colnames(fitcoef)<-colnames(fitsum$coefficients)
          if (longitudinal=="yes"){
            fitcoef[,"Pr(>|t|)"]<-1.96*stats::pnorm(-abs(fitcoef[,"Estimate"]/fitcoef[,"Std. Error"]))
          }
          fitcoef[,"varname"]<-rownames(fitcoef)
          fitcoef[,"id"]<-taxname[i]
          fitcoefw<-stats::reshape(fitcoef, idvar="id", timevar="varname", direction="wide")
        }
        # when there is no coef
        if (length(which(rownames(fitsum$coefficients)!="(Intercept)"))==0){
          fitcoefw<-NULL
        }

        estisum<-plyr::rbind.fill(estisum,fitcoefw)
      }
    }
    #Generalized Additive Models for Location Scale and Shape (GAMLSS): Betazeroinflated (BEZI) family, mu link logit
    if (propmed.rel=="gamlss"){
      if (longitudinal == "yes") {
        testdat <- taxdat[, c(taxname[i], comvar, adjustvar, "personid")]
        testdat[, taxname[i]][testdat[, taxname[i]] >= 1] <- 0.9999
        testdat <- stats::na.omit(testdat)
        fitsum <- try(summary(gamlss::gamlss(stats::as.formula(paste(taxname[i],
                                  paste(c(comvar, adjustvar, "random(personid)"),
                                  collapse = "+"), sep = "~")), family = BEZI,
                                  data = testdat, trace = FALSE), save = TRUE))
      }
      if (longitudinal == "no") {
        testdat <- taxdat[, c(taxname[i], comvar, adjustvar)]
        testdat[, taxname[i]][testdat[, taxname[i]] >= 1] <- 0.9999
        testdat <- stats::na.omit(testdat)
        fitsum <- try(summary(gamlss::gamlss(stats::as.formula(paste(taxname[i],
                                    paste(c(comvar, adjustvar), collapse = "+"),
                                    sep = "~")), family = BEZI, data = testdat,
                                    trace = FALSE), save = TRUE))
      }
      if (class(fitsum) == "try-error") {
        warning("Error in model fit, NA introduced.\n")
        fitcoefw<-NULL
        estisum<-plyr::rbind.fill(estisum,fitcoefw)
      }
      if (class(fitsum) != "try-error") {
        if (length(which(rownames(fitsum$coef.table)!="(Intercept)"))>1){
          fitcoef<-as.data.frame(fitsum$coef.table[rownames(fitsum$coef.table)!="(Intercept)",])
          fitcoef[,"varname"]<-rownames(fitcoef)
          fitcoef[,"id"]<-taxname[i]
          fitcoefw<-stats::reshape(fitcoef, idvar="id", timevar="varname", direction="wide")
        }
        #handle the issues when there is only one row
        if (length(which(rownames(fitsum$coef.table)!="(Intercept)"))==1){
          fitcoef<-as.data.frame(matrix(fitsum$coef.table[rownames(fitsum$coef.table)!="(Intercept)",],ncol=ncol(fitsum$coef.table)))
          rownames(fitcoef)<-rownames(fitsum$coef.table)[rownames(fitsum$coef.table)!="(Intercept)"]
          colnames(fitcoef)<-colnames(fitsum$coef.table)
          fitcoef[,"varname"]<-rownames(fitcoef)
          fitcoef[,"id"]<-taxname[i]
          fitcoefw<-stats::reshape(fitcoef, idvar="id", timevar="varname", direction="wide")
        }
        # when there is no coef
        if (length(which(rownames(fitsum$coef.table)!="(Intercept)"))==0){
          fitcoefw<-NULL
        }
        estisum<-plyr::rbind.fill(estisum,fitcoefw)
      }
    }
  }
  np<-length(colnames(estisum)[grep("Pr(>|t|)",colnames(estisum))])
  if(np>1){
    estisum[,sub('.*\\.', 'pval.adjust.',colnames(estisum)[grep("Pr(>|t|)",colnames(estisum))])]<-apply(estisum[,colnames(estisum)[grep("Pr(>|t|)",colnames(estisum))]],2,stats::p.adjust,method = p.adjust.method)
  }else{
    estisum[,sub('.*\\.', 'pval.adjust.',colnames(estisum)[grep("Pr(>|t|)",colnames(estisum))])]<-stats::p.adjust(estisum[,colnames(estisum)[grep("Pr(>|t|)",colnames(estisum))]],method = p.adjust.method)
  }
  return(estisum)
}

