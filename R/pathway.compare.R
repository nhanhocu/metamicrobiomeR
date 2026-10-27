#' Compare (kegg) pathway abundance
#'
#' This is a slightly modified version of the taxa.compare function.
#' It compares pathway abundance generated from PICRUSt analysis between groups using different methods (apply to pathway summary tables already merged to mapping file and put in a list (e.g.level 1, 2 and 3)).
#' Specifically, it compares relative abundances of bacterial functional pathways at all levels using GAMLSS or LM/LMEM and compares of log(absolute abundances) of bacterial functional pathways at all levels using LM/LMEM.
#' @param pathtab list of pathway abundance table of all levels.
#' @param mapfile mapping file or file containing covariates.
#' @param sampleid variable containing sample id to be matched between pathway abundance table and mapping file.
#' @param pathsum type of abundance to be compared. Options are "rel" for relative abundance or "log" for log of absolute abundance. Default is "rel".
#' @param stat.med statistical method for comparison. stat.med can be "lm" for LM/LMEM (usable for both pathsum="rel" or "log") or "gamlss" for GAMLSS with BEZI family (gamlss only make sense if pathsum="rel" ).
#' @param transform transformation of relative abundance data. Options are "none" for no transformation, "asin.sqrt" for arcsine transformation, "logit" for logit transformation. Default is "none".
#' @param comvar main variable for comparison.
#' @param adjustvar variables to be adjusted.
#' @param personid name of variable for person id in mapping file (applicable for longitudinal data)
#' @param longitudinal whether the data is longitudinal. Default is "yes".
#' @param p.adjust.method method for multiple testing adjustment. Available options are those of the p.adjust function. Default is "fdr".
#' @param percent.filter prevalence threshold (the percentage of number of samples the taxa/pathway available). Default is 0.05.
#' @param relabund.filter relative abundance threshold (the minimum of the average relative abundance for a taxa/pathway to be retained). Default is 0.00005.
#' @param pooldata whether the data is pooled from multiple studies. Default is FALSE.
#' @return matrice of coefficients, standard errors, p-values and multiple testing adjusted p-values of all variables in the models.
#' @keywords pathway abundance comparison.
#' @import gamlss
#' @rawNamespace import(lme4, except = refit)
#' @export
#' @examples
#' #Load Bangladesh pathway and metadata
#' data(kegg.12)
#' data(covar.rm)
#' # Comparison of pathway relative abundances for some first pathways of level 1 only
#' # and assuming crosssectional data (to save running time)
#' path1<-pathway.compare(pathtab=list(kegg.12[[1]][, 1:2]),
#' mapfile=covar.rm,sampleid="sampleid",pathsum="rel", stat.med="gamlss",
#' comvar="gender",adjustvar=c("age.sample","bf"), longitudinal="no",
#' p.adjust.method="fdr", percent.filter=0.05,relabund.filter=0.00005)
#' taxcomtab.show(taxcomtab=path1$l1, sumvar="path",tax.lev="l2",
#' tax.select="none", showvar="genderMale", p.adjust.method="fdr",p.cutoff=1)


pathway.compare<-function(pathtab,mapfile,sampleid="sampleid",pathsum="rel",
                          stat.med="gamlss",transform="none", comvar,adjustvar,personid="personid",
                          longitudinal="yes",p.adjust.method="fdr",
                          percent.filter=0.05,relabund.filter=0.00005,pooldata=FALSE){
  pathlev<-paste("l",1:length(pathtab),sep="")
  estilist<-list()
  for (j in 1:length(pathlev)){
    print(j)
    samlist<-rownames(pathtab[[j]])
    pathtab[[j]]<-as.data.frame(lapply(pathtab[[j]],as.character))
    pathtab[[j]]<-as.data.frame(lapply(pathtab[[j]],as.numeric)) #if original data are interger (! be careful if the original data are factor => can be very wrong)
    rownames(pathtab[[j]])<-samlist
    pathlist<-colnames(pathtab[[j]])
    #filter using percent.filter
    pathtest<-apply(pathtab[[j]],2,function(x){length(x[!is.na(x)&x>0])})
    pathget<-pathtest[pathtest>=percent.filter*(nrow(pathtab[[j]]))]
    #filter using relabund.filter
    pathtests<-apply(pathtab[[j]],2,function(x){sum(x,na.rm=T)/sum(pathtab[[j]])})
    pathgets<-pathtests[pathtests>relabund.filter]
    pathname<-names(pathget)[names(pathget) %in% names(pathgets)]
    mapfile[,sampleid]<-tolower(mapfile[,sampleid])
    if (pathsum=="rel"){
      #calculate relative abundance
      pathrel<-as.data.frame(t(apply(pathtab[[j]], 1, function(x) x / sum(x)))[,pathname])
      pathrel[,sampleid]<-tolower(rownames(pathrel))
      pathdat<-merge(mapfile,pathrel,by=sampleid)
      #transformation of relative abundance
      if (stat.med=="gamlss" &transform!="none"){
        stop("gamlss with beta zero-inflated family should only be used for relative abundance without transformation")
      }
      if (stat.med=="lm" &transform=="asin.sqrt"){
        asintransform <- function(p) { asin(sqrt(p)) }
        pathdat[,pathname]<-apply(pathdat[,pathname],2,asintransform)
      }
      if (stat.med=="lm" &transform=="logit"){
        logittransform <- function(p) { log(p/(1-p)) }
        pathdat[,pathname]<-apply(pathdat[,pathname],2,logittransform )
      }
    }
    if (pathsum=="log"){
      # log2 transform
      pathlog<-log2(pathtab[[j]][,pathname]+1) #dirty handling of zero values
      pathlog[,sampleid]<-tolower(rownames(pathlog))
      pathdat<-merge(mapfile,pathlog,by=sampleid)
    }
    pathdat[,comvar]<-gdata::drop.levels(pathdat[,comvar],reorder=FALSE) #drop missing/unused level and keep level order
    if (longitudinal=="yes"){
      pathdat$personid<-as.factor(pathdat[,personid])
    }
    estisum<-NULL
    for (i in 1: length(pathname)){
      #print(i)
      #linear mixed model
      if (stat.med=="lm" & (pathsum=="rel"|pathsum=="log")){
        if (longitudinal=="yes"){
          fitsum<-try(summary(lme4::lmer(stats::as.formula(paste(pathname[i],paste(c(comvar,adjustvar,"(1|personid)"),collapse="+"),sep="~")), data=pathdat)))
          #fitsum<-try(summary(lme4::glmer(stats::as.formula(paste(pathname[i],paste(c(comvar,adjustvar,"(1|personid)"),collapse="+"),sep="~")), data=pathdat,family="gaussian")))
        }
        if (longitudinal=="no"){
          fitsum<-try(summary(stats::glm(stats::as.formula(paste(pathname[i],paste(c(comvar,adjustvar),collapse="+"),sep="~")), data=pathdat,family="gaussian")))
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
            fitcoef[,"id"]<-pathname[i]
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
            fitcoef[,"id"]<-pathname[i]
            fitcoefw<-stats::reshape(fitcoef, idvar="id", timevar="varname", direction="wide")
          }
          # when there is no coef
          if (length(which(rownames(fitsum$coefficients)!="(Intercept)"))==0){
            fitcoefw<-NULL
          }
          estisum<-plyr::rbind.fill(estisum,fitcoefw)
        }
      }
      #Generalized Additive Models for Location Scale and Shape: Betazeroinflated family, mu link logit
      if (stat.med=="gamlss" &(pathsum=="log")){
        stop("gamlss with beta zero-inflated family should only be used for relative abundance data")
      }
      if (stat.med=="gamlss" &(pathsum=="rel")){
        if (longitudinal=="yes"){
          testdat<-pathdat[,c(pathname[i],comvar,adjustvar,"personid")]
          testdat[,pathname[i]][testdat[,pathname[i]]==1]<-0.9999 # dirty fix for 1 value of relative abundance
          testdat<-stats::na.omit(testdat)
          fitsum<-try(summary(gamlss::gamlss(stats::as.formula(paste(pathname[i],paste(c(comvar,adjustvar,"random(personid)"),collapse="+"),sep="~")), family = BEZI, data = testdat, trace = FALSE),save=TRUE))
          #To check
          #fitsum<-try(summary(gamlss::gamlss(stats::as.formula(paste(pathname[i],paste(c(comvar,adjustvar,"gamlss::random(personid)"),collapse="+"),sep="~")), gamlss.family = BEZI, data = testdat, trace = FALSE),save=TRUE))
        }
        if (longitudinal=="no"){
          testdat<-pathdat[,c(pathname[i],comvar,adjustvar)]
          testdat[,pathname[i]][testdat[,pathname[i]]==1]<-0.9999 # dirty fix for 1 value of relative abundance
          testdat<-stats::na.omit(testdat)
          fitsum<-try(summary(gamlss::gamlss(stats::as.formula(paste(pathname[i],paste(c(comvar,adjustvar),collapse="+"),sep="~")), family = BEZI, data = testdat, trace = FALSE),save=TRUE))
          #TO check
          #fitsum<-try(summary(gamlss::gamlss(stats::as.formula(paste(pathname[i],paste(c(comvar,adjustvar),collapse="+"),sep="~")), gamlss.family = BEZI, data = testdat, trace = FALSE),save=TRUE))
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
            fitcoef[,"id"]<-pathname[i]
            fitcoefw<-stats::reshape(fitcoef, idvar="id", timevar="varname", direction="wide")
          }
          #handle the issues when there is only one row
          if (length(which(rownames(fitsum$coef.table)!="(Intercept)"))==1){
            fitcoef<-as.data.frame(matrix(fitsum$coef.table[rownames(fitsum$coef.table)!="(Intercept)",],ncol=ncol(fitsum$coef.table)))
            rownames(fitcoef)<-rownames(fitsum$coef.table)[rownames(fitsum$coef.table)!="(Intercept)"]
            colnames(fitcoef)<-colnames(fitsum$coef.table)
            fitcoef[,"varname"]<-rownames(fitcoef)
            fitcoef[,"id"]<-pathname[i]
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
    estisum[,sub('.*\\.', 'pval.adjust.',colnames(estisum)[grep("Pr(>|t|)",colnames(estisum))])]<-apply(estisum[,colnames(estisum)[grep("Pr(>|t|)",colnames(estisum))]],2,stats::p.adjust,method = p.adjust.method)
    estilist[[j]]<-estisum
  }
  names(estilist)<-pathlev
  return(estilist)
}
