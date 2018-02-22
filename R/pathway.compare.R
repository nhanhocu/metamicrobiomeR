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
#' @param comvar main variable for comparison.
#' @param adjustvar variables to be adjusted.
#' @param personid name of variable for person id in mapping file (applicable for longitudinal data)
#' @param longitudinal whether the data is longitudinal. Default is "yes".
#' @param p.adjust.method method for multiple testing adjustment. Available options are those of the p.adjust function. Default is "fdr".
#' @param percent.filter prevalence threshold (the percentage of number of samples the taxa/pathway available). Default is 0.05.
#' @param relabund.filter relative abundance threshold (the minimum of the average relative abundance for a taxa/pathway to be retained). Default is 0.00005.
#' @param pooldata whether the data is pooled from multiple studies. Default is FALSE.
#' @param age.limit upper age limit for data to be included. Default is 1000000 months (~no upper age limit).
#' @param age.lowerlimit lower age limit for data to be included. Default is 0 month.
#' @return matrice of coefficients, standard errors, p-values and multiple testing adjusted p-values of all variables in the models.
#' @keywords pathway abundance comparison.
#' @export
#' @examples
#' #Load Bangladesh extra metadata
#' data(sam.rm)
#' #Read in PICRUSt output for KEGG pathways level 1-3
#' patht<-system.file("extdata/QIIME_outputs/Bangladesh/picrust", package = "metamicrobiomeR", mustWork = TRUE)
#' kegg<-read.multi(patht=patht,patternt=".txt",assignt="no")
#' kegg.rm<-list()
#' for (i in 1:length(kegg)){
#'   rownames(kegg[[i]])<-kegg[[i]][,"kegg_pathways"]
#'   kegg[[i]]<-kegg[[i]][,colnames(kegg[[i]])[!colnames(kegg[[i]]) %in% c("otu.id","kegg_pathways")]]
#'   kegg.rm[[i]]<-as.data.frame(t(kegg[[i]]))
#' }
#' covar.rm<-merge(samde, he50[,c("child.id","gender","zygosity","day.firstsample","day.lastsample","n.sample","sampling.interval.msd","month.exbf","month.food",
#'                                "n.diarrhea.yr","percent.time.diarrhea","fraction.antibiotic","subject.allocation")], by="child.id")
#' covar.rm<-dplyr::rename(covar.rm,sampleid=fecal.sample.id, personid=child.id ,age.sample=age.months)
#' covar.rm$bf<-factor(covar.rm$bf, levels=c('ExclusiveBF','Non_exclusiveBF','No_BF'))
#' covar.rm$personid<-as.factor(covar.rm$personid)
#' #Comparison of pathway relative abundances (take time to run)
#' pathcom.rm6.rel.gamlss.sexg<-pathway.compare(pathtab=kegg.rm,mapfile=covar.rm,sampleid="sampleid",pathsum="rel",stat.med="gamlss",comvar="gender",adjustvar=c("age.sample","bf"),longitudinal="yes",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,age.limit=6)
#' taxcomtab.show(taxcomtab=pathcom.rm6.rel.gamlss.sexg$l3, sumvar="path",tax.lev="l3",tax.select="none",showvar="genderMale", p.adjust.method="fdr",p.cutoff=0.05,digit=4)

pathway.compare<-function(pathtab,mapfile,sampleid="sampleid",pathsum="rel",stat.med="gamlss",comvar="bf",adjustvar="age.sample",personid="personid",longitudinal="yes",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005,pooldata=FALSE,age.limit=1000000,age.lowerlimit=0,...){
  #sapply(c("lme4","lmerTest","gamlss","gdata","plyr"), require, character.only = TRUE)
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
      pathdat<-merge(subset(mapfile,age.sample<=age.limit & age.sample>=age.lowerlimit),pathrel,by=sampleid)
    }
    if (pathsum=="log"){
      # log2 transform
      pathlog<-log2(pathtab[[j]][,pathname]+1) #dirty handling of zero values
      pathlog[,sampleid]<-tolower(rownames(pathlog))
      pathdat<-merge(subset(mapfile,age.sample<=age.limit & age.sample>=age.lowerlimit),pathlog,by=sampleid)
    }
    pathdat[,comvar]<-gdata::drop.levels(pathdat[,comvar],reorder=FALSE) #drop missing/unused level and keep level order
    if (longitudinal=="yes"){
      pathdat$personid<-as.factor(pathdat[,personid])
    }
    estisum<-NULL
    for (i in 1: length(pathname)){
      print(i)
      #linear mixed model
      if (stat.med=="lm" & (pathsum=="rel"|pathsum=="log")){
        if (longitudinal=="yes"){
          fitsum<-try(summary(lme4::glmer(as.formula(paste(pathname[i],paste(c(comvar,adjustvar,"(1|personid)"),collapse="+"),sep="~")), data=pathdat,family=gaussian(link="identity"))))
        }
        if (longitudinal=="no"){
          fitsum<-try(summary(glm(as.formula(paste(pathname[i],paste(c(comvar,adjustvar),collapse="+"),sep="~")), data=pathdat,family="gaussian")))
        }
        if (class(fitsum) == "try-error") {
          cat("Error in model fit, NA introduced.\n")
          fitcoefw<-NULL
          estisum<-plyr::rbind.fill(estisum,fitcoefw)
        }
        if (class(fitsum) != "try-error") {
          fitcoef<-as.data.frame(fitsum$coefficients[rownames(fitsum$coefficients)!="(Intercept)",]) #remove intercept
          if (longitudinal=="yes"){
            fitcoef[,"Pr(>|t|)"]<-2*pnorm(-abs(fitcoef[,"Estimate"]/fitcoef[,"Std. Error"]))
          }
          fitcoef[,"varname"]<-rownames(fitcoef)
          fitcoef[,"id"]<-pathname[i]
          fitcoefw<-reshape(fitcoef, idvar="id", timevar="varname", direction="wide")
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
          testdat<-na.omit(testdat)
          fitsum<-try(summary(gamlss::gamlss(as.formula(paste(pathname[i],paste(c(comvar,adjustvar,"gamlss::random(personid)"),collapse="+"),sep="~")), gamlss.family = BEZI, data = testdat, trace = FALSE),save=TRUE))
        }
        if (longitudinal=="no"){
          testdat<-pathdat[,c(pathname[i],comvar,adjustvar)]
          testdat[,pathname[i]][testdat[,pathname[i]]==1]<-0.9999 # dirty fix for 1 value of relative abundance
          testdat<-na.omit(testdat)
          fitsum<-try(summary(gamlss::gamlss(as.formula(paste(pathname[i],paste(c(comvar,adjustvar),collapse="+"),sep="~")), gamlss.family = BEZI, data = testdat, trace = FALSE),save=TRUE))
        }
        if (class(fitsum) == "try-error") {
          cat("Error in model fit, NA introduced.\n")
          fitcoefw<-NULL
          estisum<-plyr::rbind.fill(estisum,fitcoefw)
        }
        if (class(fitsum) != "try-error") {
          if (length(which(rownames(fitsum$coef.table)!="(Intercept)"))>1){
            fitcoef<-as.data.frame(fitsum$coef.table[rownames(fitsum$coef.table)!="(Intercept)",])
            fitcoef[,"varname"]<-rownames(fitcoef)
            fitcoef[,"id"]<-pathname[i]
            fitcoefw<-reshape(fitcoef, idvar="id", timevar="varname", direction="wide")
          }
          #handle the issues when there is only one row
          if (length(which(rownames(fitsum$coef.table)!="(Intercept)"))==1){
            fitcoef<-as.data.frame(matrix(fitsum$coef.table[rownames(fitsum$coef.table)!="(Intercept)",],ncol=ncol(fitsum$coef.table)))
            rownames(fitcoef)<-rownames(fitsum$coef.table)[rownames(fitsum$coef.table)!="(Intercept)"]
            colnames(fitcoef)<-colnames(fitsum$coef.table)
            fitcoef[,"varname"]<-rownames(fitcoef)
            fitcoef[,"id"]<-pathname[i]
            fitcoefw<-reshape(fitcoef, idvar="id", timevar="varname", direction="wide")
          }
          # when there is no coef
          if (length(which(rownames(fitsum$coef.table)!="(Intercept)"))==0){
            fitcoefw<-NULL
          }
          estisum<-plyr::rbind.fill(estisum,fitcoefw)
        }
      }
    }
    estisum[,sub('.*\\.', 'pval.adjust.',colnames(estisum)[grep("Pr(>|t|)",colnames(estisum))])]<-apply(estisum[,colnames(estisum)[grep("Pr(>|t|)",colnames(estisum))]],2,p.adjust,method = p.adjust.method)
    estilist[[j]]<-estisum
  }
  names(estilist)<-pathlev
  return(estilist)
}
