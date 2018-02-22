#' Summarize abundance by group
#'
#' This function summarizes taxa/pathway abundance tables to provide mean, sd, count by groups.
#' @param taxtab taxa/pathway  abundance table from phylum to species or any preferred highest taxa level.
#' @param sumvar main variable for summary
#' @param groupvar variable to be stratified.
#' @param percent.filter prevalence threshold (the percentage of number of samples the taxa/pathway available). Default is 0.05.
#' @param relabund.filter relative abundance threshold (the minimum of the average relative abundance for a taxa/pathway to be retained). Default is 0.00005.
#' @param othervar list of variables that are not abundance variables to be summarized. Default is "none".
#' @return table of mean, sd, count by group.
#' @keywords abundance summary
#' @export
#' @examples
#' #Load summary tables of bacterial taxa relative abundance from Bangladesh data
#' data(taxtab.rm7)
#' library(magrittr)
#' taxa.meansdn.rm<-taxa.meansdn(taxtab=taxtab.rm[[5]],sumvar="bf",groupvar="age.sample")

taxa.meansdn<-function(taxtab, sumvar, groupvar,percent.filter=0.05,relabund.filter=0.00005,othervar="none"){
  #require(dplyr)
  taxdat<-as.data.frame(taxtab)
  taxdat$sumvar<-taxdat[,sumvar]
  taxdat$groupvar<-taxdat[,groupvar]
  # get assigned taxa only
  if (othervar!="none"){
    taxlist<-colnames(taxdat)[!colnames(taxdat) %in% othervar]
  }
  if (othervar=="none") {
    taxlist<-colnames(taxdat)[grep("k__",colnames(taxdat))]
  }
  #filter using percent.filter
  taxtest<-apply(taxdat[,taxlist],2,function(x){length(x[!is.na(x)&x>0])})
  taxget<-taxtest[taxtest>=percent.filter*(nrow(taxdat))]
  #filter using relabund.filter
  taxtestm<-apply(taxdat[,taxlist],2,mean,na.rm=T)
  taxgetm<-taxtestm[taxtestm>relabund.filter]
  taxname<-names(taxget)[names(taxget) %in% names(taxgetm)]
  sumdat<-taxdat[,c("sumvar", "groupvar",taxname)]
  sumdat[,taxname]<-lapply(sumdat[,taxname],as.character)
  sumdat[,taxname]<-lapply(sumdat[,taxname],as.numeric)
  sumdat<-sumdat[!is.na(sumdat[,"sumvar"]),]
  if (is.numeric(sumdat$groupvar)){
    estisum<-sumdat %>%
      dplyr::mutate(groupvar=as.factor(as.character(round(as.numeric(as.character(groupvar)),0)))) %>%
      dplyr::group_by(sumvar, groupvar) %>%
      dplyr::summarise_all(dplyr::funs(mean(.), sd(.), n()))
    estisum<-estisum%>%
      dplyr::mutate(groupvar=as.numeric(as.character(groupvar))) %>%
      dplyr::arrange(sumvar,groupvar)
    estisum<-na.omit(estisum)
  }
  else {
    estisum<-sumdat %>%
      dplyr::group_by(sumvar, groupvar) %>%
      dplyr::summarise_all(dplyr::funs(mean(.), sd(.), n()))
  }
  colnames(estisum)[colnames(estisum) %in% c("sumvar","groupvar")]<-c(sumvar,groupvar)
  return(estisum)
}
