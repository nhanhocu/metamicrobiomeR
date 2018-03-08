#' Filter relative abundance data
#'
#' This function filters bacterial taxa or pathway relative abundance tables based on the percentage of samples with their availability (prevalence) and relative abundance thresholds.
#' It will remove taxa/pathway with relative abundance <relabund.filter and available in <percent.filter of number of samples.
#' @param taxtab taxa/pathway relative abundance table.
#' @param percent.filter prevalence threshold (the percentage of number of samples the taxa/pathway available). Default is 0.05.
#' @param relabund.filter relative abundance threshold (the minimum of the average relative abundance for a taxa/pathway to be retained). Default is 0.00005.
#' @return list of all taxa/pathways retained after filtering.
#' @keywords relative abundance filter
#' @export
#' @examples
#' #Load summary tables of bacterial taxa relative abundance from Bangladesh data
#' data(taxtab.rm7)
#' taxlist.rm<-taxa.filter(taxtab=taxtab.rm[[5]],percent.filter = 0.05, relabund.filter = 0.00005)

taxa.filter<-function(taxtab, percent.filter=0.05, relabund.filter=0.00005){
  # filter (remove) taxa with relative abundance <relabund.filter and available in <percent.filter of number of samples
  taxdat<-as.data.frame(taxtab)
  # get assigned taxa only
  taxlist<-colnames(taxdat)[grep("k__",colnames(taxdat))]
  #filter using percent.filter
  taxtest<-apply(taxdat[,taxlist],2,function(x){length(x[!is.na(x)&x>0])})
  taxget<-taxtest[taxtest>=percent.filter*(nrow(taxdat))]
  #filter using relabund.filter
  taxtestm<-apply(taxdat[,taxlist],2,mean,na.rm=T)
  taxgetm<-taxtestm[taxtestm>relabund.filter]
  taxlistf<-c(names(taxget)[names(taxget) %in% names(taxgetm)])
  return(taxlistf)
}
