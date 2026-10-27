#' Taxonomic relative abundance data.
#'
#' Monthly longitudinal relative abundance data from phylum to genus level of
#' 50 infants from birth to 2 year of life.
#' Mapping file is merged to the data for ready use.
#'
#' @docType data
#'
#' @usage data(taxtab6)
#'
#' @format A data frame with 322 row (samples) and 803 variables (including mapping varilable and bacterial taxonomies from phylum to genus level).
#'
#' @keywords datasets
#'
#' @references Subramanian et al. Nature. 2014 Jun 19; 510(7505): 417–421.
#' (\href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4189846/}{PubMed})
#'
#' @source \href{https://gordonlab.wustl.edu/supplemental-data/supplemental-data-portal/subramanian-et-al-2014/}{Gordon Lab}
#'
#' @examples
#' #Load summary tables of bacterial taxa relative abundance from Bangladesh data
#' data(taxtab6)
#' tab6<-as.data.frame(taxtab6)
#' tl<-colnames(taxtab6)[grep("k__bacteria.p__fusobacteria",colnames(taxtab6))]
#' taxacom.ex<-taxa.compare(taxtab=tab6[,c("personid","x.sampleid","bf","age.sample",tl)],
#' propmed.rel="gamlss",comvar="bf",adjustvar="age.sample",
#' longitudinal="yes",p.adjust.method="fdr")
"taxtab6"
