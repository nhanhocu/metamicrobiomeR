#' Alpha diversity data.
#'
#' Alpha diversity data from alpha_rarefaction output from QIIME.
#'
#'
#' @docType data
#'
#' @usage data(alphadat)
#'
#' @format A list of 4 dataframes for four indexes: Chao1, Observed_species, PD_whole_tree, Shannon.
#'
#' @keywords datasets
#'
#' @references Subramanian et al. Nature. 2014 Jun 19; 510(7505): 417–421.
#' (\href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4189846/}{PubMed})
#'
#' @source \href{https://gordonlab.wustl.edu/supplemental-data/supplemental-data-portal/subramanian-et-al-2014/}{Gordon Lab}
#'
#' @examples
#' data(alphadat)
#' # Load covariate data
#' data(covar.rm)
#' covar.rm$sampleid<-tolower(covar.rm$sampleid)
#' #comparison of standardized alpha diversity indexes between genders adjusting for
#' #breastfeeding and infant age at sample collection in infants <=6 months of age
#' alphacom<-alpha.compare(datlist=alphadat,depth=3,mapfile=covar.rm,
#' mapsampleid="sampleid", comvar="gender",adjustvar=c("age.sample","bf"),
#' longitudinal="yes", age.limit=6,standardize=TRUE)
#' alphacom$alphasum[,1:5]
"alphadat"
