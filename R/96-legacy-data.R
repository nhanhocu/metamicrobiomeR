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
#' (\href{https://pmc.ncbi.nlm.nih.gov/articles/PMC4189846/}{PubMed Central})
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

#' Combined alpha diversity data for meta-analysis.
#'
#' Result outputs of differential analysis of alpha diversity
#' using linear or linear mixed effects model
#' from "alpha.compare" function combined from
#' 4 studies for meta-analysis. The comparison was between gender
#' adjusted for age of infants at sample collection.
#'
#' @docType data
#'
#' @usage data(asum4)
#'
#' @format A dataframe with 16 rows and 18 variables.
#'
#' @keywords datasets
#'
#' @references Subramanian et al. Nature. 2014 Jun 19; 510(7505): 417–421.
#' (\href{https://pmc.ncbi.nlm.nih.gov/articles/PMC4189846/}{PubMed Central})
#'
#' Bender et al. Sci Transl Med. 2016 Jul 27; 8(349): 349ra100.
#' (\href{https://pmc.ncbi.nlm.nih.gov/articles/PMC5301310/}{PubMed Central})
#'
#' Pannaraj et al. JAMA Pediatr. 2017;90095(7):647–54.
#' (\doi{10.1001/jamapediatrics.2017.0378})
#'
#' Thompson et al. Front Cell Infect Microbiol. 2015;5:3.
#' (\doi{10.3389/fcimb.2015.00003})
#'
#' @source \href{https://gordonlab.wustl.edu/supplemental-data/supplemental-data-portal/subramanian-et-al-2014/}{Gordon Lab}
#'
#' @examples
#' data(asum4)
"asum4"

#' Covariate data.
#'
#' Monthly longitudinal clinical data of 50 infants from
#' birth to 2 years of life.
#'
#' @docType data
#'
#' @usage data(covar.rm)
#'
#' @format A dataframe with 996 rows and 32 variables.
#'
#' @keywords datasets
#'
#' @references Subramanian et al. Nature. 2014 Jun 19; 510(7505): 417–421.
#' (\href{https://pmc.ncbi.nlm.nih.gov/articles/PMC4189846/}{PubMed Central})
#'
#' @source \href{https://gordonlab.wustl.edu/supplemental-data/supplemental-data-portal/subramanian-et-al-2014/}{Gordon Lab}
#'
#' @examples
#' data(covar.rm)
#' # Load KEGG pathway data
#' data(kegg.12)
#' # Comparison of pathway relative abundances for some first pathways of level 1 only
#' # and assuming crosssectional data (to save running time)
#' path1<-pathway.compare(pathtab=list(kegg.12[[1]][, 1:2]),
#' mapfile=covar.rm,sampleid="sampleid",pathsum="rel", stat.med="gamlss",
#' comvar="gender",adjustvar=c("age.sample","bf"), longitudinal="no",
#' p.adjust.method="fdr", percent.filter=0.05,relabund.filter=0.00005)
#' taxcomtab.show(taxcomtab=path1$l1, sumvar="path",tax.lev="l2",
#' tax.select="none", showvar="genderMale", p.adjust.method="fdr",p.cutoff=1)
"covar.rm"

#' Pathway abundance data.
#'
#' KEGG pathway abundance data from PICRUSt analysis.
#' This is monthly longitudinal data of 50 infants from birth to
#' 2 years of life.
#'
#' @docType data
#'
#' @usage data(kegg.12)
#'
#' @format A list of 2 dataframes for level 1 and level 2 of KEGG pathways.
#'
#' @keywords datasets
#'
#' @references Subramanian et al. Nature. 2014 Jun 19; 510(7505): 417–421.
#' (\href{https://pmc.ncbi.nlm.nih.gov/articles/PMC4189846/}{PubMed Central})
#'
#' @source \href{https://gordonlab.wustl.edu/supplemental-data/supplemental-data-portal/subramanian-et-al-2014/}{Gordon Lab}
#'
#' @examples
#' data(kegg.12)
#' # Load covariate data
#' data(covar.rm)
#' # Comparison of pathway relative abundances for some first pathways of level 1 only
#' # and assuming crosssectional data (to save running time)
#' path1<-pathway.compare(pathtab=list(kegg.12[[1]][, 1:2]),
#' mapfile=covar.rm,sampleid="sampleid",pathsum="rel", stat.med="gamlss",
#' comvar="gender",adjustvar=c("age.sample","bf"), longitudinal="no",
#' p.adjust.method="fdr", percent.filter=0.05,relabund.filter=0.00005)
#' taxcomtab.show(taxcomtab=path1$l1, sumvar="path",tax.lev="l2",
#' tax.select="none", showvar="genderMale", p.adjust.method="fdr",p.cutoff=1)
"kegg.12"

#' Combined data for meta-analysis.
#'
#' Result outputs of differential abundance analysis using GAMLSS_BEZI
#' from "taxa.compare" function combined from
#' 4 studies for meta-analysis. The comparison was between gender
#' adjusted for age of infants at sample collection.
#'
#' @docType data
#'
#' @usage data(tabsex4)
#'
#' @format A dataframe with 701 rows and 23 variables.
#'
#' @keywords datasets
#'
#' @references Subramanian et al. Nature. 2014 Jun 19; 510(7505): 417–421.
#' (\href{https://pmc.ncbi.nlm.nih.gov/articles/PMC4189846/}{PubMed Central})
#'
#' Bender et al. Sci Transl Med. 2016 Jul 27; 8(349): 349ra100.
#' (\href{https://pmc.ncbi.nlm.nih.gov/articles/PMC5301310/}{PubMed Central})
#'
#' Pannaraj et al. JAMA Pediatr. 2017;90095(7):647–54.
#' (\doi{10.1001/jamapediatrics.2017.0378})
#'
#' Thompson et al. Front Cell Infect Microbiol. 2015;5:3.
#' (\doi{10.3389/fcimb.2015.00003})
#'
#' @source \href{https://gordonlab.wustl.edu/supplemental-data/supplemental-data-portal/subramanian-et-al-2014/}{Gordon Lab}
#'
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
#' #phylum table
#' metatab.show(metatab=metab.sex$random,com.pooled.tab=tabsex4[tabsex4$id %in% tlm,],
#' tax.lev="l2",showvar="genderMale",p.cutoff.type="p", p.cutoff=1,display="table")
"tabsex4"

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
#' @format A data frame with 322 rows (samples) and 803 variables (including mapping variables and bacterial taxonomies from phylum to genus level).
#'
#' @keywords datasets
#'
#' @references Subramanian et al. Nature. 2014 Jun 19; 510(7505): 417–421.
#' (\href{https://pmc.ncbi.nlm.nih.gov/articles/PMC4189846/}{PubMed Central})
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
