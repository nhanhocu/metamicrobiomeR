## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,cache=FALSE,message=FALSE,warning=FALSE,echo = TRUE,
  comment = "#>"
)

## ----pack, message=FALSE, comment=FALSE,results='hide'------------------------
library(metamicrobiomeR) 

## ----bacrelmeanba, message=FALSE, comment=FALSE, fig.width=10, fig.height=8----
data(taxtab6)
taxlist.rm<-taxa.filter(taxtab=taxtab6,percent.filter = 0.05, relabund.filter = 0.00005)
taxa.meansdn.rm<-taxa.meansdn(taxtab=taxtab6,sumvar="bf",groupvar="age.sample")
taxa.meansdn.rm<-taxa.meansdn.rm[taxa.meansdn.rm$bf!="No_BF",] #&taxa.meansdn.rm$age.sample<=6,
taxa.meansdn.rm$bf<-gdata::drop.levels(taxa.meansdn.rm$bf,reorder=FALSE)
#phylum
p.bf.l2<-taxa.mean.plot(tabmean=taxa.meansdn.rm,tax.lev="l2", comvar="bf", groupvar="age.sample",mean.filter=0.005, show.taxname="short")
p.bf.l2$p

## ----bacrelcomgamlss, results="hide"------------------------------------------
# Note: running time is not long in regular laptop for both GAMLSS-BEZI analysis (~10s) and meta-analysis (~5s).  
# However, to save running time, only taxonomies of one small phylum are selected for differential analysis example. 
tab6<-as.data.frame(taxtab6)
tl<-colnames(taxtab6)[grep("k__bacteria.p__fusobacteria",colnames(taxtab6))]
taxacom.ex<-taxa.compare(taxtab=tab6[,c("personid","x.sampleid","bf","age.sample",tl)],propmed.rel="gamlss",comvar="bf",adjustvar="age.sample",
longitudinal="yes",p.adjust.method="fdr")

## ----bacrelcomgamlssshow------------------------------------------------------
#phylum
taxcomtab.show(taxcomtab=taxacom.ex,tax.select="none", showvar="bfNon_exclusiveBF", tax.lev="l2",readjust.p=TRUE,p.adjust.method="fdr",p.cutoff = 1)
#genus
taxcomtab.show(taxcomtab=taxacom.ex,tax.select="none", showvar="bfNon_exclusiveBF", tax.lev="l6", readjust.p=TRUE,p.adjust.method="fdr",p.cutoff = 1)

## ----bacrellmeas--------------------------------------------------------------
taxacom.lmas<-taxa.compare(taxtab=tab6[,c("personid","x.sampleid","bf","age.sample",tl)], propmed.rel="lm",transform="asin.sqrt",comvar="bf",adjustvar="age.sample", longitudinal="yes",p.adjust.method="fdr")
#phylum
taxcomtab.show(taxcomtab=taxacom.lmas,tax.select="none", showvar="bfNon_exclusiveBF", tax.lev="l2",readjust.p=TRUE,p.adjust.method="fdr",p.cutoff = 1, digit=5,p.digit=5) 
#family
taxcomtab.show(taxcomtab=taxacom.lmas,tax.select="none", showvar="bfNon_exclusiveBF", tax.lev="l5",readjust.p=TRUE,p.adjust.method="fdr",p.cutoff = 1, digit=5,p.digit=5)

## ----keggba, results="hide"---------------------------------------------------
data(kegg.12)
data(covar.rm)
# Comparison of pathway relative abundances for level 1 only (to save running time)
path1<-pathway.compare(pathtab=list(kegg.12[[1]]),mapfile=covar.rm,sampleid="sampleid", pathsum="rel",stat.med="gamlss",comvar="gender",adjustvar=c("age.sample","bf"), longitudinal="yes",p.adjust.method="fdr",percent.filter=0.05,relabund.filter=0.00005)

## ----keggbar------------------------------------------------------------------
taxcomtab.show(taxcomtab=path1$l1, sumvar="path", tax.lev="l2",tax.select="none",showvar="genderMale", p.adjust.method="fdr",p.cutoff=1)

## ----bacrelmeta, fig.width=10, fig.height=8-----------------------------------
# load saved GAMLSS-BEZI results of four studies for the comparison of bacterial taxa relative abundance between genders adjusted for breastfeeding and infant age at sample collection 
data(tabsex4)

#select only taxonomies of a small phylum for meta-analysis example (to save running time)
tlm<-tabsex4$id[grep("k__bacteria.p__fusobacteria",tabsex4$id)]
# meta-analysis 
metab.sex<-meta.taxa(taxcomdat=tabsex4[tabsex4$id %in% tlm,], summary.measure="RR", pool.var="id", studylab="study", backtransform=FALSE, percent.meta=0.5, p.adjust.method="fdr")
#show results by table and plot
#phylum
#table
metatab.show(metatab=metab.sex$random,com.pooled.tab=tabsex4[tabsex4$id %in% tlm,], tax.lev="l2",showvar="genderMale",p.cutoff.type="p", p.cutoff=1,display="table")
#plot
metadat<-metatab.show(metatab=metab.sex$random,com.pooled.tab=tabsex4, tax.lev="l2",showvar="genderMale",p.cutoff.type="p", p.cutoff=1,display="data")
meta.niceplot(metadat=metadat,sumtype="taxa",level="main",p="p", p.adjust="p.adjust",phyla.col="rainbow",p.sig.heat="yes", heat.forest.width.ratio =c(1.5,1), leg.key.size=0.8, leg.text.size=10, heat.text.x.size=10, heat.text.x.angle=0, forest.axis.text.y=8,forest.axis.text.x=10, point.ratio = c(4,2),line.ratio = c(2,1))

## ----alpharm, fig.width=10, fig.height=8--------------------------------------
data(alphadat)
data(covar.rm)
covar.rm$sampleid<-tolower(covar.rm$sampleid)

#comparison of standardized alpha diversity indexes between genders adjusting for breastfeeding and infant age at sample collection in infants <=6 months of age 
alphacom<-alpha.compare(datlist=alphadat,depth=3,mapfile=covar.rm, mapsampleid="sampleid",comvar="gender",adjustvar=c("age.sample","bf"), longitudinal="yes",age.limit=6,standardize=TRUE)
alphacom$alphasum[,1:5] 

## ----shameta, fig.width=10, fig.height=5--------------------------------------
# load saved results of 4 studies 
data(asum4)
asum4[,c(colnames(asum4)[1:5],"pop")]
#Shannon index 
shannon.sex <- meta::metagen(Estimate.genderMale, `Std. Error.genderMale`, studlab=pop,data=subset(asum4,id=="shannon"),sm="RD", backtransf=FALSE)
meta::forest(shannon.sex,smlab="Standardized \n diversity difference",sortvar=subset(asum4,id=="shannon")$pop,lwd=2)
shannon.sex
cbind(study=shannon.sex$studlab,pval=shannon.sex$pval)

