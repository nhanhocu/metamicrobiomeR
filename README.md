


This is an initial version of the package *metamicrobiomeR*.   

To start, please download and view the steps/examples/workflow in the tutorial file "metamicrobiomeR_Supplementary_file". 

The *metamicrobiomeR* package implements Generalized Additive Model for Location, Scale and Shape (GAMLSS) 
    with zero inflated beta (BEZI) family for analysis of microbiome relative abundance data (with various         options for data transformation/normalization to address compositional effects) and 
    random effect meta-analysis models for meta-analysis pooling estimates across microbiome studies.
    Random Forest model to predict microbiome age based on relative abundances of  
    shared bacterial genera with the Bangladesh data (Subramanian et al 2014), 
    comparison of multiple diversity indexes using linear/linear mixed effect models 
    and some data display/visualization are also implemented.

The initial version of the manuscript describing the methods implemented in the *metamicrobiomeR* package as well as simulations, and example analyses/results are now available in preprint at: <https://www.biorxiv.org/content/early/2018/04/04/294678>. 


## Some snapshots

### Example 1: Comparison between breastfeeding statuses in infants < 6 months of age
#### Plot of mean relative abundance by breastfeeding statuses and age at phylum level  


```r
rm(list=ls()) # clear all
library(devtools)
#install and load package metamicrobiomeR
install_github("nhanhocu/metamicrobiomeR")
library(metamicrobiomeR) 
#Load other needed packages 
library(knitr)
library(plyr)
library(dplyr)
library(gdata)
library(gridExtra)
library(ggplot2)
library(lme4) 
library(lmerTest)
library(mgcv) 
library(meta) 

data(taxtab.rm7)
taxlist.rm<-taxa.filter(taxtab=taxtab.rm[[5]],percent.filter = 0.05, relabund.filter = 0.00005)
taxa.meansdn.rm<-taxa.meansdn(taxtab=taxtab.rm[[5]],sumvar="bf",groupvar="age.sample")
taxa.meansdn.rm<-taxa.meansdn.rm[taxa.meansdn.rm$bf!="No_BF" &taxa.meansdn.rm$age.sample<=6,]
taxa.meansdn.rm$bf<-drop.levels(taxa.meansdn.rm$bf,reorder=FALSE)
#phylum
p.bf.l2<-taxa.mean.plot(tabmean=taxa.meansdn.rm,tax.lev="l2", comvar="bf", groupvar="age.sample",mean.filter=0.005, show.taxname="short")
p.bf.l2$p
```

#### Comparison between breastfeeding statuses adjusting for age of infants at sample collection using GAMLSS 

```r
# Comparison of bacterial taxa relative abundance using GAMLSS 
taxacom6.zi.rmg<-taxa.compare(taxtab=taxtab6.rm,propmed.rel="gamlss",comvar="bf",adjustvar="age.sample",longitudinal="yes",p.adjust.method="fdr")
#phylum
kable(taxcomtab.show(taxcomtab=taxacom6.zi.rmg,tax.select=p.bf.l2$taxuse.rm, showvar="bfNon_exclusiveBF", tax.lev="l2",readjust.p=TRUE,p.adjust.method="fdr",p.cutoff = 1))
```


### Example 2: Comparison of bacterial taxa relative abundance in infants < 6 months between gender adjusting for breastfeeding statuses and age of infants at sample collection with GAMLSS 

```r
# Comparison of bacterial taxa relative abundance up to genus level 
taxacom6.zi.rm.sex.adjustbfage<-taxa.compare(taxtab=taxtab6.rm[[5]],propmed.rel="gamlss",comvar="gender",adjustvar=c("bf","age.sample"),longitudinal="yes")
#phylum
kable(taxcomtab.show(taxcomtab=taxacom6.zi.rm.sex.adjustbfage,tax.select="none", showvar="genderMale", tax.lev="l2",p.adjust.method="fdr"))
#order
kable(taxcomtab.show(taxcomtab=taxacom6.zi.rm.sex.adjustbfage,tax.select="none", showvar="genderMale", tax.lev="l4",p.adjust.method="fdr"))
#family
kable(taxcomtab.show(taxcomtab=taxacom6.zi.rm.sex.adjustbfage,tax.select="none", showvar="genderMale", tax.lev="l5",p.adjust.method="fdr"))
#genus
kable(taxcomtab.show(taxcomtab=taxacom6.zi.rm.sex.adjustbfage,tax.select="none", showvar="genderMale", tax.lev="l6",p.adjust.method="fdr"))
```

The analysis for other studies was done similarly. 

#### Meta-analysis of four studies (Bangladesh, Haiti, USA(CA_FL), USA(NC))
##### Load and combine saved results of four studies for the comparison of bacterial taxa relative abundance between genders adjusted for breastfeeding and infant age at sample collection.   

```r
data(taxacom.rm.sex.adjustbfage)
data(taxacom.ha.sex.adjustbfage)
data(taxacom6.zi.usbmk.sex.adjustbfage)
data(taxacom6.unc.sex.adjustedbfage)
taxacom6.zi.rm.sex.adjustbfage$study<-"Subramanian et al 2014 (Bangladesh)"
taxacom6.zi.rm.sex.adjustbfage$pop<-"Bangladesh"
taxacom.zi.ha.sex.adjustbfage$study<-"Bender et al 2016 (Haiti)"
taxacom.zi.ha.sex.adjustbfage$pop<-"Haiti"
taxacom6.zi.usbmk.sex.adjustbfage$study<-"Pannaraj et al 2017 (USA(CA_FL))"
taxacom6.zi.usbmk.sex.adjustbfage$pop<-"USA(CA_FL)"
taxacom6.zi.unc.sex.adjustedbfage$study<-"Thompson et al 2015 (USA(NC))"
taxacom6.zi.unc.sex.adjustedbfage$pop<-"USA(NC)"
tabsex4<-rbind.fill(taxacom6.zi.rm.sex.adjustbfage,taxacom.zi.ha.sex.adjustbfage,taxacom6.zi.usbmk.sex.adjustbfage,taxacom6.zi.unc.sex.adjustedbfage)
```

##### Meta-analysis 
```r
metab.sex<-meta.taxa(taxcomdat=tabsex4,summary.measure="RR",pool.var="id",studylab="study",backtransform=FALSE,percent.meta=0.5,p.adjust.method="fdr")
```
##### Display results as tables and figures. 
Phylum
```r
#table 
kable(metatab.show(metatab=metab.sex$random,com.pooled.tab=tabsex4,tax.lev="l2",showvar="genderMale",p.cutoff.type="p", p.cutoff=0.05,display="table"))
#plot
metadat<-metatab.show(metatab=metab.sex$random,com.pooled.tab=tabsex4,tax.lev="l2",showvar="genderMale",p.cutoff.type="p", p.cutoff=1,display="data")
meta.niceplot(metadat=metadat,sumtype="taxa",level="main",p="p",p.adjust="p.adjust",phyla.col="rainbow",forest.col="by.estimate",leg.key.size=1,leg.text.size=8,heat.text.x.size=7,heat.text.x.angle=0,forest.axis.text.y=8,forest.axis.text.x=7, point.ratio = c(4,2),line.ratio = c(2,1))
```
Family
```r
#table 
kable(metatab.show(metatab=metab.sex$random,com.pooled.tab=tabsex4,tax.lev="l5",showvar="genderMale",p.cutoff.type="p", p.cutoff=0.05,display="table"))
#Plot with different color palette, heatmap-forest width ratio
metadat<-metatab.show(metatab=metab.sex$random,com.pooled.tab=tabsex4,tax.lev="l5",showvar="genderMale",p.cutoff.type="p", p.cutoff=1,display="data")
meta.niceplot(metadat=metadat,sumtype="taxa",level="sub",p="p",p.adjust="p.adjust",phyla.col="rainbow",heat.forest.width.ratio =c(1,1.5), neg.palette = "Greens",pos.palette = "Purples",p.sig.heat="yes",leg.key.size=1,leg.text.size=8,heat.text.x.size=7,forest.axis.text.y=8,forest.axis.text.x=7)
```
Genus 
```r
#table 
kable(metatab.show(metatab=metab.sex$random,com.pooled.tab=tabsex4,tax.lev="l6",showvar="genderMale",p.cutoff.type="p", p.cutoff=0.05,display="table"))
# Plot with some different options: pooled estimates in forest plot with the same color scales as heatmap, those with p-values<0.05 in bold, FDR adjusted p-values<0.1 in triangles
metadat<-metatab.show(metatab=metab.sex$random,com.pooled.tab=tabsex4,tax.lev="l6",showvar="genderMale",p.cutoff.type="p", p.cutoff=1,display="data")
meta.niceplot(metadat=metadat,sumtype="taxa",level="sub",p="p",p.adjust="p.adjust",phyla.col="rainbow",p.sig.heat="yes",heat.forest.width.ratio =c(1,1.3),forest.col="by.estimate",leg.key.size=0.8,leg.text.size=10,heat.text.x.size=6,forest.axis.text.y=7,forest.axis.text.x=6, point.ratio = c(4,2),line.ratio = c(2,1))
```
### Alpha diversity 
#### Calculate mean alpha diversity indexes for a selected rarefaction depth, standardize and compare standardized alpha diversity indexes between groups adjusting for covariates using Bangladesh data 

```r
data(sam.rm)
patht<-system.file("extdata/QIIME_outputs/Bangladesh/alpha_div_collated", package = "metamicrobiomeR", mustWork = TRUE)
alpha.rm<-read.multi(patht=patht,patternt=".txt",assignt="no",study="Bangladesh")
names(alpha.rm)<-sub(patht,"",names(alpha.rm))
samfile<-merge(samde, he50[,c("child.id","gender","month.exbf","month.food")],by="child.id")
samfile$age.sample<-samfile$age.months
samfile$bf<-factor(samfile$bf,levels=c("ExclusiveBF","Non_exclusiveBF","No_BF"))
samfile$personid<-samfile$child.id
samfile$sampleid<-tolower(samfile$fecal.sample.id)
#comparison of standardized alpha diversity indexes between genders adjusting for breastfeeding and infant age at sample collection in infants <=6 months of age 
alphacom6.rm.sexsg<-alpha.compare(datlist=alpha.rm,depth=3,mapfile=samfile,mapsampleid="fecal.sample.id",comvar="gender",adjustvar=c("age.sample","bf"),longitudinal="yes",age.limit=6,standardize=TRUE)
kable(alphacom6.rm.sexsg$alphasum[,1:5])
```

### Microbiome age
#### Predicting microbiome age, checking model performance, and replicate the results of the Bangladesh study
The *microbiomeage* function get the shared genera list between the Bangladesh study and all other included studies,  get the training and test sets from Bangladesh data based on the shared genera list, fit the train Random Forest model and predict microbiome age in the test set of Bangladesh data and data from all included studies, check for performance of the model based on the shared genera list on Bangladesh healthy cohort data, reproduce the findings of the Bangladesh malnutrition study.   

```r
#load Bangladesh taxa relative abundance summary up to genus level merged with mapping file (output from QIIME)
bal6<-read.delim(system.file("extdata/QIIME_outputs/Bangladesh/tax_mapping7", "Subramanian_et_al_mapping_file_L6.txt", package = "metamicrobiomeR", mustWork = TRUE))
colnames(bal6)<-tolower(colnames(bal6))
#View(bal6)
#format for data of other studies should be similar to Bangladesh data, must have 'age.sample' variable as age of infant at stool sample collection 
# Load data of 3 other studies 
data(gtab.3stud)
names(gtab.3stud)
#predict microbiome age on Bangladesh data and data of other three studies based on shared genera across 4 studies  
#(take time to run)
miage<-microbiomeage(l6.relabundtab=gtab.3stud)
# list of shared genera that are available in the Bangladesh study and other included studies 
kable(miage$sharedgenera.importance)
#check performance
grid.arrange(miage$performanceplot$ptrain, miage$performanceplot$ptest,nrow=1)
```


