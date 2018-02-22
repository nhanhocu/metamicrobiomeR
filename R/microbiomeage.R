#' Predict microbiome age.
#'
#' This function predicts microbiome age using Random Forest model based on relative abundances of bacterial genera shared with the Bangladesh study (Subramanian et al 2014).
#' This function gets the shared genera list between the Bangladesh study and all other included studies,
#' get the training and test sets from Bangladesh data based on the shared genera list,
#' fit the train Random Forest model and predict microbiome age in the test set of Bangladesh data and data from all included studies,
#' check for performance of the model based on the shared genera list on Bangladesh healthy cohort data, reproduce the findings of the Bangladesh malnutrition study.
#' @param l6.relabundtab list of taxa summary table from phylum up to genus level merged to mapping file outputed from QIIME of all included studies.
#' @return list of training and test sets of Bangladesh data, shared genera list, relative abundance data of shared genera, randomforest fit, RF model performance plot,predicted microbiome age of Bangladesh data and data of other included studies.
#' @keywords microbiome age Random Forest.
#' @export
#' @examples
#' #Load data from each study and put in a list
#' bal6<-read.delim(system.file("extdata/QIIME_outputs/Bangladesh/tax_mapping7", "Subramanian_et_al_mapping_file_L6.txt", package = "metamicrobiomeR", mustWork = TRUE))
#' colnames(bal6)<-tolower(colnames(bal6))
#' uncl6<-read.delim(system.file("extdata/QIIME_outputs/USA(NC)/tax_mapping7", "mapping_fileall_L6.txt", package = "metamicrobiomeR", mustWork = TRUE))
#' colnames(uncl6)<-tolower(colnames(uncl6))
#' uncl6$age.sample<-uncl6$agemo
#' usl6<-read.delim(system.file("extdata/QIIME_outputs/USA(CA_FL)/tax_mapping7", "USBMK_Mapping_122414.filtered.with_pairing_status_L6.txt", package = "metamicrobiomeR", mustWork = TRUE))
#' colnames(usl6)<-tolower(colnames(usl6))
#' usl6$age.sample<-usl6$bbage/30
#' hal6<-read.delim(system.file("extdata/QIIME_outputs/Haiti/tax_mapping7", "Haiti_Mapping.complete_L6.txt", package = "metamicrobiomeR", mustWork = TRUE))
#' colnames(hal6)<-tolower(colnames(hal6))
#' hal6$age.sample<-hal6$bbage/30
#' gtab<-list(unc=uncl6,ca_fl=usl6[usl6$sampletype=="STL"&usl6$mombb=="Baby",],haiti=hal6[hal6$sampletype=="STL"&hal6$mombb=="B",])
#' #predict microbiome age (take time to run)
#' miage<-microbiomeage(l6.relabundtab=gtab)
#' #list of shared genera that are available in the Bangladesh study and other included studies
#' miage$sharedgenera.importance
#' #check performance
#' gridExtra::grid.arrange(miage$performanceplot$ptrain, miage$performanceplot$ptest,nrow=1)
#' #replicate the findings of Subramanian et al paper
#' ggplot2::ggplot() +geom_point(data=miage$microbiomeage.bangladesh$all,aes(x=age.sample, y=age.predicted, colour=health_analysis_groups))

microbiomeage<-function(l6.relabundtab){
  #require(caret); require(randomForest); require(ggplot2)
  bal6<-read.delim(system.file("extdata/QIIME_outputs/Bangladesh/tax_mapping7", "Subramanian_et_al_mapping_file_L6.txt", package = "metamicrobiomeR", mustWork = TRUE))
  colnames(bal6)<-tolower(colnames(bal6))
  gba<-colnames(bal6)[grep("g__",colnames(bal6))]
  #get shared genera list of Bangladesh and other studies
  glist<-list()
  gshare<-gba
  for (i in 1: length(l6.relabundtab)){
    glist[[i]]<-colnames(l6.relabundtab[[i]])[grep("g__",colnames(l6.relabundtab[[i]]))]
    gshare<-gshare[gshare %in% glist[[i]]]
  }
  bal6.share<-bal6[,c("x.sampleid","personid","familyid","ena.libraryname","health_analysis_groups","age_in_months","sex",gshare)]
  train<-as.character(bal6.share$x.sampleid[bal6.share$ena.libraryname=="BANG_HLTHY" &bal6.share$health_analysis_groups=="Healthy Singletons"])
  test<-as.character(bal6.share$x.sampleid[bal6.share$ena.libraryname!="BANG_HLTHY" & (bal6.share$health_analysis_groups=="Healthy Singletons" | bal6.share$health_analysis_groups=="Healthy Twins Triplets") |bal6.share$health_analysis_groups=="Severe Acute Malnutrition Study"])
  age.sample<-bal6.share$age_in_months
  names(age.sample)<-bal6.share$x.sampleid
  SdataMatrix <- cbind(age.sample, bal6.share[,gshare])
  traindat<-SdataMatrix[train,]
  testdat<-SdataMatrix[test,]
  #randomForest
  set.seed(123)
  rffit<- caret::train(age.sample ~ ., data = traindat, method = "rf",preProc = "center", proximity = TRUE)
  #predict on the Bangladesh data
  testage <- predict(rffit, newdata = testdat)
  testdat1<-cbind(sampleid=rownames(testdat),age.sample=testdat[,"age.sample"],age.predicted=testage)
  testdat2<-merge(testdat1,bal6.share, by.x="sampleid",by.y="x.sampleid")
  testdat2[,c("age.sample","age.predicted")]<-lapply(testdat2[,c("age.sample","age.predicted")],as.character)
  testdat2[,c("age.sample","age.predicted")]<-lapply(testdat2[,c("age.sample","age.predicted")],as.numeric)
  # list of share genera and relative importance
  taxim<-as.data.frame(randomForest::importance(rffit$finalModel))
  taxim$genera<-rownames(taxim)
  taxim$importance<-taxim[,"IncNodePurity"]
  taxim<-taxim[order(taxim[,"IncNodePurity"],decreasing = TRUE),]
  taxim<-taxim[,c("genera","importance")]
  rownames(taxim)<-NULL
  # performance on Bangladesh data
  traindat$age.predicted <- predict(rffit, newdata = traindat)
  actual<-traindat$age.sample
  predict<-traindat$age.predicted
  R2 <- 1 - (sum((actual-predict )^2)/sum((actual-mean(actual))^2))
  R2<-round(R2,2)
  ptrain<-ggplot2::ggplot() +geom_point(data=traindat,aes(x=age.sample, y=age.predicted))+
    theme(legend.text = element_text(colour="black", size = 10))+
    annotate("text", x=15, y=5,label=paste("R squared =",R2,sep=" ")) +
    labs(title="Training set")+
    theme(legend.key.size = unit(0.5, "cm"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())+
    xlab("Chronological age (month)") +ylab("Microbiome age (month)")
  bahealthy<-testdat2[testdat2$health_analysis_groups %in% c("Healthy Singletons","Healthy Twins Triplets"),]
  actual<-bahealthy$age.sample
  predict<-bahealthy$age.predicted
  R2 <- 1 - (sum((actual-predict )^2)/sum((actual-mean(actual))^2))
  R2<-round(R2,2)
  ptest<-ggplot2::ggplot() +geom_point(data=bahealthy,aes(x=age.sample, y=age.predicted))+
    theme(legend.text = element_text(colour="black", size = 10))+
    annotate("text", x=15, y=5,label=paste("R squared =",R2,sep=" ")) +
    labs(title="Test set")+
    theme(legend.key.size = unit(0.5, "cm"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())+
    xlab("Chronological age (month)") +ylab("Microbiome age (month)")
  # predict on other data
  datshare<-list()
  predictdat<-list()
  for (i in 1: length(l6.relabundtab)){
    datshare[[i]]<-l6.relabundtab[[i]][,c("age.sample",gshare)]
    rownames(datshare[[i]])<-l6.relabundtab[[i]]$x.sampleid
    age.predicted <- predict(rffit, newdata = datshare[[i]])
    predictdat[[i]]<-merge(cbind(x.sampleid=rownames(datshare[[i]]),age.predicted=age.predicted),l6.relabundtab[[i]], by="x.sampleid")
    predictdat[[i]][,c("age.sample","age.predicted")]<-lapply(predictdat[[i]][,c("age.sample","age.predicted")],as.character)
    predictdat[[i]][,c("age.sample","age.predicted")]<-lapply(predictdat[[i]][,c("age.sample","age.predicted")],as.numeric)
  }
  names(datshare)<-names(predictdat)<-names(l6.relabundtab)
  return(list(traindat.bangledesh=traindat,testdat.bangladesh=testdat,datshare=list(bangladesh=bal6.share,otherstudy=datshare),randomforestfit=rffit,sharedgenera.importance=taxim, performanceplot=list(ptrain=ptrain,ptest=ptest),microbiomeage.bangladesh=list(all=testdat2,healthy=bahealthy),microbiomeage.otherstudy=predictdat))
}
