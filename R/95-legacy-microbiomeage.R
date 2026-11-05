#' Predict microbiome age.
#'
#' This function predicts microbiome age using Random Forest model based on relative abundances of bacterial genera shared with the Bangladesh study (Subramanian et al 2014).
#' This function gets the shared genera list between the Bangladesh study and all other included studies,
#' get the training and test sets from Bangladesh data based on the shared genera list,
#' fit the train Random Forest model and predict microbiome age in the test set of Bangladesh data and data from all included studies,
#' check for performance of the model based on the shared genera list on Bangladesh healthy cohort data, reproduce the findings of the Bangladesh malnutrition study.
#' @param l6.relabundtab list of taxa summary table from phylum up to genus level merged to mapping file outputed from QIIME of all included studies.
#' @param bal6 reference data for model training (taxa summary table from phylum up to genus level merged to mapping file outputed from QIIME of the Bangladesh study).
#' @param seed optional random seed, so that the Random Forest fit can be repeated. Default is NULL, as in earlier versions.
#' @return list of training and test sets of Bangladesh data, shared genera list, relative abundance data of shared genera, randomforest fit, RF model performance plot,predicted microbiome age of Bangladesh data and data of other included studies.
#' @keywords microbiome age Random Forest.
#' @export
#' @examples
#' \donttest{
#' if (requireNamespace("caret", quietly = TRUE) && requireNamespace("randomForest", quietly = TRUE)) {
#'   # A small simulated example. Real use needs the Bangladesh reference table.
#'   set.seed(1)
#'   genera <- paste0("k__bacteria.p__firmicutes.c__c1.o__o1.f__f1.g__g", 1:6)
#'   sim <- function(n) {
#'     age <- runif(n, 0, 24)
#'     tab <- sapply(seq_along(genera), function(j) {
#'       stats::plogis(-3 + age / (6 * j) + rnorm(n, sd = 0.3))
#'     })
#'     colnames(tab) <- genera
#'     data.frame(age = age, tab, check.names = FALSE)
#'   }
#'   b <- sim(60)
#'   bal6 <- data.frame(x.sampleid = paste0("s", 1:60), personid = paste0("p", 1:60),
#'                      familyid = paste0("f", 1:60),
#'                      ena.libraryname = rep(c("BANG_HLTHY", "OTHER"), each = 30),
#'                      health_analysis_groups = "Healthy Singletons",
#'                      age_in_months = b$age, sex = "F", b[, genera], check.names = FALSE)
#'   o <- sim(20)
#'   other <- list(study1 = data.frame(x.sampleid = paste0("o", 1:20), age.sample = o$age,
#'                                     o[, genera], check.names = FALSE))
#'   miage <- microbiomeage(l6.relabundtab = other, bal6 = bal6, seed = 1)
#'   miage$sharedgenera.importance
#' }
#' }

microbiomeage<-function(l6.relabundtab, bal6, seed=NULL){
  .need("caret", "microbiomeage()")
  .need("randomForest", "microbiomeage()")
  if (!is.null(seed)) set.seed(seed)
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
  rffit<- caret::train(age.sample ~ ., data = traindat, method = "rf",preProc = "center", proximity = TRUE)
  #predict on the Bangladesh data
  testage <- stats::predict(rffit, newdata = testdat)
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
  traindat$age.predicted <- stats::predict(rffit, newdata = traindat)
  actual<-traindat$age.sample
  predict<-traindat$age.predicted
  R2 <- 1 - (sum((actual-predict )^2)/sum((actual-mean(actual))^2))
  R2<-round(R2,2)
  ptrain<-ggplot2::ggplot() +ggplot2::geom_point(data=traindat,ggplot2::aes(x=age.sample, y=age.predicted))+
    ggplot2::theme(legend.text = ggplot2::element_text(colour="black", size = 10))+
    ggplot2::annotate("text", x=15, y=5,label=paste("R squared =",R2,sep=" ")) +
    ggplot2::labs(title="Training set")+
    ggplot2::theme(legend.key.size = ggplot2::unit(0.5, "cm"),
          axis.line = ggplot2::element_line(colour = "black"),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.border = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank())+
    ggplot2::xlab("Chronological age (month)") +ggplot2::ylab("Microbiota age (month)")
  bahealthy<-testdat2[testdat2$health_analysis_groups %in% c("Healthy Singletons","Healthy Twins Triplets"),]
  actual<-bahealthy$age.sample
  predict<-bahealthy$age.predicted
  R2 <- 1 - (sum((actual-predict )^2)/sum((actual-mean(actual))^2))
  R2<-round(R2,2)
  ptest<-ggplot2::ggplot() +ggplot2::geom_point(data=bahealthy,ggplot2::aes(x=age.sample, y=age.predicted))+
    ggplot2::theme(legend.text = ggplot2::element_text(colour="black", size = 10))+
    ggplot2::annotate("text", x=15, y=5,label=paste("R squared =",R2,sep=" ")) +
    ggplot2::labs(title="Test set")+
    ggplot2::theme(legend.key.size = ggplot2::unit(0.5, "cm"),
          axis.line = ggplot2::element_line(colour = "black"),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.border = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank())+
    ggplot2::xlab("Chronological age (month)") +ggplot2::ylab("Microbiome age (month)")
  # predict on other data
  datshare<-list()
  predictdat<-list()
  for (i in 1: length(l6.relabundtab)){
    datshare[[i]]<-l6.relabundtab[[i]][,c("age.sample",gshare)]
    rownames(datshare[[i]])<-l6.relabundtab[[i]]$x.sampleid
    age.predicted <- stats::predict(rffit, newdata = datshare[[i]])
    predictdat[[i]]<-merge(cbind(x.sampleid=rownames(datshare[[i]]),age.predicted=age.predicted),l6.relabundtab[[i]], by="x.sampleid")
    predictdat[[i]][,c("age.sample","age.predicted")]<-lapply(predictdat[[i]][,c("age.sample","age.predicted")],as.character)
    predictdat[[i]][,c("age.sample","age.predicted")]<-lapply(predictdat[[i]][,c("age.sample","age.predicted")],as.numeric)
  }
  names(datshare)<-names(predictdat)<-names(l6.relabundtab)
  return(list(traindat.bangledesh=traindat,testdat.bangladesh=testdat,datshare=list(bangladesh=bal6.share,otherstudy=datshare),randomforestfit=rffit,sharedgenera.importance=taxim, performanceplot=list(ptrain=ptrain,ptest=ptest),microbiomeage.bangladesh=list(all=testdat2,healthy=bahealthy),microbiomeage.otherstudy=predictdat))
}
