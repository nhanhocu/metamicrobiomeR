
#' Read multiple files
#'
#' This function reads multiple files of the same pattern in a directory into R.
#' @param patht path to files.
#' @param patternt pattern of files. Options are ".txt" or ".csv" or ".sav" or ".sas7bdat".
#' @param assignt assign files. Default is "no".
#' @param study name of the study. Default is NULL.
#' @param tolower.colnames turn all column names to lower case. Default is TRUE.
#' @keywords multiple files
#' @export
#' @return list of all data files in the path
#' @examples
#' patht<-system.file("extdata/QIIME_outputs/Bangladesh/alpha_div_collated", package = "metamicrobiomeR", mustWork = TRUE)
#' alpha.ba<-read.multi(patht=patht,patternt=".txt",assignt="no",study="Bangladesh")

read.multi<-function (patht,patternt=".txt",assignt="no",study=NULL,tolower.colnames=TRUE){
  filenames <- paste(patht,list.files(path=patht, pattern=patternt),sep="/")
  if (patternt==".txt"){
    tmp <- lapply(filenames, function(x) read.delim(file=x))
  }
  if (patternt==".csv"){
    tmp <- lapply(filenames, function(x) read.csv(file=x))
  }
  if (patternt==".sav"){
    tmp<-lapply(filenames,function(x) read.spss(file=x,use.value.labels=FALSE, to.data.frame=TRUE))
  }
  if (patternt==".sas7bdat"){
    tmp<-lapply(filenames,function(x) read.sas7bdat(file=x, debug=FALSE))
  }
  names(tmp)<-tolower(gsub(patternt,"", filenames))
  names(tmp)<-basename(tools::file_path_sans_ext(names(tmp)))
  for (i in 1:length(names(tmp))){
    if (tolower.colnames==TRUE){
      colnames(tmp[[i]])<-tolower(colnames(tmp[[i]]))
    }
    tmp[[i]][,"study"]<-study
    if (assignt=="yes"){
      assign(names(tmp[i]),tmp[[names(tmp[i])]])
    }
  }
  return(tmp)
  rm(i)
}

