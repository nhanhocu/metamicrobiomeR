
#' Read multiple files
#'
#' This function reads multiple files of the same pattern in a directory into R.
#' @param patht path to files.
#' @param patternt pattern of files. Options are ".txt" or ".csv".
#' @param assignt assign files. Default is "no".
#' @param study name of the study. Default is NULL.
#' @param tolower.colnames turn all column names to lower case. Default is TRUE.
#' @keywords multiple files
#' @export
#' @return list of all data files in the path
#' @examples
#' \donttest{
#' # The data used for this example are available
#' # in the "metamicrobiomeR" package version in Github.
#' # Download example data from the package github repo
#' #setwd("your directory") #put your working directory inside the quotation marks
#' download.file(url = "https://github.com/nhanhocu/metamicrobiomeR/archive/master.zip",
#' destfile = "metamicrobiomeR-master.zip")
#' # unzip the .zip file
#' unzip(zipfile = "metamicrobiomeR-master.zip")
#' patht<-paste(getwd(),
#' "metamicrobiomeR-master/inst/extdata/QIIME_outputs/Bangladesh/alpha_div_collated", sep="/")
#' alpha.ba<-read.multi(patht=patht,patternt=".txt", assignt="no",study="Bangladesh")
#' }

read.multi<-function (patht,patternt=".txt",assignt="no",study=NULL,tolower.colnames=TRUE){
  filenames <- paste(patht,list.files(path=patht, pattern=patternt),sep="/")
  if (patternt==".txt"){
    tmp <- lapply(filenames, function(x) utils::read.delim(file=x))
  }
  if (patternt==".csv"){
    tmp <- lapply(filenames, function(x) utils::read.csv(file=x))
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

