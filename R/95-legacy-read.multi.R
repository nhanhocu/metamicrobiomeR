
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
#' d <- tempfile(); dir.create(d)
#' utils::write.table(data.frame(Sample = "a", Shannon = 1.2), file.path(d, "shannon.txt"),
#'                    sep = "\t", row.names = FALSE, quote = FALSE)
#' utils::write.table(data.frame(Sample = "a", Chao1 = 40), file.path(d, "chao1.txt"),
#'                    sep = "\t", row.names = FALSE, quote = FALSE)
#' alpha <- read.multi(patht = d, patternt = ".txt", assignt = "no", study = "example")
#' names(alpha)

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

