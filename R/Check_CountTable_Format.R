#' Check the Input Count Table Format
#' 
#' This function takes the CSV file path of the count table provided by user. The input count table must
#' be a CSV file, includes integer count table, first column are the GENE SYMBOL or GENE ENSEMBL, and first row are SAMPLE NAME.
#' The count table is a rounded up after normalization and batch correction.
#' Please check TrendCatcher QC functions for normalization and batch correction details. The row names must be in the format of
#' "ProjectName_Time_Rep1" format. ProjectName is a string, can be single letter. Time is a integer. Rep is the replicate ID.
#' This function checked the count table format, and order the sample columns based on its time order.
#' It will return a right formatted count table to run TrendCatcher.
#'
#'
#' @param count.table.path, string contain the absolute path of the CSV file count table, 
#' with first column as GENE SYMBOL or GENE ENSEMBL and first row as SAMPLE NAME (with format composed by project name,time and 
#' replicateID, such as "Lung_0_Rep1")
#' @param min.low.count, one numeric variable, the minimal count threshold for filtering low count within each time group. 
#' By default it is 1.
#' @return A list object, contains "raw.df", original count table ordered by time and replicate ID; 
#' "count.table", filtered out low genes for further fitting; "removed.genes", low count genes.
#'
#' @examples
#' example.file.path<-system.file("extdata", "Lung_DemoCountTable.csv", package = "TrendCatcher")
#' \dontrun{
#' count<-Check_CountTable_Format(example.file.path, min.low.count = 1)
#' }
#' @export
#'
#'
Check_CountTable_Format<-function(count.table.path, min.low.count = 1){
  ### 1.Check the count table path
  if(!file.exists(count.table.path)) stop("Count table file doesn't exist!")

  ### 2. Read in count table csv file
  raw.df<-read.csv(count.table.path, row.names = 1)

  ### 3. Check column name format is in the format of "ProjectName_Time_Rep1" and Reorder it by time
  pass.check.n<-sum(grepl("[[:alpha:]]_[[:digit:]]*_Rep[[:digit:]]", colnames(raw.df)))
  if(pass.check.n!=ncol(raw.df)){
    stop("The column name must be in the format of Project_Time_Rep1!")
  }

  ### 4. Check row name are string
  row.name.string<-is.character(rownames(raw.df)[1]) & length(rownames(raw.df)[1])
  if(row.name.string!=1){
    stop("The row name must be gene name!")
  }

  ### 5. Count table must be non-negative integer
  is.interger.table<-floor(raw.df)==raw.df
  if(!all(is.interger.table)){
    stop("The count table must be an integer table!")
  }
  if(any(is.na(raw.df))){
    stop("There is NA value in the count table!")
  }
  is.positive.table<-raw.df>=0
  if(!all(is.positive.table)){
    stop("The count table must be non-negative table!")
  }
  
  ### 6. Must be 1 single project
  prj.arr<-unique(as.data.frame(str_split(colnames(raw.df), "_", simplify = T))[,1])
  if(length(prj.arr)>1){
    stop("The project name must be the same!!!")
  }
  

  ##################### If passed all the check, now order the count table by time ##########

  ### 7. Get time array & Print
  time.arr<-get_time_array(raw.count.df = raw.df)
  t.unique.arr<-unique(time.arr)
  t.unique.arr<-t.unique.arr[order(as.numeric(t.unique.arr))]
  
  #### The baseline time point requires more than one sample !!!
  t.min<-min(t.unique.arr)
  baseline.idx<-grep(paste0(prj.arr,"_",t.min, "_"), colnames(raw.df))
  if(length(baseline.idx)==1){stop("TrendCatcher requires more than 1 replicate at baseline time point!!!")}

  ### 8. Get rep array
  rep.arr<-get_rep_array(raw.count.df = raw.df)

  ### 9. Arrange the samples by time order and replicate order
  original.order<-colnames(raw.df)
  colnames.ordered<-NULL
  for(t in t.unique.arr){
    sample.t.index<-grep(paste0(prj.arr,"_",t, "_"), colnames(raw.df))
    tmp.df<-as.data.frame(raw.df[,sample.t.index])
    colnames(tmp.df)<-colnames(raw.df)[sample.t.index]
    t.rep.arr<-get_rep_array(raw.count.df = tmp.df)
    t.rep.arr.order<-t.rep.arr[order(as.numeric(t.rep.arr))]
    colnames.ordered<-c(colnames.ordered, paste0(prj.arr, "_", t, "_", "Rep",t.rep.arr.order))
  }
  raw.df<-raw.df[colnames.ordered]

  ### 10. Filter out low count gene from each time group for later model fitting
  ######## Low count genes in each time group will cause model fitting fail
  t.min.list<-list()
  for(i in 1:length(t.unique.arr)){
    t<-t.unique.arr[i]
    sample.t.index<-grep(paste0(prj.arr,"_",t, "_"), colnames(raw.df))
    tmp.df<-as.data.frame(raw.df[,sample.t.index])
    colnames(tmp.df)<-colnames(raw.df)[sample.t.index]
    t.min.list[[i]]<-apply(tmp.df, MARGIN = 1, min)
  }
  t.min.df<-do.call(cbind, t.min.list)
  remove.idx<-which(rowSums(t.min.df)<=min.low.count)
  if(length(remove.idx)!=0){
    count.table<-raw.df[-remove.idx,]
  } else{count.table<-raw.df}
  message(c("Passed Format Check. Count table is ", paste0(ncol(raw.df), " Samples, ", nrow(raw.df), " Genes")))
  message(c("Found ", length(remove.idx), " low count genes."))
  remove.genes<-names(remove.idx)
  raw.list<-list(raw.df = raw.df, count.table = count.table, remove.genes = remove.genes)
  return(raw.list)
}
