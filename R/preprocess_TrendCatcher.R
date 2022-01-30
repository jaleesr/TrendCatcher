#' Preprocessing for TrendCatcher
#'
#' This is the preprocessing function to prepare the input count table for run_TrendCatcher function.
#' It takes the CSV file count table and logic variables to check if normalization and batch correction needed for preprocessing.
#' It creates pdf figure report to show before and after of normalization and batch correction to assess quality control (QC).
#'
#' @param count.table.path, string contain the absolute path of the CSV file count table,
#' with first column as GENE SYMBOL or GENE ENSEMBL and first row as SAMPLE NAME (with format composed by project name,time and
#' replicateID, such as "Lung_0_Rep1")
#' @param need.batch.correction, logic variable. If batch correction is needed. By default is TRUE.
#' @param need.normalization, logic variable. If normalization is needed. By default is TRUE.
#' @param batch.arr, a numeric vector of batch number. Need to be the same length as the number of samples.
#' @param pdf.file.path, an absolute file path for save the QC report file. If not need, set if to NA. The report will be printed.
#' By default is NA.
#' @param pdf.width, a numeric variable. The PDF file width size. By default is 8.
#' @param pdf.height, a numeric variable. The PDF file height size. By default is 10.
#' @param n.low.count, a numeric variable. The minimal number to filter low count genes. By default is 10.
#'
#'
#' @return A matrix array object.
#'
#' @examples
#' example.file.path<-system.file("extdata", "Brain_DemoCountRawTable.csv", package = "TrendCatcher")
#' \dontrun{
#' count.table<-preprocess_TrendCatcher(count.table.path = example.file.path,
#' need.batch.correction = TRUE,
#' need.normalization = TRUE,
#' batch.arr ="",
#' pdf.file.path = NA,
#' pdf.width=8, pdf.height=10,
#' n.low.count = 10)
#' }
#' @export
#'

preprocess_TrendCatcher<-function(count.table.path = "",
                                  need.batch.correction = TRUE,
                                  need.normalization = TRUE,
                                  batch.arr ="",
                                  pdf.file.path = NA,
                                  pdf.width=8, pdf.height=10,
                                  n.low.count = 10){

  ### 1.Check the count table path
  if(!file.exists(count.table.path)) stop("Count table file doesn't exist!")

  ### 2. Read in count table csv file
  raw.count<-read.csv(count.table.path, row.names = 1)

  ### 3. Check column name format is in the format of "ProjectName_Time_Rep1" and Reorder it by time
  pass.check.n<-sum(grepl("[[:alpha:]]_[[:digit:]]*_Rep[[:digit:]]", colnames(raw.count)))
  if(pass.check.n!=ncol(raw.count)){
    stop("The column name must be in the format of Project_Time_Rep1!")
  }

  ### 4. Check row name are string
  row.name.string<-is.character(rownames(raw.count)[1]) & length(rownames(raw.count)[1])
  if(row.name.string!=1){
    stop("The row name must be gene name!")
  }

  ### 5. Count table must be non-negative integer
  is.interger.table<-floor(raw.count)==raw.count
  if(!all(is.interger.table)){
    stop("The count table must be an integer table!")
  }
  if(any(is.na(raw.count))){
    stop("There is NA value in the count table!")
  }
  is.positive.table<-raw.count>=0
  if(!all(is.positive.table)){
    stop("The count table must be non-negative table!")
  }
  

  ### 6. Must be 1 single project
  prj.arr<-unique(as.data.frame(str_split(colnames(raw.count), "_", simplify = T))[,1])
  if(length(prj.arr)>1){
    stop("The project name must be the same!!!")
  }


  ######## Batch Correction ###############
  if(need.batch.correction){
    if(length(batch.arr)!=ncol(raw.count)){stop("The batch.arr should be the same as the sample number!!!!")}
    ### Prepare pheno data
    pheno.df<-data.frame(sample = colnames(raw.count), batch = batch.arr)
  } else{
    pheno.df<-data.frame(sample = colnames(raw.count), batch = rep(1, ncol(raw.count)))
  }
  rownames(pheno.df)<-pheno.df$sample
  pheno.df$time<-paste0("T",str_split(colnames(raw.count), "_", simplify = T)[,2])
  pheno.df$time<-as.factor(pheno.df$time)
  pheno.df$batch<-as.factor(pheno.df$batch)

  ### Create DGElist object ####
  dge.origin <- DGEList(counts = raw.count, genes = rownames(raw.count))
  group<-factor(pheno.df$time)
  keep<-filterByExpr(dge.origin, group = group)
  dge.origin <- dge[keep,,keep.lib.sizes=FALSE]
  adjusted <- ComBat_seq(dge.origin$counts, batch=pheno.df$batch, group=group)

  ######## Normaliztion ###############
  if(need.normalization){
    dge.batch.corrected <- DGEList(counts = adjusted, genes = rownames(adjusted))
    dge.batch.corrected<-calcNormFactors(dge.batch.corrected)
    logCPM <- cpm(dge.batch.corrected, log=TRUE, prior.count=3)
  } else{
    dge.batch.corrected <- DGEList(counts = adjusted, genes = rownames(adjusted))
    logCPM <- cpm(dge.batch.corrected, log=TRUE, prior.count=3)
  }

  ############## Generate QC plot ################
  if(!is.na(pdf.file.path)){
    pdf(pdf.file.path, width = pdf.width, height = pdf.height)
    par(mfrow=c(3,2))
    boxplot(log2(adjusted+1), las=2, main="Before Normlization")
    boxplot(logCPM, las=2, main="After Normalization")

    # Set color
    t<-as.numeric(str_split(pheno.df$time, "T", simplify = T)[,2])
    n.time<-length(unique(t))
    col.1<-brewer.pal(9, "Set1")[1:n.time]
    col.time <- rep(col.1, table(t))

    batch <- pheno.df$batch
    n.batch <- length(unique(batch))
    col.2<-brewer.pal(8, "Set2")[1:n.batch]
    col.batch <- col.2[batch]

    plotMDS(cpm(dge.origin, log=TRUE, prior.count=3), col = col.time, pch=16, labels = colnames(dge.origin))
    title("Before batch correction (color by sample name)")
    plotMDS(cpm(dge.origin, log=TRUE, prior.count=3), col = col.batch, pch=16, labels = pheno.df$batch)
    title("Before batch correction (color by batch)")


    plotMDS(logCPM, col=col.time, pch = 16, labels = colnames(logCPM))
    title("After batch correction (color by sample name)")
    plotMDS(logCPM,col=col.batch, pch=16, labels = raw.data$pheno$batch)
    title("After batch correction (color by batch)")
    dev.off()
  } else{
    par(mfrow=c(3,2))
    boxplot(log2(adjusted+1), las=2, main="Before Normlization")
    boxplot(logCPM, las=2, main="After Normalization")

    # Set color
    t<-as.numeric(str_split(pheno.df$time, "T", simplify = T)[,2])
    n.time<-length(unique(t))
    col.1<-brewer.pal(9, "Set1")[1:n.time]
    col.time <- rep(col.1, table(t))

    batch <- pheno.df$batch
    n.batch <- length(unique(batch))
    col.2<-brewer.pal(8, "Set2")[1:n.batch]
    col.batch <- col.2[batch]

    plotMDS(cpm(dge.origin, log=TRUE, prior.count=3), col = col.time, pch=16, labels = colnames(dge.origin))
    title("Before batch correction (color by sample name)")
    plotMDS(cpm(dge.origin, log=TRUE, prior.count=3), col = col.batch, pch=16, labels = pheno.df$batch)
    title("Before batch correction (color by batch)")


    plotMDS(logCPM, col=col.time, pch = 16, labels = colnames(logCPM))
    title("After batch correction (color by sample name)")
    plotMDS(logCPM,col=col.batch, pch=16, labels = raw.data$pheno$batch)
    title("After batch correction (color by batch)")
  }
  ##### Return filtered count table
  count<-round(2^(logCPM))
  count<-count[rowSums(count)>n.low.count,]
  return(count)
}
