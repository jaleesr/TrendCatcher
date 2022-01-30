#' Run TrendCatcher Main Algorithm
#'
#' This is the main function to run TrendCatcher to identify Dynamic Differentially Expressed Genes (DDEGs).
#' This function loads a rounded count matrix CSV file after the normalization and batch correction, run the core
#' algorithmand output a list object contains all the genes dynamic information.
#'
#' @param count.table.path, string contain the absolute path of the CSV file count table,
#' with first column as GENE SYMBOL or GENE ENSEMBL and first row as SAMPLE NAME (with format composed by project name,time and
#' replicateID, such as "Lung_0_Rep1")
#' @param baseline.t, one numeric variable, the baseline time of the longitudinal study. By default it is 0.

#' @param time.unit, one character variable, the time unit of longitudinal study. If choose hour,
#' please transform all sample collecting time into hour.
#' @param min.low.count, one numeric variable, the minimal count threshold for filtering low count within each time group.
#' By default it is 1.
#' @param para.core.n, one numeric variable, number of cores will be used for running TrendCatcher.
#' By default it is NA, which will use N-1 cores from computer.
#' @param dyn.p.thres, one numeric variable, the threshold of p-value of the dynamic gene. By default 0.05.
#' @param show.verbose, logic variable. If gssanova fitting failed, users can set this to TRUE, it will print out which gene failed the fitting. 
#' This process takes only one CPU, so it may be slower than the multi-core version. 
#' Normally the fitting failure is caused by low count genes. Users can manually remove it from your count table. By default set to FALSE.
#' 
#'
#' @return A list object, including "time.unit", "baseline.t", "t.arr", "Project.name", "raw.df",
#' "fitted.count" and "master.table".
#'
#' @examples
#' example.file.path<-system.file("extdata", "Brain_DemoCountTable.csv", package = "TrendCatcher")
#' \dontrun{
#' master.list<-run_TrendCatcher(count.table.path = example.file.path,
#' baseline.t = 0,
#' time.unit = "h",
#' min.low.count = 1,
#' para.core.n = NA,
#' dyn.p.thres = 0.05,
#' show.verbose = FALSE)
#' }
#' @export
#'
run_TrendCatcher<-function(count.table.path = "~/Documents/TrendCatcher/inst/extdata/Lung_DemoCountTable.csv",
                              baseline.t = 0,
                              time.unit = "h",
                              min.low.count = 1,
                              para.core.n = NA,
                              dyn.p.thres = 0.05,
                              show.verbose = F){
  # Version 1.0.0
  start_time <- Sys.time()

  ######### Create the master.list ######
  master.list<-list()
  master.list[["time.unit"]]<-time.unit

  ##### Step 1, check count table format, and filter out low count gene
  message("Read count table.")
  raw.list<-Check_CountTable_Format(count.table.path = count.table.path, min.low.count = 1)
  raw.df<-raw.list$count.table
  remove.genes<-raw.list$remove.genes
  origin.count.table<-raw.list$raw.df

  message("Count table format correct, finished loading.")
  
  t.arr<- unique(get_time_array(raw.count.df = raw.df))
  if(min(t.arr)!=baseline.t){stop("baseline.t is not your smallest time point!!! Please change it into your smallest time point.")}
  master.list[["t.arr"]]<-t.arr
  master.list[["baseline.t"]]<-baseline.t
  master.list[["Project.name"]]<-unique(as.data.frame(str_split(colnames(raw.df), "_", simplify = T))[,1])
  master.list[["raw.df"]]<-origin.count.table

  ##### Step 2, convert column data into colData
  colData<- data.frame(Sample = colnames(raw.df), Time = get_time_array(raw.df))
  colData$Time<-as.factor(colData$Time)

  ##### Step 3, estimate gene-wise dispersion
  dds<-DESeqDataSetFromMatrix(countData = raw.df, colData = colData, design = ~Time)
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds, fitType = 'local')
  gene.dispersion<-1/dispersions(dds)

  ##### Step 4. Loop through each gene to generate master table
  message("Run TrendCatcher.")
  gene.spline.list.table<-run_TrendCatcherAllGenes(raw.df = raw.df, gene.dispersion = gene.dispersion,
                                                   baseline.t = baseline.t, cores = para.core.n, show.verbose = show.verbose)

  ##### Step 5, Calculate meta p for each gene
  message("Calculate dynamic p-value for each gene.")

  dynamic.gene.df<-plyr::ddply(gene.spline.list.table, .(Gene), function(df){
    p.combine<-combine_p_single_gene(p.arr = df$t.p.val)
    p.adj<-p.combine$p.adj
    if(p.adj==0){p.adj<-1e-10}
    data.frame(dyn.p.val = formatC(p.adj, format = "e", digits = 2))
  })
  dynamic.gene.df$dyn.p.val.adj<-p.adjust(p = dynamic.gene.df$dyn.p.val, method = "BH")
  dynamic.gene.df<-dynamic.gene.df %>% dplyr::arrange(dyn.p.val.adj)
  gene.spline.list.table.merge<-merge(gene.spline.list.table, dynamic.gene.df, by = "Gene")
  rm(gene.spline.list.table)

  ##### Attach low count gene
  if(length(remove.genes)>0){
    attach.fitted.count<-data.frame(Gene = rep(remove.genes, each = length(master.list$t.arr)),
                                    Time = rep(master.list$t.arr, length(remove.genes)),
                                    Fit.Count = 0,
                                    mu = 0,
                                    disp = 0,
                                    t.p.val = 1,
                                    dyn.p.val = 1,
                                    dyn.p.val.adj = 1)
    gene.spline.list.table.all<-rbind(gene.spline.list.table.merge, attach.fitted.count)
    attach.dynamic.gene.df<-data.frame(Gene = remove.genes, dyn.p.val = 1, dyn.p.val.adj=1)
    dynamic.gene.df.all<-rbind(dynamic.gene.df, attach.dynamic.gene.df)
    rm(dynamic.gene.df)
    rm(gene.spline.list.table.merge)
  } else{
    gene.spline.list.table.all<-gene.spline.list.table.merge
    dynamic.gene.df.all<-dynamic.gene.df
  }

  master.list[["fitted.count"]]<-gene.spline.list.table.all
  master.list[["dynamic.gene.df"]]<-dynamic.gene.df.all

  ##### Step 6, Assign Dynamic Pattern
  message("Assign trajectory dynamic pattern to each gene.")

  gene.traj.pattern.df<-get_GeneTrajPattern(master.list = master.list, dyn.p.thres = dyn.p.thres, time.unit = time.unit)
  master.list[["master.table"]]<-merge(master.list$dynamic.gene.df, gene.traj.pattern.df, by = "Gene") %>% dplyr::arrange(p.adj)
  master.list$master.table<-master.list$master.table[, c("Gene", "pattern", "start.idx", "end.idx", "dynTime", "dynSign", "start.t", "end.t", "pattern_str", "dyn.p.val", "dyn.p.val.adj")]
  master.list$dynamic.gene.df<-NULL


  ##### Calculate running time
  end_time <- Sys.time()
  run.time<-difftime(end_time, start_time, units='mins')
  message("Finished.")
  cat("Running time is", round(run.time,2), "mins \n")

  ##### Return master.list
  return(master.list)
}
