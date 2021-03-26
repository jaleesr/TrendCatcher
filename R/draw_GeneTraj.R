#' Draw gene(s) trajectory with observed data and fitted data
#'
#' This function takes the master.list object output from run_TrendCatcher function, and an array of gene(s). 
#' It will draw gene(s) trajectory with observed data and fitted data.
#'  
#' @param master.list.new, a list object. Output from the run_TrendCatcher with ID conversion to add Symbol column to master table.
#' @param gene.symbol.arr, a character array. It must be a subset of row names from the master.list$master.table$Symbol. 
#' The Symbol column need get_GeneEnsembl2Symbol function to convert original ensembl ID into gene symbol.
#' @param savepdf.path, an obsolute file path to save the figure as PDF file. By default is NA, it will be printed.
#' @param ncol, an integer variable. If more than one gene need to be plotted, it will layout as grid structure. 
#' This represents the number of column of the grid layout.
#' @param nrow, an integer variable. If more than one gene need to be plotted, it will layout as grid structure. 
#' This represents the number of row of the grid layout.
#' @param fig.width, a numeric variable. If save figure as PDF file, the width of the PDF file. By default is 15.
#' @param fig.height, a numeric variable. If save figure as PDF file, the height of the PDF file. By default is 10.
#' 
#' 
#' 
#' @param count.table, the count table of mRNA.
#' @param gene.name, name of the gene you would like to plot.
#' @param master.list, the list object returned from run_TrendCatcher.
#' @return "arrangelist" "list" object.
#' 
# 
#' @examples
#' gene.symbol.arr <-  c("Cxcl1", "Cxcl5", "Cxcl10", "Cxcl11", "Abcb1b", "Icam1", "Ifitm1", "Ifitm2", "Ifitm3")
#' example.file.path<-system.file("extdata", "BrainMasterList.rda", package = "TrendCatcher")
#' load(file= example.file.path)
#' gene.symbol.df<-get_GeneEnsembl2Symbol(ensemble.arr = master.list$master.table$Gene)
#' master.table.new<-cbind(master.list$master.table, gene.symbol.df[match(master.list$master.table$Gene, gene.symbol.df$Gene), c("Symbol", "description")])
#' master.list$master.table<-master.table.new
#' gplots<-draw_GeneTraj(master.list, gene.symbol.arr, savepdf.path = NA, ncol = 3, nrow = 3, fig.width =15, fig.height=10)
#' @export
#'
draw_GeneTraj<-function(master.list, gene.symbol.arr, savepdf.path=NA, ncol = 5, nrow = 3, fig.width = 15, fig.height = 10){
  
  # Subset master.table
  sub<-master.list$master.table %>% filter(Symbol %in% gene.symbol.arr)
  count.table<-master.list$raw.df[rownames(master.list$raw.df) %in% sub$Gene,]


  # plot list
  p<-list()
  for(i in 1:length(gene.symbol.arr)){

    gene.symbol<-gene.symbol.arr[i]
    gene.dyn.p<-sub$dyn.p.val.adj[match(gene.symbol, sub$Symbol)]
    gene.pattern<-sub$pattern_str[match(gene.symbol, sub$Symbol)]
    gene.ensembl<-sub$Gene[match(gene.symbol, sub$Symbol)]

    figure.title<-paste0(gene.symbol, "\n","dyn.p.adj=", signif(gene.dyn.p,digits = 3), "\n", gene.pattern)

    ### Draw plot, find count data and fit data
    idx<-which(rownames(count.table) == gene.ensembl)
    gene.row.info<-as.numeric(count.table[idx,])
    t.arr<-get_time_array(raw.count.df = count.table)
    df<-data.frame(Time = t.arr, Count = gene.row.info)
    p[[i]]<-ggplot()
    p[[i]]<-p[[i]]+geom_point(data = df, aes(x = Time, y = Count))
    fit.df<-master.list$fitted.count %>% filter(Gene == gene.ensembl)
    disp<-fit.df$disp[1]
    mu<-fit.df$mu[1]
    p[[i]]<-p[[i]]+geom_point(data = fit.df, aes(x = Time, y = Fit.Count), colour = "red")
    p[[i]]<-p[[i]]+geom_line(data = fit.df, aes(x=Time, y = Fit.Count), colour = "red")
    p[[i]]<-p[[i]]+geom_hline(yintercept = qnbinom(0.05, disp,mu =  mu), colour = "grey")
    p[[i]]<-p[[i]]+geom_hline(yintercept = qnbinom(0.95, disp,mu =  mu), colour = "grey")
    p[[i]]<-p[[i]]+ggtitle(figure.title)
  }
  ml <- marrangeGrob(p, ncol=ncol, nrow =nrow)
  if(!is.na(savepdf.path)){
    ggsave(savepdf.path, ml,
           width = fig.width, height = fig.height)
    return(p)
  }else{
    return(marrangeGrob(grobs=p, nrow=nrow, ncol=ncol))
  }
}
