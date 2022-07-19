#' Draw GOHeatmap containing Terms from TimeHeatmap and included Genes
#'
#' This funcitons takes the master.list output from run_TrendCatcher, and merge.df output from draw_TimeHeatmap_GO and draw_TimeHeatmap_enrichR.
#' And showing all the genes used for enrichment analysis and their logFC compared to previous break point.
#'
#' @param master.list, a list object. The output from run_TrendCatcher function, contains master.table element.
#' @param time.window, a character. Must be one of the merge.df$t.name.
#' @param go.terms, a character array. Must be an array of go terms from the merge.df$Description.
#' @param merge.df, a dataframe. The output dataframe from output list of draw_TimeHeatmap_GO or draw_TimeHeatmap_enrichR. Use $merge.df to obtain it.
#' @param logFC.thres, a numeric variable. The logFC threshold compared to each genes previous break point expression level.
#' By default is 2, meaning for each gene, the current time window's expression level is 2-fold compared to previous
#' break point's expression level.
#' @param figure.title, character
#' @param save.tiff.path, by default is NA
#' @param tiff.res, resolution
#' @param tiff.width, figure width
#' @param tiff.height, figure height
#'
#' @return A list object, including elements named GOheatmap and GOheatmapDat. GOheatmap is a ComplexHeatmap object for figure. GOheatmapDat is a data.frame include log2FC 
#' value of each gene's expression change compared to the previous break point. 
#'
#' @examples
#' \dontrun{
#' example.file.path<-system.file("extdata", "BrainMasterList.rda", package = "TrendCatcher")
#' load(example.file.path)
#' time_heatmap<-draw_TimeHeatmap_GO(master.list = master.list)
#' merge.df<-time_heatmap$merg.df
#' time.window<-"0h-6h"
#' go.terms<-c("regulation of defense response", "leukocyte migration", "myeloid leukocyte migration", "leukocyte chemotaxis",
#' "granulocyte chemotaxis", "cellular response to chemokine", "chemokine-mediated signaling pathway", "angiogenesis", "sprouting angiogenesis",  "respone to bacterium", "leukocyte mediated immunity")
#' go.df<-draw_GOHeatmap(master.list = master.list, time.window = "0h-6h", go.terms = go.terms, merge.df = merge.df, logFC.thres = 2)
#' }
#' @export
#'
#'

draw_GOHeatmap<-function(master.list, time.window="", go.terms="",
                         merge.df=NA, logFC.thres=2, figure.title = "", save.tiff.path = NA, tiff.res = 100, tiff.width = 1500, tiff.height =1500){
  
  ####### Check If there is the Symbol column exist!!!! ######
  idx<-grep("Symbol", colnames(master.list$master.table))
  if(length(idx)==0){stop("Please add Symbol column to your master.list$master.table. By default, it should be a column of gene SYMBOLs!")}
  
  
  if(!is.data.frame(merge.df)){stop("Merge.df missing!!!")}
  if(!time.window %in% merge.df$t.name){stop("Time window format is wrong!!!")}
  if(length(go.terms)==0){stop("Please enter multiple go.terms!!!")}
  
  #### Add symbol to fitted.count
  master.list$fitted.count$Symbol<-master.list$master.table$Symbol[match(master.list$fitted.count$Gene, master.list$master.table$Gene)]

  # get each term gene
  sub.merge.df<-merge.df %>% dplyr::filter(Description %in% go.terms & t.name == time.window)
  end.t.str<-str_split(time.window, "-", simplify = T)[2]
  end.t<-as.numeric(substring(end.t.str, first = 1, last = nchar(end.t.str)-1))

  term.gene.df<-ddply(sub.merge.df, .(type, Description), function(df){
    id.arr<-as.character(str_split(df$geneID, "/", simplify = T))
    logFC.arr<-NULL
    prev.bk.t.arr<-NULL
    # For each gene, go find the prev. bk's value
    for(i in 1:length(id.arr)){
      gene.name<-id.arr[i]
      gene.dyn.info<-master.list$master.table %>% dplyr::filter(Symbol == gene.name)
      dyn.t.arr<-str_split(master.list$master.table$dynTime, "_", simplify = T)
      dyn.t.arr<-as.numeric(dyn.t.arr[-length(dyn.t.arr)])
      idx<-which(dyn.t.arr < end.t)
      if(length(idx)==0){prev.bk.t<-0}else{prev.bk.t<-dyn.t.arr[max(idx)]}
      prev.bk.t.arr<-c(prev.bk.t.arr, prev.bk.t)
      # Get the count value for end.t and prev.t
      gene.count.df<-master.list$fitted.count %>% dplyr::filter(Symbol == gene.name)
      prev.bk.count<-gene.count.df$Fit.Count[which(gene.count.df$Time == prev.bk.t)]
      end.t.count<-gene.count.df$Fit.Count[which(gene.count.df$Time == end.t)]
      logFC <- log(end.t.count, base = 2) - log(prev.bk.count, base = 2)
      logFC.arr<-c(logFC.arr, logFC)
    }
    data.frame(gene = id.arr, logFC.prev.bk = logFC.arr, prev.bk.t = prev.bk.t)
  })

    term.gene.df.symbol <- cbind(term.gene.df, Symbol = term.gene.df$gene)
    term.gene.df.symbol<-term.gene.df.symbol %>% dplyr::filter(abs(logFC.prev.bk) > logFC.thres)



  # Construct the heatmap matrix
  n.col<-length(unique(term.gene.df.symbol$Description))
  n.row<-length(unique(term.gene.df.symbol$Symbol))
  mat<-matrix(data = 0, nrow = n.row, ncol = n.col)
  colnames(mat)<-unique(term.gene.df.symbol$Description)
  rownames(mat)<-unique(term.gene.df.symbol$Symbol)
  for(i in 1:n.col){
    term.name<-colnames(mat)[i]
    for(j in 1:n.row){
      gene.name<-rownames(mat)[j]
      sub<-term.gene.df.symbol %>% dplyr::filter(Description == term.name & Symbol == gene.name)
      if(nrow(sub)>0){
        mat[j,i]<-sub$logFC.prev.bk
      }
    }
  }

  # Define paint color
  col_fun = colorRamp2(c(-10, 0, 10), c("blue", "white", "red"))
  col.split = data.frame(cutree(hclust(dist(mat)), k = 3))
  row.split = data.frame(cutree(hclust(dist(t(mat))), k = 3))
  colnames(mat)<-str_wrap(colnames(mat), width = 50)


  ht<-Heatmap(mat,name = "logFC vs. prevBK", cluster_row_slices = T,
              cluster_column_slices = T,
              cluster_rows = T,
              cluster_columns = T,
              row_split = col.split,
              column_split = row.split,
              col = col_fun,
              show_column_dend = F,
              show_row_dend = F,
              width = 5, height = 3,
              rect_gp = gpar(col = "white", lwd = 1),
              row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), border = T,
              row_km = 1, column_km = 1, row_title = NULL, column_title = time.window,
              show_heatmap_legend = T, column_names_rot = 45, column_names_side = "top", row_names_side = "left",
              row_names_gp = gpar(fontsize = 10))

  if(is.na(save.tiff.path)){
    print(ht)
  } else{
    tiff(filename = save.tiff.path, res = tiff.res, width = tiff.width, height = tiff.height)
    print(ht)
    dev.off()
  }
  return(list(GOheatmap = ht, GOheatmapDat = term.gene.df.symbol))
}

