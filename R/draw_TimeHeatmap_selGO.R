#' Subset TimeHeatmap by providing a manually selected non-redundant GO terms
#' 
#' Some GO terms are redundant. Users can manually select GO terms that are shown in the GO.df element in the TimeHeatmap object and show the TimeHeatmap figure.
#' @param time_heatmap, a list, the output of draw_TimeHeatmap_GO function. A TimeHeatmap object, with GO.df element included.
#' @param sel.go, a character variable. An array of character names of GO terms, that match the GO terms from the Description column of GO.df.
#' @param master.list, a list, the output of run_TrendCatcher function, a master.list object.
#' @param GO.perc.thres, a numeric variable. A threshold to filter out GOs that only a little percentage of the genes are DDEGs. By default is 0.
#' @param nDDEG.thres, an integer variable. A threshold to filter out GOs that only a small number of genes included. By default is 0.
#' @param term.width, an integer variable. The character length for each GO term. If one GO term is super long, we can wrap
#' it into term.width of strings into multiple rows. By default is 80.
#' @param figure.title, a character variable. The main title of TimeHeatmap.
#' @param save.tiff.path, a character variable, the file path to save the TIFF figure. If set to NA, it will plot it out. By default is NA.
#' @param tiff.res, a numeric variable, the resolution of the TIFF figure. By default is 100. 
#' @param tiff.width, a numeric variable, the width of the TIFF figure. By default is 1500. 
#' @param tiff.height, a numeric variable, the height of the TIFF figure. By default is 1500. 
#'
#' @return A list object, including elements names time.heatmap, merge.df and GO.df.
#' time.heatmap is the ComplexHeatmap object. merge.df includes all the GO enrichment result and their activation/deactivation time window.
#' GO.df includes GO enrichment used for plot TimeHeatmap and all the individual genes within each time window. 
#' 
#' @examples
#' \dontrun{
#' example.file.path<-system.file("extdata", "BrainMasterList.rda", package = "TrendCatcher")
#' load(example.file.path)
#' gene.symbol.df<-get_GeneEnsembl2Symbol(ensemble.arr = master.list$master.table$Gene)
#  master.table.new<-cbind(master.list$master.table, gene.symbol.df[match(master.list$master.table$Gene, gene.symbol.df$Gene), c("Symbol", "description")])
#  master.list$master.table<-master.table.new
#' time_heatmap<-draw_TimeHeatmap_GO(master.list = master.list)
#' go.terms<-unique(time_heatmap$GO.df$Description)[1:5]
#' time_heatmap_selGO<-draw_TimeHeatmap_selGO(time_heatmap = time_heatmap, sel.go = go.terms, master.list = master.list, GO.perc.thres = 0, nDDEG.thres = 0, save.tiff.path = NA)
#' }
#' 
#' @export
#'
#'


draw_TimeHeatmap_selGO<-function(time_heatmap, sel.go, master.list, GO.perc.thres =0, nDDEG.thres = 0, 
                                 term.width = 80, figure.title = "",
                                 save.tiff.path = NA, tiff.res = 100, tiff.width = 1500, tiff.height =1500){
  if(FALSE){
    # for test only
    sel.go<-c("response to lipopolysaccharide",
              "response to interferon-beta",
              "cytokine-mediated signaling pathway",
              "response to interferon-gamma",
              "response to virus",
              "leukocyte migration",
              "mitotic nuclear division",
              "regulation of vasculature development",
              "extracellular structure organization",
              "regulation of epithelial cell proliferation")
    nDDEG.thres=10
    GO.perc.thres=0.1
  }
  
  ### 1. Get the time array
  t.arr <- master.list$t.arr
  ### 2. Get the time unit
  t.unit <- master.list$time.unit
  #### 3. order time window
  start.t.arr<-paste0(t.arr[1:(length(t.arr)-1)], t.unit)
  end.t.arr<-paste0(t.arr[2:length(t.arr)],t.unit)
  col.name.order<-paste0(start.t.arr, "-", end.t.arr)
  GO.df<-time_heatmap$GO.df
  GO.df<-GO.df %>% dplyr::filter(Description %in% sel.go & nDDEG > nDDEG.thres & perc > GO.perc.thres)
  ################# Prepare mat1 
  sub.GO.df<-GO.df[,c("Description", "t.name", "Avg_log2FC", "nDDEG", "n_background")]
  sub.GO.mat<-dcast(sub.GO.df, formula = Description~t.name, value.var = "Avg_log2FC")
  sub.GO.mat<-sub.GO.mat[,c("Description",col.name.order)]
  sub.GO.mat<-sub.GO.mat[match(unique(GO.df$Description), sub.GO.mat$Description),]
  rownames(sub.GO.mat)<-sub.GO.mat$Description
  sub.GO.mat$Description<-NULL
  sub.GO.mat<-as.matrix(round(sub.GO.mat,2))
  sub.GO.mat<-sub.GO.mat[order(sub.GO.mat[,1], decreasing = T),]
  
  text.GO.mat<-sub.GO.mat
  text.GO.mat<-replace(text.GO.mat, is.na(text.GO.mat), "")
  
  col_fun<-colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  
  ################# Prepare mat2
  mat2<-sub.GO.df$nDDEG[match(rownames(sub.GO.mat), sub.GO.df$Description)]/sub.GO.df$n_background[match(rownames(sub.GO.mat), sub.GO.df$Description)]
  mat2<-as.matrix(round(mat2*100,1))
  colnames(mat2)<-"%GO"
  rownames(mat2)<-rownames(sub.GO.mat)
  col_fun2 = colorRamp2(c(min(mat2),max(mat2)), c("white", "grey"))
  ################# Prepare mat3
  mat3<-as.matrix(sub.GO.df$nDDEG[match(rownames(sub.GO.mat), sub.GO.df$Description)])
  colnames(mat3)<-"nDDEG"
  rownames(mat3)<-rownames(sub.GO.mat)
  col_fun3 = colorRamp2(c(min(mat3),max(mat3)), c("white", "grey"))
  ############### Wrap super long row names
  rownames(sub.GO.mat)<-str_wrap(rownames(sub.GO.mat), width = term.width)
  rownames(mat2)<-str_wrap(rownames(mat2), width = term.width)
  rownames(mat3)<-str_wrap(rownames(mat3), width = term.width)
  
  h1<-Heatmap(sub.GO.mat, na_col = "transparent", cluster_rows = F, cluster_columns = F, 
              row_names_side = "left", column_names_side = "top", column_names_rot = 45, column_names_centered = T,
              rect_gp = gpar(col = "black", lwd = 0.5), row_names_max_width = unit(80, "cm"),
              name = "Ave_log2FC", col = col_fun,
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(text.GO.mat[i,j], x, y)
              })
  h2<-Heatmap(mat2, na_col = "transparent", cluster_rows = F, cluster_columns = F, show_row_names = F,
              column_names_side = "top", column_names_rot = 0, column_names_centered = T,
              rect_gp = gpar(col = "black", lwd = 0.5), 
              show_heatmap_legend = F, col = col_fun2,
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(mat2[i,j], x, y)
              })
  
  h3<-Heatmap(mat3, na_col = "transparent", cluster_rows = F, cluster_columns = F, show_row_names = F,
              column_names_side = "top", column_names_rot = 0, column_names_centered = T,
              rect_gp = gpar(col = "black", lwd = 0.5), 
              show_heatmap_legend = F, col = col_fun3,
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(mat3[i,j], x, y)
              })
  
  p<-h1+h2+h3
  if(is.na(save.tiff.path)){
    print(p)
  } else{
    tiff(filename = save.tiff.path, res = tiff.res, width = tiff.width, height = tiff.height)
    print(p)
    dev.off()
  }
  return(list(time.heatmap = p, merge.df = time_heatmap$merge.df, GO.df = GO.df))
}
