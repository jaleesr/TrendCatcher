#' Subset TimeHeatmap
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
  GO.df<-GO.df %>% filter(Description %in% sel.go & nDDEG > nDDEG.thres & perc > GO.perc.thres)
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
