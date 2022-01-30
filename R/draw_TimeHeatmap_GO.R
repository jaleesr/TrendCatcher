#' Draw Time-Heatmap Using Gene Ontology (GO) Enrichment
#'
#' This funcitons takes the master.list output from run_TrendCatcher. And apply a time window sliding strategy
#' to capture all the genes increased/decreased compared to its previous break point, and apply GO enrichment
#' analysis.
#'
#' @param master.list, a list object. The output from run_TrendCatcher function, contains master.table element.
#' @param logFC.thres, a numeric variable. The logFC threshold compared to each genes previous break point expression level.
#' By default is 1, meaning for each gene, the current time window's expression level is 2-fold compared to previous
#' break point's expression level.
#' @param top.n, an integer variable. The top N GO enrichment term need to be shown in the Time-Heatmap for up and down
#' regulated pathway. By default is 10. Top 20 GO terms, 10 from up-regulated pathway and 10 from down-regulated pathway
#' will shown in Time-Heatmap.
#' @param dyn.gene.p.thres, a numeric variable. The DDEGs dynamic p-value threshold. By default is 0.05.
#' @param keyType, must be either ENSEMBL or SYMBOL. The row names of your master.list$master.table.
#' @param OrgDb, must be either "org.Mm.eg.db" or "org.Hg.eg.db". Currently only support mouse and human GO annotation database.
#' @param ont, one of "BP", "MF", and "CC" subontologies, or "ALL" for all three. By default is "BP".
#' @param term.width, an integer variable. The character length for each GO term. If one GO term is super long, we can wrap
#' it into term.width of strings into multiple rows. By default if 80.
#' @param GO.enrich.p, an numeric variable. The GO enrichment p-value threshold. By default if 0.05.
#' @param figure.title, a character variable. The main title of Time-Heatmap.
#' @param save.pdf.path, a character variable. If need to save the figure into PDF file. This must be an absolute
#' file path. If not needed save as PDF file, set it to NA. By defualt iS NA.
#' @param pdf.width, a numeric variable. The width of PDF file. By default is 15.
#' @param pdf.height, a numeric variable. The height of PDF file. By default is 15.
#'
#' @return A list object, including elements names merge.df and time.heatmap.
#' time.heatmap is the ggplot object. merge.df includes all the GO enrichment result and activation/deactivation time.
#'
#' @examples
#' \dontrun{
#' example.file.path<-system.file("extdata", "BrainMasterList.rda", package = "TrendCatcher")
#' load(example.file.path)
#' gene.symbol.df<-get_GeneEnsembl2Symbol(ensemble.arr = master.list$master.table$Gene)
#  master.table.new<-cbind(master.list$master.table, gene.symbol.df[match(master.list$master.table$Gene, gene.symbol.df$Gene), c("Symbol", "description")])
#  master.list$master.table<-master.table.new
#' th.obj<-draw_TimeHeatmap_GO(master.list = master.list)
#' print(th.obj$time.heatmap)
#' head(th.obj$merge.df)
#' }
#' @export
#'
#'
draw_TimeHeatmap_GO<-function(master.list, logFC.thres = 1, top.n = 10, dyn.gene.p.thres = 0.05,
                              keyType = "SYMBOL", OrgDb = "org.Mm.eg.db", ont = "BP", term.width = 80,
                              GO.enrich.p = 0.05, figure.title = "", save.tiff.path = NA, tiff.res = 100, tiff.width = 1500, tiff.height =1500){
  if(FALSE){
    # for testing
    logFC.thres = 0.0
    top.n = 10
    dyn.gene.p.thres = 0.05
    keyType = "SYMBOL"
    OrgDb = "org.Mm.eg.db"
    ont = "BP"
    GO.enrich.p = 0.05
    figure.title = ""
    term.width<-80
  }
  
  ####### Check If there is the Symbol column exist!!!! ######
  idx<-grep("Symbol", colnames(master.list$master.table))
  if(length(idx)==0){stop("Please add Symbol column to your master.list$master.table. By default, it should be a column of gene SYMBOLs!")}
  
  
  ############ If keyType is SYMBOL, but Symbol column is ENSEMBL, need to stop!!!
  # Check Symbol column has many genes start with EN
  txt.arr<-master.list$master.table$Symbol
  txt.initial.2<-substr(x = txt.arr, start = 1, stop = 2)
  en.num<-sum(txt.initial.2 == "EN", na.rm = T)
  if(en.num>100 & keyType == "SYMBOL"){stop("Your Symbol column contains ENSEMBL, please change keyType to ENSEMBL!!!") }
  if(en.num==0 & keyType != "SYMBOL"){stop("Your Symbol column contains SYMBOL, please change keyType to SYMBOL!!!") }
  
  
  ### 1. Get the time array
  t.arr <- master.list$t.arr
  ### 2. Get the time unit
  t.unit <- master.list$time.unit
  ### 3. Filter out only dyn-DDEGs
  if (is.na(t.arr) || is.na(t.unit)) 
    stop("Master.list needs time unit and time array.")
  dyn.gene.pattern <- master.list$master.table %>% filter(dyn.p.val.adj <= 
                                                            dyn.gene.p.thres)

  ### 4. Up-regulated pathways genes
  # Loop by time window to filter out genes going up within each time window
  act.list <- list()
  for (i in 1:(length(t.arr) - 1)) {
    start.t.thres <- t.arr[i]
    end.t.thres <- t.arr[i + 1]
    list.name <- paste0(start.t.thres, t.unit, "-", end.t.thres, 
                        t.unit)
    cat("Processing up-regulated genes for time window ", 
        list.name, "\n")
    act.list[[list.name]] <- plyr::ddply(dyn.gene.pattern, 
                                         .(Gene), function(df) {
                                           idx <- grep("up", str_split(df$pattern, "_", 
                                                                       simplify = T))
                                           if (length(idx) != 0) {
                                             start.idx <- as.numeric(str_split(df$start.idx, 
                                                                               "_", simplify = T)[idx])
                                             end.idx <- as.numeric(str_split(df$end.idx, 
                                                                             "_", simplify = T)[idx])
                                             start.t <- t.arr[start.idx]
                                             end.t <- t.arr[end.idx]
                                             act.flag <- 0
                                             for (j in 1:length(start.t)) {
                                               if (start.t[j] <= start.t.thres & end.t[j] >= 
                                                   end.t.thres) {
                                                 data.trans <- master.list$fitted.count %>% 
                                                   filter(Gene == df$Gene)
                                                 end.t.thres.count <- data.trans$Fit.Count[which(data.trans$Time == 
                                                                                                   end.t.thres)]
                                                 bk.arr <- str_split(df$start.t, "_", simplify = T)[1, 
                                                 ]
                                                 bk.arr <- as.numeric(bk.arr[-length(bk.arr)])
                                                 previous.bk.t <- max(bk.arr[bk.arr <= 
                                                                               start.t.thres])
                                                 prev.bk.count <- data.trans$Fit.Count[which(data.trans$Time == 
                                                                                               previous.bk.t)]
                                                 if (log(end.t.thres.count, base = 2) - 
                                                     log(prev.bk.count, base = 2) > logFC.thres) {
                                                   act.flag <- 1
                                                 }
                                               }
                                             }
                                             if (act.flag == 1) {
                                               return(df)
                                             }
                                           }
                                         })
  }
  act.df.pattern <- do.call(rbind, act.list)
  ### 5. Up-regulated pathway enrichment
  act.go.list <- list()
  for (i in 1:length(act.list)) {
    list.name <- names(act.list)[i]
    cat("Processing up-regulated pathways for time window ", 
        list.name, "\n")
    act.genes <- act.list[[i]]$Symbol
    act.genes <- unique(act.genes[act.genes!=""])
    act.go.list[[list.name]] <- enrichGO(act.genes, keyType = keyType, 
                                         OrgDb = OrgDb, ont = ont)
  }
  ### 6. Get enrichment GOs within GO.enrich.p threshold
  act.top.go.list <- list()
  for (i in 1:length(act.go.list)) {
    t.go.df <- act.go.list[[i]]@result %>% filter(p.adjust <= 
                                                    GO.enrich.p)
    list.name <- names(act.go.list)[i]
    if (nrow(t.go.df) > 0) {
      act.top.go.list[[list.name]] <- data.frame(t.go.df[, 
                                                         c("ID", "Description", "p.adjust", "GeneRatio", 
                                                           "BgRatio", "geneID")], t.name = list.name, 
                                                 type = "Activation")
    }
  }
  act.df <- do.call(rbind, act.top.go.list)
  act.df <- act.df[!is.na(act.df$ID), ]
  
  ### 7. Down-regulated pathways genes
  deact.list <- list()
  for (i in 1:(length(t.arr) - 1)) {
    start.t.thres <- t.arr[i]
    end.t.thres <- t.arr[i + 1]
    list.name <- paste0(start.t.thres, t.unit, "-", end.t.thres, 
                        t.unit)
    cat("Processing down-regulated genes for time window ", 
        list.name, "\n")
    deact.list[[list.name]] <- plyr::ddply(dyn.gene.pattern, 
                                           .(Gene), function(df) {
                                             idx <- grep("down", str_split(df$pattern, "_", 
                                                                           simplify = T))
                                             if (length(idx) != 0) {
                                               start.idx <- as.numeric(str_split(df$start.idx, 
                                                                                 "_", simplify = T)[idx])
                                               end.idx <- as.numeric(str_split(df$end.idx, 
                                                                               "_", simplify = T)[idx])
                                               start.t <- t.arr[start.idx]
                                               end.t <- t.arr[end.idx]
                                               act.flag <- 0
                                               for (j in 1:length(start.t)) {
                                                 if (start.t[j] <= start.t.thres & end.t[j] >= 
                                                     end.t.thres) {
                                                   data.trans <- master.list$fitted.count %>% 
                                                     filter(Gene == df$Gene)
                                                   end.t.thres.count <- data.trans$Fit.Count[which(data.trans$Time == 
                                                                                                     end.t.thres)]
                                                   bk.arr <- str_split(df$start.t, "_", simplify = T)[1, 
                                                   ]
                                                   bk.arr <- as.numeric(bk.arr[-length(bk.arr)])
                                                   previous.bk.t <- max(bk.arr[bk.arr <= 
                                                                                 start.t.thres])
                                                   prev.bk.count <- data.trans$Fit.Count[which(data.trans$Time == 
                                                                                                 previous.bk.t)]
                                                   if (log(prev.bk.count, base = 2) - log(end.t.thres.count, 
                                                                                          base = 2) > logFC.thres) {
                                                     act.flag <- 1
                                                   }
                                                 }
                                               }
                                               if (act.flag == 1) {
                                                 return(df)
                                               }
                                             }
                                           })
  }
  deact.df.pattern <- do.call(rbind, deact.list)
  ### 8. Down regultaed genes enrichment pathway
  deact.go.list <- list()
  for (i in 1:length(deact.list)) {
    list.name <- names(deact.list)[i]
    cat("Processing down-regulated pathways for time window ", 
        list.name, "\n")
    deact.genes <- deact.list[[i]]$Symbol
    deact.genes <- unique(deact.genes[deact.genes!=""])
    deact.go.list[[list.name]] <- enrichGO(deact.genes, 
                                           keyType = keyType, OrgDb = OrgDb, ont = ont)
  }
  ### 9. Down pathways within GO.enrich.p
  deact.top.go.list <- list()
  for (i in 1:length(deact.go.list)) {
    list.name <- names(deact.go.list)[i]
    t.go.df <- deact.go.list[[i]]@result %>% filter(p.adjust <= 
                                                      GO.enrich.p)
    if (nrow(t.go.df) > 0) {
      deact.top.go.list[[list.name]] <- data.frame(t.go.df[, 
                                                           c("ID", "Description", "p.adjust", "GeneRatio", 
                                                             "BgRatio", "geneID")], t.name = list.name, 
                                                   type = "Deactivation")
    }
  }
  deact.df <- do.call(rbind, deact.top.go.list)
  deact.df <- deact.df[!is.na(deact.df$ID), ]
  merge.df <- rbind(act.df, deact.df) ## Contains all the enrichment GOs!!!!!!!!!!!!
  
  ##################################### Get top 10 up and down GOs within each time window, 
  ############################# the total can be less than N time window * 10 because overlap
  act.topn.term <- NULL
  for (i in 1:length(unique(act.df$t.name))) {
    t.name.i <- as.character(unique(act.df$t.name)[i])
    sub.df <- act.df %>% filter(t.name == t.name.i)
    if (nrow(sub.df) < top.n) {
      top.n <- nrow(sub.df)
    }
    term.i <- act.df$Description[which(act.df$t.name == 
                                         t.name.i)][1:top.n]
    act.topn.term <- c(act.topn.term, term.i)
  }
  act.topn.term <- unique(act.topn.term)
  deact.topn.term <- NULL
  for (i in 1:length(unique(deact.df$t.name))) {
    t.name.i <- as.character(unique(deact.df$t.name)[i])
    sub.df <- deact.df %>% filter(t.name == t.name.i)
    if (nrow(sub.df) < top.n) {
      top.n <- nrow(sub.df)
    }
    term.i <- deact.df$Description[which(deact.df$t.name == 
                                           t.name.i)][1:top.n]
    deact.topn.term <- c(deact.topn.term, term.i)
  }
  deact.topn.term <- unique(deact.topn.term)
  go.term <- unique(c(act.topn.term, deact.topn.term))
  
  ##################### For each go, calculate average logFC t-t-1 to define the break point of GO #########
  logFC.mean.arr<-NULL
  sub.merge.df<-merge.df %>% filter(Description %in% go.term)
  #### Calculate log2FC within each time window for each GO
  GO.list<-list()
  counter<-1
  for(i in 1:length(go.term)){
    # each GO
    go.i<-go.term[i]
    # for each GO, get candidate up and down genes
    sub.merge.df<-merge.df %>% filter(Description == go.i)
    for(j in 1:(length(t.arr)-1)){
      # each time window
      start.t.thres <- t.arr[j]
      end.t.thres <- t.arr[j + 1]
      list.name <- paste0(start.t.thres, t.unit, "-", end.t.thres, 
                          t.unit)
      sub.t<-sub.merge.df %>% filter(t.name == list.name) # 0 row, 1 row up/down, 2 row mix
      if(nrow(sub.t)!=0){
        sub.t.up<-sub.t %>% filter(type == "Activation")
        sub.t.down<-sub.t %>% filter(type == "Deactivation")
        if(nrow(sub.t.up)!=0){
          n_up<-as.numeric(str_split(sub.t.up$GeneRatio, "/", simplify = T)[1])
          geneID_up<-sub.t.up$geneID
          p.adjust.up<-sub.t.up$p.adjust
          sel.genes.up<-paste0(str_split(sub.t.up$geneID, "/", simplify = T))
          sel.genes.up<-unique(sel.genes.up[sel.genes.up!=""])
          master.list$fitted.count$Symbol<-master.list$master.table$Symbol[match(master.list$fitted.count$Gene, master.list$master.table$Gene)]
          logFC.arr.up<-NULL
          for(k in 1:length(sel.genes.up)){
            # fitted count change
            gene.i<-sel.genes.up[k]
            count.df<-master.list$fitted.count %>% filter(Symbol == gene.i)
            logFC<-log(count.df$Fit.Count[which(count.df$Time == end.t.thres)], 2) - log(count.df$Fit.Count[which(count.df$Time == start.t.thres)],2)
            logFC.arr.up<-c(logFC.arr.up, logFC)
          }
          logFC.arr.up<-logFC.arr.up[logFC.arr.up>0]
        }else{
          n_up<-0
          geneID_up<-""
          p.adjust.up<-""
          logFC.arr.up<-NULL
        }
        if(nrow(sub.t.down)!=0){
          n_down<-as.numeric(str_split(sub.t.down$GeneRatio, "/", simplify = T)[1])
          geneID_down<-sub.t.down$geneID
          p.adjust.down<-sub.t.down$p.adjust
          sel.genes.down<-paste0(str_split(sub.t.down$geneID, "/", simplify = T))
          sel.genes.down<-unique(sel.genes.down[sel.genes.down!=""])
          master.list$fitted.count$Symbol<-master.list$master.table$Symbol[match(master.list$fitted.count$Gene, master.list$master.table$Gene)]
          logFC.arr.down<-NULL
          for(k in 1:length(sel.genes.down)){
            # fitted count change
            gene.i<-sel.genes.down[k]
            count.df<-master.list$fitted.count %>% filter(Symbol == gene.i)
            logFC<-log(count.df$Fit.Count[which(count.df$Time == end.t.thres)], 2) - log(count.df$Fit.Count[which(count.df$Time == start.t.thres)],2)
            logFC.arr.down<-c(logFC.arr.down, logFC)
          }
          logFC.arr.down<-logFC.arr.down[logFC.arr.down<0]
        }else{
          n_down<-0
          geneID_down<-""
          p.adjust.down<-""
          logFC.arr.down<-NULL
        }
        logFC.arr<-c(logFC.arr.up, logFC.arr.down)
        
        logFC.mean<-mean(logFC.arr)
        direction<-ifelse(logFC.mean>0, "Activation", "Deactivation")
        GO.list[[counter]]<-data.frame(ID = sub.merge.df$ID[1], Description = go.i, t.name = list.name, direction = direction, 
                                       Avg_log2FC = logFC.mean, n_total = n_up+n_down, 
                                       n_background = as.numeric(str_split(sub.merge.df$BgRatio, "/", simplify = T)[1]),
                                       n_up = n_up, n_down = n_down, geneID_up = geneID_up, geneID_down = geneID_down, p.adjust.up = p.adjust.up, p.adjust.down = p.adjust.down)
        counter<-counter+1
      }
    }
  }
  GO.df<-do.call(rbind, GO.list)
  GO.df.info<-ddply(GO.df, .(Description), function(df){
    up.genes<-as.character(str_split(df$geneID_up, "/", simplify = T))
    up.genes<-up.genes[up.genes!=""]
    down.genes<-as.character(str_split(df$geneID_down, "/", simplify = T))
    down.genes<-down.genes[down.genes!=""]
    nDDEG<-length(unique(c(up.genes, down.genes)))
    DDEGs<-paste0(unique(c(up.genes, down.genes)), "/", collapse = "")
    return(data.frame(nDDEG = nDDEG, DDEGs = DDEGs))
  })
  GO.df$nDDEG<-GO.df.info$nDDEG[match(GO.df$Description, GO.df.info$Description)]
  GO.df$DDEGs<-GO.df.info$DDEGs[match(GO.df$Description, GO.df.info$Description)]
  GO.df$perc<-GO.df$nDDEG/GO.df$n_background  ####### GO.df contains log2FC for each selected GO!!!!!!
  
  ################# Prepare for complex heatmap ###############
  start.t.arr<-paste0(t.arr[1:(length(t.arr)-1)], t.unit)
  end.t.arr<-paste0(t.arr[2:length(t.arr)],t.unit)
  col.name.order<-paste0(start.t.arr, "-", end.t.arr)
  
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
  p<-draw(p, column_title = figure.title,
          column_title_gp = gpar(fontsize = 16))
  if(is.na(save.tiff.path)){
    print(p)
  } else{
    tiff(filename = save.tiff.path, res = tiff.res, width = tiff.width, height = tiff.height)
    print(p)
    dev.off()
  }
  return(list(time.heatmap = p, merge.df = merge.df, GO.df = GO.df))
}
