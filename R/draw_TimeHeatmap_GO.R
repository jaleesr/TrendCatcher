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
#' th.obj<-draw_TimeHeatmap_GO(master.list = master.list)
#' print(th.obj$time.heatmap)
#' head(th.obj$merge.df)
#' }
#' @export
#'
#'
draw_TimeHeatmap_GO<-function(master.list, logFC.thres = 1, top.n = 10, dyn.gene.p.thres = 0.05,
                              keyType = "ENSEMBL", OrgDb = "org.Mm.eg.db", ont = "BP", term.width = 80,
                              GO.enrich.p = 0.05, figure.title = "", save.pdf.path = NA, pdf.width = 15, pdf.height =15){
  if(FALSE){
    # For testing
    logFC.thres = 1
    top.n = 10
    dyn.gene.p.thres = 0.05
    keyType = "ENSEMBL"
    OrgDb = "org.Mm.eg.db"
    ont = "BP"
    GO.enrich.p=0.05
    figure.title = "Brain EC"
  }
  ### 1. Get the time array
  t.arr<-master.list$t.arr
  ### 2. Get the time unit
  t.unit<-master.list$time.unit
  if(is.na(t.arr) || is.na(t.unit)) stop("Master.list needs time unit and time array.")
  ### 3. Filter out only dyn-DEGs
  dyn.gene.pattern<-master.list$master.table %>% filter(dyn.p.val.adj<=dyn.gene.p.thres)
  ### 4. Up-regulated pathways genes
  # Loop by time window to filter out genes going up within each time window
  act.list<-list()
  for(i in 1:(length(t.arr)-1)){
    # For each time window
    start.t.thres<-t.arr[i]
    end.t.thres<-t.arr[i+1]
    list.name<-paste0(start.t.thres,t.unit,"-", end.t.thres, t.unit)
    cat("Processing up-regulated genes for time window ", list.name, "\n")
    # For each gene check if it activated within the time range
    act.list[[list.name]]<-plyr::ddply(dyn.gene.pattern, .(Gene), function(df){
      # Select pattern include "up"
      idx<-grep("up", str_split(df$pattern,"_", simplify = T))
      if(length(idx)!=0){
        # Select within the range
        start.idx<-as.numeric(str_split(df$start.idx, "_", simplify = T)[idx])
        end.idx<-as.numeric(str_split(df$end.idx, "_", simplify = T)[idx])
        start.t<-t.arr[start.idx]
        end.t<-t.arr[end.idx]
        act.flag<-0
        # For each range need to check
        for(j in 1:length(start.t)){
          # Make sure within the range
          if(start.t[j]<=start.t.thres & end.t[j]>=end.t.thres){
            # Make sure the end.t.thres logFC is larger than start.t.thres
            data.trans<-master.list$fitted.count %>% filter(Gene == df$Gene)
            end.t.thres.count<-data.trans$Fit.Count[which(data.trans$Time == end.t.thres)]
            # start.t.thres.count need to be the previous bk value
            bk.arr<-str_split(df$start.t, "_", simplify = T)[1,]
            bk.arr<-as.numeric(bk.arr[-length(bk.arr)])
            previous.bk.t<-max(bk.arr[bk.arr<=start.t.thres])
            prev.bk.count<-data.trans$Fit.Count[which(data.trans$Time == previous.bk.t)]
            if(log(end.t.thres.count, base = 2)-log(prev.bk.count, base = 2)>logFC.thres){
              act.flag<-1
            }
          }
        }
        if(act.flag==1){
          return(df)
        }
      }
    })
  }
  act.df.pattern<-do.call(rbind,act.list)
  ### 5. Up-regulated pathway enrichment
  act.go.list<-list()
  for(i in 1:length(act.list)){
    list.name<-names(act.list)[i]
    cat("Processing up-regulated pathways for time window ", list.name, "\n")
    act.genes<-act.list[[i]]$Gene
    act.go.list[[list.name]]<-enrichGO(act.genes, keyType = keyType,
                                       OrgDb = OrgDb,
                                       ont = ont,
                                       universe = master.list$master.table$Gene)
  }
  # Only keep the enriched terms, with GO.enrich.p as threshold
  act.top.go.list<-list()
  for(i in 1:length(act.go.list)){
    t.go.df<-act.go.list[[i]]@result %>% filter(p.adjust<=GO.enrich.p)
    list.name<-names(act.go.list)[i]
    if(nrow(t.go.df)>0){
      act.top.go.list[[list.name]]<-data.frame(t.go.df[,c("ID","Description", "p.adjust", "GeneRatio","BgRatio","geneID")],
                                               t.name = list.name, type = "Activation")
    }
  }
  act.df<-do.call(rbind,act.top.go.list)
  # Remove NA term
  act.df<-act.df[!is.na(act.df$ID),]

  ### 6. Down-regulated pathway genes
  # Loop by time window to filter out genes going up within each time window
  deact.list<-list()
  for(i in 1:(length(t.arr)-1)){
    # Get time range
    start.t.thres<-t.arr[i]
    end.t.thres<-t.arr[i+1]
    list.name<-paste0(start.t.thres,t.unit,"-", end.t.thres, t.unit)
    cat("Processing down-regulated genes for time window ", list.name, "\n")
    # For each gene check if it activated within the time range
    deact.list[[list.name]]<-plyr::ddply(dyn.gene.pattern, .(Gene), function(df){
      idx<-grep("down", str_split(df$pattern,"_", simplify = T))
      if(length(idx)!=0){
        start.idx<-as.numeric(str_split(df$start.idx, "_", simplify = T)[idx])
        end.idx<-as.numeric(str_split(df$end.idx, "_", simplify = T)[idx])
        start.t<-t.arr[start.idx]
        end.t<-t.arr[end.idx]
        act.flag<-0
        # For each range need to check
        for(j in 1:length(start.t)){
          if(start.t[j]<=start.t.thres & end.t[j]>=end.t.thres){
            # Make sure the end.t.thres logFC is larger than start.t.thres
            data.trans<-master.list$fitted.count %>% filter(Gene == df$Gene)
            end.t.thres.count<-data.trans$Fit.Count[which(data.trans$Time == end.t.thres)]
            # start.t.thres.count need to be the previous bk value
            bk.arr<-str_split(df$start.t, "_", simplify = T)[1,]
            bk.arr<-as.numeric(bk.arr[-length(bk.arr)])
            previous.bk.t<-max(bk.arr[bk.arr<=start.t.thres])
            prev.bk.count<-data.trans$Fit.Count[which(data.trans$Time == previous.bk.t)]
            if(log(prev.bk.count, base = 2)-log(end.t.thres.count, base = 2)>logFC.thres){
              act.flag<-1
            }
          }
        }
        if(act.flag==1){
          return(df)
        }
      }
    })
  }
  deact.df.pattern<-do.call(rbind,deact.list)
  ### 7. Down-regulated pathway enrichment
  deact.go.list<-list()
  for(i in 1:length(deact.list)){
    list.name<-names(deact.list)[i]
    cat("Processing down-regulated pathways for time window ", list.name, "\n")
    deact.genes<-deact.list[[i]]$Gene
    deact.go.list[[list.name]]<-enrichGO(deact.genes, keyType = keyType,
                                         OrgDb = OrgDb,
                                         ont = ont, universe = master.list$master.table$Gene)
  }
  # Only keep the enriched terms
  deact.top.go.list<-list()
  for(i in 1:length(deact.go.list)){
    list.name<-names(deact.go.list)[i]
    t.go.df<-deact.go.list[[i]]@result %>% filter(p.adjust<=GO.enrich.p)
    if(nrow(t.go.df)>0){
      deact.top.go.list[[list.name]]<-data.frame(t.go.df[,c("ID","Description", "p.adjust", "GeneRatio","BgRatio", "geneID")],
                                                 t.name = list.name, type = "Deactivation")
    }
  }
  deact.df<-do.call(rbind,deact.top.go.list)
  # Remove NA term
  deact.df<-deact.df[!is.na(deact.df$ID),]

  ################ Combine enriched terms #####################
  merge.df<-rbind(act.df, deact.df)

  act.topn.term<-NULL
  for(i in 1:length(unique(act.df$t.name))){
    t.name.i<-as.character(unique(act.df$t.name)[i])
    sub.df<-act.df %>% filter(t.name == t.name.i)
    if(nrow(sub.df)<top.n){
      top.n<-nrow(sub.df)
    }
    term.i<-act.df$Description[which(act.df$t.name==t.name.i)][1:top.n]
    act.topn.term<-c(act.topn.term, term.i)
  }
  act.topn.term<-unique(act.topn.term)

  deact.topn.term<-NULL
  for(i in 1:length(unique(deact.df$t.name))){
    t.name.i<-as.character(unique(deact.df$t.name)[i])
    sub.df<-deact.df %>% filter(t.name == t.name.i)
    if(nrow(sub.df)<top.n){
      top.n<-nrow(sub.df)
    }
    term.i<-deact.df$Description[which(deact.df$t.name==t.name.i)][1:top.n]
    deact.topn.term<-c(deact.topn.term, term.i)
  }
  deact.topn.term<-unique(deact.topn.term)

  go.term<-unique(c(act.topn.term, deact.topn.term))
  heatmap.go.n<-length(go.term)
  cat("Identified", heatmap.go.n, "GO terms for Time Heatmap.", "\n")
  dynamic.df<-merge.df %>% filter(Description %in% go.term)

  ### 8. Draw Time-heatmap ####################################################
  bg.arr<-dynamic.df$BgRatio
  term.n.total<-str_split(bg.arr, "/", simplify = T)[,1]
  found.arr<-dynamic.df$GeneRatio
  found.n.total<-str_split(found.arr, "/", simplify = T)[,1]

  dynamic.df$NewRatio<-paste0(found.n.total,"/",term.n.total)
  dynamic.df$Description<-str_wrap(dynamic.df$Description,width = term.width)
  dynamic.df$NewRatio.num<-paste0(round(as.numeric(found.n.total)*100/as.numeric(term.n.total),1), "%")

  dynamic.df$Description<-factor(dynamic.df$Description, levels = unique(dynamic.df$Description))
  dynamic.df$t.name<-factor(dynamic.df$t.name, levels = unique(dynamic.df$t.name))
  dynamic.df$p.adjust.round<-signif(dynamic.df$p.adjust, digits = 3)

  data1<-dynamic.df %>% filter(type == "Activation")
  data2<-dynamic.df %>% filter(type == "Deactivation")

  data1$Description<-as.character(data1$Description)
  data1$Description<-factor(data1$Description, levels = unique(dynamic.df$Description))
  data1$t.name<-factor(data1$t.name, levels = unique(dynamic.df$t.name))


  data2$Description<-as.character(data2$Description)
  data2$Description<-factor(data2$Description, levels = unique(dynamic.df$Description))
  data2$t.name<-factor(data2$t.name, levels=unique(dynamic.df$t.name))

  time.heatmap<-ggplot() +
    geom_tile(data = data1, aes(x=t.name, y = reorder(Description, plyr::desc(Description)), fill =  p.adjust.round), size = 0.2, color = "grey") +
    geom_text(data = data1, aes(x=t.name, y = reorder(Description, plyr::desc(Description)), label = NewRatio.num)) +
    scale_fill_gradient2(low="red",high="white", name = "Up pathway p.adj", trans="log") +
    ggnewscale::new_scale_fill() +
    geom_tile(data = data2, aes(x=t.name, y = reorder(Description, plyr::desc(Description)), fill = p.adjust.round), size = 0.2, color = "grey") +
    geom_text(data = data2, aes(x=t.name, y = reorder(Description, plyr::desc(Description)),label = NewRatio.num)) +
    scale_fill_gradient2(low="blue",high="white", name = "Down pathway p.adj", trans = "log") +
    scale_x_discrete(position = "top", limits=levels(dynamic.df$t.name)) + xlab(paste0(figure.title, "\n")) + ylab("") +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold")) +
    scale_y_discrete(limits=rev(levels(dynamic.df$Description)))
  if(is.na(save.pdf.path)){
    print(time.heatmap)
  } else{
    ggsave(time.heatmap, filename = save.pdf.path, width = pdf.width, height = pdf.height)
  }
  return(list(time.heatmap = time.heatmap, merge.df = merge.df))
}
