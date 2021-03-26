#' Draw Time-Heatmap Using enrichR
#'
#' This funcitons takes the master.list output from run_TrendCatcher. And apply a time window sliding strategy
#' to capture all the genes increased/decreased compared to its previous break point, and apply enrichR enrichment
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
#' @param dbs, must one of the enrichR supported database name. To check the list, run dbs <- listEnrichrDbs() command.
#' By default is "BioPlanet_2019".
#' @param id.ensembl, a logic variable. If using ensembl as the keytype. This must match the row name of master.table.
#' By default is TRUE.
#' @param ont, one of "BP", "MF", and "CC" subontologies, or "ALL" for all three. By default is "BP".
#' @param term.width, an integer variable. The character length for each GO term. If one GO term is super long, we can wrap
#' it into term.width of strings into multiple rows. By default if 80.
#' @param GO.enrich.p, an numeric variable. The GO enrichment p-value threshold. By default if 0.05.
#' @param figure.title, a character variable. The main title of Time-Heatmap.
#' @param save.pdf.path, a character variable. If need to save the figure into PDF file. This must be an absolute
#' file path. If not needed save as PDF file, set it to NA. By defualt iS NA.
#' @param pdf.width, a numeric variable. The width of PDF file. By default is 13.
#' @param pdf.height, a numeric variable. The height of PDF file. By default is 15.
#'
#' @return A list object, including elements names merge.df and time.heatmap.
#' time.heatmap is the ggplot object. merge.df includes all the enrichR enrichment result and activation/deactivation time.
#'
#' @examples
#' \dontrun{
#' example.file.path<-system.file("extdata", "BrainMasterList.rda", package = "TrendCatcher")
#' load(example.file.path)
#' th.obj<-draw_TimeHeatmap_enrichR(master.list = master.list)
#' print(th.obj$time.heatmap)
#' head(th.obj$merge.df)
#' }
#' @export
#'
#'
draw_TimeHeatmap_enrichR<-function(master.list, logFC.thres = 1, top.n = 10, dyn.gene.p.thres = 0.05,
                                   dbs = "BioPlanet_2019", id.ensembl = TRUE, term.width = 80, OrgDb = "org.Mm.eg.db",
                                   GO.enrich.p = 0.05, figure.title = "", save.pdf.path = NA, pdf.width = 13, pdf.height =15){

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

  # Apply EnrichR for the up-regulated genes
  enrichr.act.t.genes.go<-list()
  for(i in 1:length(act.list)){
    list.name<-names(act.list)[i]
    cat("Processing up-regulated pathways (enrichR) for time window ", list.name, "\n")
    act.genes<-act.list[[i]]$Gene
    # Convert ID to symbol
    if(OrgDb == "org.Mm.eg.db" & id.ensembl == TRUE){
      act.genes <- mapIds(org.Mm.eg.db, keys = act.genes, keytype = "ENSEMBL", column="SYMBOL")
      act.genes <- as.character(act.genes[!is.na(act.genes)])
    }else if(OrgDb == "org.Hg.eg.db"& id.ensembl == TRUE){
      act.genes <- mapIds(org.Hg.eg.db, keys = act.genes, keytype = "ENSEMBL", column="SYMBOL")
      act.genes <- as.character(act.genes[!is.na(act.genes)])
    } else{stop("TrendCatcher only support Human and Mouse ID conversion from ENSEMBL to SYMBOL. \n")}
    enriched <- enrichr(act.genes, dbs)
    enrichr.act.t.genes.go[[list.name]]<-enriched[[1]]
    rm(enriched)
  }
  enrichr.act.top.go.list<-list()
  for(i in 1:length(enrichr.act.t.genes.go)){
    t.go.df<-enrichr.act.t.genes.go[[i]] %>% filter(Adjusted.P.value<=GO.enrich.p) %>% arrange(desc(Combined.Score))
    if(nrow(t.go.df)>0){
      list.name<-names(act.list)[i]
      enrichr.act.top.go.list[[list.name]]<-data.frame(t.go.df[,c("Term","Adjusted.P.value", "Overlap", "Combined.Score", "Genes")],
                                                       t.name = list.name, type = "Activation")
    }
  }
  enrichr.act.top.go<-do.call(rbind, enrichr.act.top.go.list)
  enrichr.act.top.go<-enrichr.act.top.go[!is.na(enrichr.act.top.go$Term),]

  ### 5. Down-regulated pathways
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

  # Apply EnrichR for the down-regulated genes
  enrichr.deact.t.genes.go<-list()
  for(i in 1:length(deact.list)){
    list.name<-names(deact.list)[i]
    cat("Processing down-regulated pathways (enrichR) for time window ", list.name, "\n")
    deact.genes<-deact.list[[i]]$Gene

    # Convert ID to symbol
    if(OrgDb == "org.Mm.eg.db"){
      deact.genes <- mapIds(org.Mm.eg.db, keys = deact.genes, keytype = "ENSEMBL", column="SYMBOL")
      deact.genes <- as.character(deact.genes[!is.na(deact.genes)])
    }else if(OrgDb == "org.Hg.eg.db"){
      deact.genes <- mapIds(org.Hg.eg.db, keys = deact.genes, keytype = "ENSEMBL", column="SYMBOL")
      deact.genes <- as.character(deact.genes[!is.na(deact.genes)])
    } else{stop("TrendCatcher only support Human and Mouse ID conversion from ENSEMBL to SYMBOL. \n")}

    enriched <- enrichr(deact.genes, dbs)
    enrichr.deact.t.genes.go[[list.name]]<-enriched[[1]]
  }

  enrichr.deact.top.go.list<-list()
  for(i in 1:length(enrichr.deact.t.genes.go)){
    t.go.df<-enrichr.deact.t.genes.go[[i]] %>% filter(Adjusted.P.value<=GO.enrich.p) %>% arrange(desc(Combined.Score))
    if(nrow(t.go.df)>0){
      list.name<-names(deact.list)[i]
      enrichr.deact.top.go.list[[list.name]]<-data.frame(t.go.df[,c("Term","Adjusted.P.value", "Overlap", "Combined.Score", "Genes")],
                                                         t.name = list.name, type = "Deactivation")
    }
  }
  enrichr.deact.top.go<-do.call(rbind, enrichr.deact.top.go.list)
  enrichr.deact.top.go<-enrichr.deact.top.go[!is.na(enrichr.deact.top.go$Term),]

  # Combine the go
  enrichr.merge.df<-rbind(enrichr.act.top.go, enrichr.deact.top.go)
  act.topn.term<-NULL
  for(i in 1:length(unique(enrichr.act.top.go$t.name))){
    t.name.i<-as.character(unique(enrichr.act.top.go$t.name)[i])
    sub.df<-enrichr.act.top.go %>% filter(t.name == t.name.i)
    if(nrow(sub.df)<top.n){
      top.n<-nrow(sub.df)
    }
    term.i<-enrichr.act.top.go$Term[which(enrichr.act.top.go$t.name==t.name.i)][1:top.n]
    act.topn.term<-c(act.topn.term, term.i)
  }
  act.topn.term<-unique(act.topn.term)

  deact.topn.term<-NULL
  for(i in 1:length(unique(enrichr.deact.top.go$t.name))){
    t.name.i<-as.character(unique(enrichr.deact.top.go$t.name)[i])
    sub.df<-enrichr.deact.top.go %>% filter(t.name == t.name.i)
    if(nrow(sub.df)<top.n){
      top.n<-nrow(sub.df)
    }
    term.i<-enrichr.deact.top.go$Term[which(enrichr.deact.top.go$t.name==t.name.i)][1:top.n]
    deact.topn.term<-c(deact.topn.term, term.i)
  }
  deact.topn.term<-unique(deact.topn.term)

  go.term<-unique(c(act.topn.term, deact.topn.term))
  cat("Identified", length(go.term), "GO terms for Time Heatmap (enrichR).", "\n")

  dynamic.df<-enrichr.merge.df %>% filter(Term %in% go.term)
  dynamic.df$Term<-str_wrap(dynamic.df$Term, width = term.width)

  dynamic.df$Term<-factor(dynamic.df$Term, levels = unique(dynamic.df$Term))
  dynamic.df$t.name<-factor(dynamic.df$t.name, levels = unique(dynamic.df$t.name))
  nom<-as.numeric(str_split(dynamic.df$Overlap, "/", simplify = T)[,1])
  denom<-as.numeric(str_split(dynamic.df$Overlap, "/", simplify = T)[,2])
  dynamic.df$NewRatio.num<-paste0(round(nom*100/denom,1), "%")

  data1<-dynamic.df %>% filter(type == "Activation")
  data2<-dynamic.df %>% filter(type == "Deactivation")

  data1$Term<-as.character(data1$Term)
  data1$Term<-factor(data1$Term, levels = unique(dynamic.df$Term))
  data1$t.name<-factor(data1$t.name, levels = unique(dynamic.df$t.name))

  data2$Term<-as.character(data2$Term)
  data2$Term<-factor(data2$Term, levels = unique(dynamic.df$Term))
  data2$t.name<-factor(data2$t.name, levels=unique(dynamic.df$t.name))

  heatmap.plot.enrichR<-ggplot() +
    geom_tile(data = data1, aes(x=t.name, y = reorder(Term, plyr::desc(Term)), fill =  Adjusted.P.value), size = 0.2, color = "grey") +
    geom_text(data = data1, aes(x=t.name, y = reorder(Term, plyr::desc(Term)), label = NewRatio.num)) +
    scale_fill_gradient2(low="red",high="white", name = "Up pathway p.adj", trans="log") +
    new_scale_fill() +
    geom_tile(data = data2, aes(x=t.name, y = reorder(Term, plyr::desc(Term)), fill = Adjusted.P.value), size = 0.2, color = "grey") +
    geom_text(data = data2, aes(x=t.name, y = reorder(Term, plyr::desc(Term)),label = NewRatio.num)) +
    scale_fill_gradient2(low="blue",high="white", name = "Down pathway p.adj", trans = "log") +
    scale_x_discrete(position = "top", limits=levels(dynamic.df$t.name)) + xlab(paste0(figure.title, "\n")) + ylab("") +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold")) +
    scale_y_discrete(limits=rev(levels(dynamic.df$Term)))

  colnames(enrichr.merge.df)<-c("Description", "Adjusted.P.value", "Overlap","Combined.Score", "geneID", "t.name", "type")

  if(is.na(save.pdf.path)){
    print(heatmap.plot.enrichR)
  } else{
    ggsave(heatmap.plot.enrichR, filename = save.pdf.path, width = pdf.width, height = pdf.height)
  }
  return(list(time.heatmap = heatmap.plot.enrichR, merge.df = enrichr.merge.df))
}
