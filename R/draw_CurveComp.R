#' Draw DDEGs trajectories from two master.list object and also show LOESS curve fitting
#'
#' This funcitons takes the master.list output from run_TrendCatcher, and merge.df output from draw_TimeHeatmap_GO and draw_TimeHeatmap_enrichR.
#' And showing all the genes used for enrichment analysis and their logFC compared to previous break point.
#'
#' @param master.list.1, a list object. The output from run_TrendCatcher function, contains master.table element.
#' @param master.list.2, a list object. The output from run_TrendCatcher function, contains master.table element.
#' @param ht.1, TimeHeatmap object. The output from draw_TimeHeatmap_GO function, contains  GO.df object.
#' @param pathway, characters. Must be a biological pathway from GO.df, Description column.
#' @param group.1.name, characters. For example, severe group.
#' @param group.2.name, characters. For example, moderate group.

#' @return A ggplot object and plot.
#' @export
#'
#'


draw_CurveComp<-function(master.list.1, master.list.2, ht.1, pathway, group.1.name, group.2.name){
  GO.df<-ht.1$GO.df
  if(is.null(GO.df)){stop("GO.df must be in the timeheatmap object!")}
  ### Check pathway exists
  if(!pathway %in% GO.df$Description){stop("Selected pathway must be in the timeheatmap object!!!")}
  
  ###  Extract DDEGs from both master.list object
  sub<-GO.df %>% filter(Description == pathway)
  sub<-sub[1,]
  gene.arr.up<-as.character(str_split(sub$geneID_up, "/", simplify = T))
  gene.arr.down<-as.character(str_split(sub$geneID_down, "/", simplify = T))
  gene.arr<-c(gene.arr.up, gene.arr.down)
  gene.arr<-unique(gene.arr[gene.arr!=""])
  n.gene<-length(gene.arr)
  cat(paste0("Found ", n.gene, " DDEGs from ", pathway))
  
  ### Center the data to baseline, so we can compare the curves
  severe<-return.center(master.list = master.list.1, gene.arr = gene.arr)
  moderate<-return.center(master.list = master.list.2, gene.arr = gene.arr)

  ### Write index order for each data

  rep.times<-as.numeric(table(severe$Symbol))
  rep.severe<-rep.int(seq(1, length(rep.times)), times = rep.times)
  
  rep.times<-as.numeric(table(moderate$Symbol))
  rep.moderate<-rep.int(seq(1, length(rep.times)), times = rep.times)
  rep.moderate<-rep.moderate + max(rep.severe)
  
  
  ############################ LOWESS curve fitting]
  ## prepare df
  ID<-c(rep.severe, rep.moderate)
  Count<-c(severe$Center.logFC, moderate$Center.logFC)
  Time<-c(severe$Time, moderate$Time)
  Group<-c(rep(group.1.name,nrow(severe)),rep(group.2.name, nrow(moderate)))
  points = seq(min(Time), max(Time), length.out = 100)
  Group = as.character(Group)
  Count = Count + 1e-8
  df = data.frame(Count = Count, Time = Time, Group = Group, ID = ID)
  df$Group<-factor(df$Group, levels = c(group.1.name, group.2.name))
  
  ## prepare dat
  group.1<-df %>% filter(Group == group.1.name)
  group.2<-df %>% filter(Group == group.2.name)
  mod.1 = loess(Count ~ Time, data = group.1)
  mod.2 = loess(Count ~ Time, data = group.2)
  est.1 = predict(mod.1, data.frame(Time = points), se = TRUE)
  est.2 = predict(mod.2, data.frame(Time = points), se = TRUE)
  dd.1 = data.frame(Time = points, Count = est.1$fit, Group = group.1.name)
  dd.2 = data.frame(Time = points, Count = est.2$fit, Group = group.2.name)
  
  ### Draw trajectories
  dat<-rbind(dd.1, dd.2)
  dat$Group<-factor(dat$Group, levels = c(group.1.name, group.2.name))
  
  g<-ggplot()+geom_path(data = df, 
              mapping = aes(x = Time, y = Count, group = ID, color = Group), alpha = 0.1) +
    geom_path(data = dat, mapping = aes(x = Time, y = Count, color = Group), size = 2) +
    theme_minimal() + ggtitle(pathway)
  return(g)
}













