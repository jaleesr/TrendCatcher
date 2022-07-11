#' Draw DDEGs trajectories from two master.list object and also show LOESS curve fitting with permuataon
#' @export
#'
#'
draw_CurveComp_Perm<-function(master.list.1, master.list.2, ht.1, pathway, 
                              group.1.name, group.2.name, n.perm = 500, parall = F, 
                              pvalue.threshold = 0.05){
  if(FALSE){
    master.list.1 = master.list.severe
    master.list.2 = master.list.moderate 
    ht.1 = ht.severe
    pathway = "neutrophil activation" 
    group.1.name = "severe" 
    group.2.name = "moderate" 
    n.perm = 10 
    parall = F 
    pvalue.threshold = 0.05
  }
  # update
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
  cat(paste0("Found ", n.gene, " DDEGs from ", pathway, "/n"))
  
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
  #df$Group<-factor(df$Group, levels = c(group.1.name, group.2.name))
  
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
  #dat$Group<-factor(dat$Group, levels = c(group.1.name, group.2.name))
  
  ############## Permutation
  area = intervalArea.new(dd.1 = dd.1, dd.2 = dd.2)
  #perm  = permutation.new(perm.dat = df, n.perm = n.perm, points = points, parall = parall)
  
  ## Start permutation
  cat("Start Permutation \n")
  pp = list() 
  perm = 0 # to be able to store the value
  n.subjects = length(unique(df$ID))
  #cat("# of Subjects = ", n.subjects, "\n")
  aggregate.df<-df
  
  ## Run in Parallel
  if(parall == TRUE) {
    max.cores = detectCores()
    cat("# cores = ", max.cores, "\n")
    desired.cores = max.cores - 1		
    cl = makeCluster(desired.cores)
    registerDoParallel(cl)
  } 
  
  size<-max(aggregate.df$ID)
  pp = llply(1:n.perm, function(j){
    perm.dat.1<-aggregate.df
    group.perm<-c(rep(1, round(size/2)), rep(2, size - round(size/2)))
    group.perm[sample(seq(1, round(size/2)), round(size*0.05), replace = F)]<-2
    
    group.perm[sample(seq(round(size/2)+1, size), round(size*0.05), replace = F)]<-1
    
    rep.times<-as.numeric(table(aggregate.df$ID))
    
    new.group.id<-rep.int(x = group.perm, times = rep.times)
    perm.dat.1$Group<-new.group.id
    g.1 = perm.dat.1[perm.dat.1$Group == 1, ]
    g.2 = perm.dat.1[perm.dat.1$Group == 2, ]
    g.min = max(sort(g.1$Time)[1], sort(g.2$Time)[1])
    g.max = min(sort(g.1$Time)[length(g.1$Time)], sort(g.2$Time)[length(g.2$Time)])
    
    
    if(g.min > min(points) | g.max < max(points))
    {
      cat("Special Case: generated permutation is out of range \n")
      assign(paste("Model", j, sep = "_"), NULL)
    } else{
      group.1<-perm.dat.1 %>% filter(Group == 1)
      group.2<-perm.dat.1 %>% filter(Group == 2)
      mod.1 = loess(Count ~ Time, data = group.1)
      mod.2 = loess(Count ~ Time, data = group.2)
      est.1 = predict(mod.1, data.frame(Time = points), se = TRUE)
      est.2 = predict(mod.2, data.frame(Time = points), se = TRUE)
      dd.1 = data.frame(Time = points, Count = est.1$fit, Group = group.1.name)
      dd.2 = data.frame(Time = points, Count = est.2$fit, Group = group.2.name)
      perm = list(dd.1, dd.2)
      #perm = curveFitting.new(perm.dat.1= perm.dat.1, points = points)
      assign(paste("Model", j, sep = "_"), perm)
    }
  }, .parallel = parall, .progress = "text", .inform = TRUE,
  .paropts = list(.export=ls(.GlobalEnv),
                  .packages=.packages(all.available=T)))
  
  
  if(parall == TRUE) {
    stopCluster(cl)
  }
  
  pp[sapply(pp, is.null)] = NULL
  
  perm<-pp
  

  area.perm = areaPermutation.new(perm)
  a1 = do.call(rbind, area.perm)
  a2 = do.call(rbind, a1[,2])
  ## Calculate AR p-value 
  pvalue.area = sapply(1:(length(points)-1), function(i){
    sum(a2[,i] >= area$ar.abs[i])/length(a2[,i])
  } )
  
  cat("p-value Adjustment Method = ", "BH", "\n")
  adjusted.pvalue = p.adjust(pvalue.area, method = "BH")
  interval = findSigInterval(adjusted.pvalue, 
                             threshold = pvalue.threshold, sign = area$ar.sign)
  st = points[interval$start]
  en = points[interval$end + 1]
  
  #### Draw ggplot and show significant seperation time intervals
  
  sub.11 = list()
  sub.12 = list()
  xx = NULL
  for(i in 1:length(st))
  {
    sub.11[[i]] = subset(dd.1, Time >= st[i] & Time <= en[i])  
    sub.12[[i]] = subset(dd.2, Time >= st[i] & Time <= en[i])
  }
  
  g<-ggplot()+geom_path(data = df, 
                        mapping = aes(x = Time, y = Count, group = ID, color = Group), alpha = 0.1) +
    geom_path(data = dat, mapping = aes(x = Time, y = Count, color = Group), size = 2) +
    theme_minimal() + ggtitle(pathway) +
    geom_ribbon(data=sub.11[[1]], aes(x = Time, ymin = sub.12[[1]]$Count, ymax = Count), 
                colour= "grey3", fill="grey69", alpha = 0.6)
  result.l<-list(adjusted.pvalue.area = adjusted.pvalue, perm = perm, st = st, en = en, plot = g)
  return(result.l)
}




