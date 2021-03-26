get_GeneTrajPattern<-function(master.list, dyn.p.thres=0.05, time.unit = "h"){

  t.arr<-master.list$t.arr
  cat("Trajectory Pattern Assignment to each gene. \n")

  GeneTrajPattern.df<-plyr::ddply(master.list$fitted.count, .(Gene), function(df){
    gene.name<-df$Gene[1]
    gene.p.adj<-master.list$dynamic.gene.df$dyn.p.val.adj[which(master.list$dynamic.gene.df$Gene == gene.name)]
    if(length(gene.p.adj)==0){
      gene.pattern.df <- data.frame(pattern = "flat", start.idx =0, end.idx =0, p.adj =1, dynTime = "_", dynSign = "-_", start.t = "_", end.t = "_", pattern_str = "flat")
    }else{
      gene.pattern.df<-gene_pattern_assignment(gene.df = df, gene.p.adj = gene.p.adj)
      gene.pattern.df$p.adj<-gene.p.adj
      # If any p-value is not <=dyn.p.thres
      if(sum(df$t.p.val<=dyn.p.thres)==0){
        # Non-dynamic significant through all time
        gene.pattern.df$dynTime<-paste0("_")
        gene.pattern.df$dynSign<-paste0(ifelse(sign(df$Fit.Count[nrow(df)] - df$Fit.Count[1])==1, "+","-"), "_", collapse = "")
      } else{
        # At least 1 dynamic time is observed
        gene.pattern.df$dynTime<-paste0(t.arr[which(df$t.p.val<=dyn.p.thres)],"_", collapse = "")
        gene.pattern.df$dynSign<-paste0(ifelse(sign(df$Fit.Count - df$Fit.Count[1])[which(df$t.p.val<=dyn.p.thres)]==1, "+","-"), "_", collapse = "")
      }
      start.idx.arr<-as.numeric(str_split(gene.pattern.df$start.idx, "_", simplify = T))
      start.idx.arr<-start.idx.arr[-length(start.idx.arr)]
      start.t<-t.arr[start.idx.arr]
      end.idx.arr<-as.numeric(str_split(gene.pattern.df$end.idx, "_", simplify = T))
      end.idx.arr<-end.idx.arr[-length(end.idx.arr)]
      end.t<-t.arr[end.idx.arr]
      gene.pattern.df$start.t<-paste0(start.t, "_", collapse = "")
      gene.pattern.df$end.t<-paste0(end.t, "_", collapse = "")

      ### For readable trajectory pattern string
      if(gene.pattern.df$pattern!="flat"){
        bk.t.arr<-str_split(paste0(gene.pattern.df$start.t, gene.pattern.df$end.t), "_", simplify = T)[1,]
        bk.t.arr<-unique(bk.t.arr[-length(bk.t.arr)])
        bk.t.arr<-paste0(bk.t.arr, time.unit)

        trend.t.arr<-str_split(gene.pattern.df$pattern, "_", simplify = T)[1,]
        trend.t.arr<-trend.t.arr[-length(trend.t.arr)]
        pattern.string<-rep("", length(c(bk.t.arr, trend.t.arr)))

        pattern.string[seq(1,length(pattern.string),2)]<-bk.t.arr
        pattern.string[seq(2,length(pattern.string),2)]<-trend.t.arr
        pattern.string<-paste0(pattern.string, "_", collapse = "")
        pattern.string<-substr(pattern.string,1,nchar(pattern.string)-1)

        gene.pattern.df$pattern_str<-pattern.string
      } else{
        gene.pattern.df$pattern_str<-"flat"
      }
    }
    return(gene.pattern.df)
  })

}


