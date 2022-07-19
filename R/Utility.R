#' Get the time array from count table
#'
#' It take the count table, grep the time element from column name.
#' The column name of the count table must satisfy \code{Prj_Time_Rep1} format.
#'
#' @param raw.count.df, a count table.
#' @return a numeric array of each sample's time.
#' @import RColorBrewer circlize compiler doSNOW
#' @import enrichR foreach ggnewscale ggplot2 gridExtra gss nlme
#' @import pracma reshape2 stringr
#' @importFrom grDevices as.graphicsAnnot col2rgb colorRampPalette dev.flush dev.hold dev.off pdf rgb tiff
#' @importFrom graphics boxplot lines par pie plot.new plot.window polygon text title
#' @importFrom stats aggregate cutree dist dnbinom hclust loess optim p.adjust pchisq pnbinom predict qnbinom
#' @importFrom utils capture.output read.csv setTxtProgressBar txtProgressBar
#' @export
#'
get_time_array<-function(raw.count.df){
  col.info.df<-as.data.frame(str_split(colnames(raw.count.df), "_", simplify = T))
  colnames(col.info.df)<-c("Project", "Time", "RepID")
  t.arr<-as.numeric(col.info.df$Time)
  return(t.arr)
}
#' Get the replicate array from count table
#'
#' It take the count table, grep the replicate element from column name. The column name of the count table must satisfy "Prj_Time_Rep1" format.
#'
#' @param raw.count.df, a count table.
#' @return a numeric array of each sample's replicate.
#' @export
#'
get_rep_array<-function(raw.count.df){
  col.info.df<-as.data.frame(str_split(colnames(raw.count.df), "_", simplify = T))
  colnames(col.info.df)<-c("Project", "Time", "RepID")
  rep.arr<-as.numeric(sub("Rep", "",col.info.df$RepID))
  return(rep.arr)
}


#' Convert a single gene's count row number into data frame with two columns.
#'
#' @param gene.row.info, a single row from count table.
#' @param gene.name, the gene name.
#' @param time.arr, the return value from get_time_array function.
#' @param rep.arr, the return value from get_rep_array function.
#' @return a dataframe object ordered by time and replicate id.
#' @export
#'
transform_single_gene_df<-function(gene.row.info, gene.name, time.arr, rep.arr){
  count.arr<-as.numeric(gene.row.info)
  data.trans = data.frame(Gene = gene.name, Count = count.arr, Time = time.arr, Rep = rep.arr)
  data.trans = data.trans[order(data.trans[,"Rep"], data.trans[,"Time"]),]
  return(data.trans)
}
#' Fit the non-baseline count data into a smoothed ANOVA model.
#'
#' @param data.trans, a data frame returned from transform_single_gene_df, without the baseline time.
#' @return a dataframe with all the fitted count from spline model.
#' @export
#'
fit_single_gene_spline<-function(data.trans){
  mod.fit <- gssanova(Count~Time, data = data.trans, family = "nbinomial")
  points <- unique(data.trans$Time)
  mod.pred <- predict(mod.fit, data.frame(Time = points), se = TRUE)
  # Prediction Data
  single.gene.spine.fit = data.frame(Gene = data.trans$Gene[1],Time = points, Count = mod.fit$nu/exp(mod.pred$fit))
  return(single.gene.spine.fit)
}

ConstNB <- function(theta.arr, count.arr, disp.var) {
  scaMu <- exp(theta.arr[1])
  if (scaMu < 10^(-10)) {
    scaMu <- 10^(-10)
  }
  log.like.var <- sum(dnbinom(
    count.arr,
    mu = scaMu,
    size = disp.var, log = TRUE))
  return(log.like.var)
}

ConstNB_comp <- compiler::cmpfun(ConstNB)

#' Fit the baseline count data into a constant negative binomial model.
#'
#' @param count.arr, a data frame returned from transform_single_gene_df, only with the baseline time.
#' @param disp.var, the dispersion value estimated from DESeq2.
#' @param MAXIT, 1000.
#' @param RELTOL, 10*-8
#' @param trace, 10
#' @return a list contain all the estimated value from NB model.
#' @export
#'
fit_single_gene_const <- function(
  # Inspired from ImpulseDE2
  count.arr, disp.var,
  MAXIT = 1000, RELTOL = 10^(-8), trace = 10) {
  param.guess.arr <- log(mean(count.arr, na.rm = TRUE) + 1)
  lsFit <- tryCatch({
    optim(par = param.guess.arr, fn = ConstNB_comp,
          count.arr = count.arr, disp.var = disp.var, method = "BFGS",
          control = list(maxit = MAXIT, reltol = RELTOL, fnscale = -1, trace = trace)
    )[c("par", "value", "convergence")]
  }, error = function(strErrorMsg) {
    print(paste0("ERROR: Fitting Const NB error!!!"))
    print(paste0("param.guess.arr ",
                 paste(param.guess.arr, collapse = " ")))
    print(paste0("count.arr ",
                 paste(count.arr, collapse = " ")))
    print(paste0("disp.var ",
                 paste(disp.var, collapse = " ")))
    print(paste0("MAXIT ", MAXIT))
    print(strErrorMsg)
    stop(strErrorMsg)
  })
  # Extract parameter estimates
  scaMu <- exp(lsFit$par[1])
  # Catch boundary of likelihood domain on mu space:
  if (scaMu < 10^(-10)) {
    scaMu <- 10^(-10)
  }
  scaNParamUsed <- 1
  return(list(scaMu = scaMu,
              disp.varParam = disp.var,
              scaLL = lsFit$value, scaConvergence = lsFit$convergence))
}

#' Calculate the significance of the dynamic signal for a single gene's single time expression.
#' Compared to baseline NB confidence interval.
#'
#' @param obs.val, estimated value from non-baseline model.
#' @param size, the size of negative binomial model.
#' @param mu, the estimated mean of constant negative binomial model
#' @return a numeric p-value for time t.
#' @export
#'

cal_p<-function(obs.val, size, mu){
  if(obs.val >= mu){
    p<-pnbinom(obs.val, size = size, mu = mu, lower.tail = F)
  } else{
    p<-pnbinom(obs.val, size = size, mu = mu, lower.tail = T)
  }
  return(p)
}
#' Calculate the significance of the dynamic signal for a single gene's over alll the time points.
#'
#' @param const.output, list returned from fit_single_gene_const function.
#' @param spline.output, dataframe returned from fit_single_gene_spline function.
#' @return p-value for a single point.
#' @export
#'
cal_time_p_single_gene<-function(const.output, spline.output){
  obs.arr<-spline.output$Count
  mu<-const.output$scaMu
  size<-const.output$disp.varParam
  p.arr<-sapply(obs.arr, cal_p, size, mu)
  p.arr[1]<-1 # To make sure baseline is not significant
  single_gene_time_p<-data.frame(Gene = spline.output$Gene, Time = spline.output$Time,
                                 Spine.Count = spline.output$Count, Const.Count.Mu = mu, Const.Disper = size,
                                 p.val = p.arr)
  return(single_gene_time_p)
}



# Due to the metaseqR package didn't support R 4.0, source code were attached here.
fisher.sum <- function(p,zero.sub=0.00001,na.rm=FALSE) {
  if(any(p>1, na.rm=TRUE)||any(p<0, na.rm=TRUE))
    stop("You provided bad p-values")
  stopifnot(zero.sub>=0 & zero.sub<=1 || length(zero.sub)!=1)
  p[p==0] <- zero.sub
  if (na.rm)
    p <- p[!is.na(p)]
  S = -2*sum(log(p))
  res <- data.frame(S=S,num.p=length(p))
  return(res)
}
fisher.method <- function(pvals,method=c("fisher"),p.corr=c("bonferroni","BH",
                                                            "none"),zero.sub=0.00001,na.rm=FALSE,mc.cores=NULL) {
  stopifnot(method %in% c("fisher"))
  stopifnot(p.corr %in% c("none","bonferroni","BH"))
  stopifnot(all(pvals>=0, na.rm=TRUE) & all(pvals<=1, na.rm=TRUE))
  stopifnot(zero.sub>=0 & zero.sub<=1 || length(zero.sub)!=1)
  if(is.null(dim(pvals)))
    stop("pvals must have a dim attribute")
  p.corr <- ifelse(length(p.corr)!=1, "BH", p.corr)
  ##substitute p-values of 0
  pvals[pvals == 0] <- zero.sub
  if(is.null(mc.cores)) {
    fisher.sums <- data.frame(do.call(rbind,apply(pvals,1,fisher.sum,
                                                  zero.sub=zero.sub,na.rm=na.rm)))
  } 
  else {
    fisher.sums <- parallel::mclapply(1:nrow(pvals), function(i) {
      fisher.sum(pvals[i,],zero.sub=zero.sub,na.rm=na.rm)
    }, mc.cores=mc.cores)
    fisher.sums <- data.frame(do.call(rbind,fisher.sums))
  }
  
  rownames(fisher.sums) <- rownames(pvals)
  fisher.sums$p.value <- 1-pchisq(fisher.sums$S,df=2*fisher.sums$num.p)
  fisher.sums$p.adj <- switch(p.corr,
                              bonferroni = p.adjust(fisher.sums$p.value,"bonferroni"),
                              BH = p.adjust(fisher.sums$p.value,"BH"),
                              none = fisher.sums$p.value
  )
  return(fisher.sums)
}

#' Combine multiple p-value using Fisher's p-value combination method.
#'
#' @param p.arr, a vector of p-values.
#' @return a numeric combined p-value.
#' @export
#'
combine_p_single_gene<-function(p.arr){
  p.arr<-p.arr[-1]
  p.combine<-fisher.method(matrix(p.arr, ncol=length(p.arr)),
                           method = c("fisher"),
                           p.corr="BH",
                           zero.sub = 1e-10,
                           na.rm = TRUE,
                           mc.cores=NULL)
  return(p.combine)
}

gene_pattern_assignment<-function(gene.df, gene.p.adj){
  if(gene.p.adj<=0.05){
    l.arr<-gene.df$t.p.val[2:(nrow(gene.df)-1)]
    idx<-which(l.arr<=0.05)+1
    idx<-c(1, idx, nrow(gene.df))
    pattern.list<-get_traj_pattern(gene.df$Fit.Count[idx])
    start.idx<-paste0(idx[pattern.list$start.idx], "_", collapse = "")
    end.idx<-paste0(idx[pattern.list$end.idx], "_", collapse = "")
    pattern.str<-paste0(pattern.list$pattern, "_", collapse = "")
    pattern.df<-data.frame(pattern = pattern.str, start.idx = start.idx, end.idx = end.idx)
  } else{
    pattern.df<-data.frame(pattern = "flat", start.idx = 0, end.idx = 0)
  }
  return(pattern.df)
}


get_traj_pattern<-function(t.arr){
  # Decide the pattern
  direction<-ifelse(t.arr[1]>=t.arr[2], "down", "up")
  segments.direction.push<-c(direction)
  segments.start.idx.push<-c(1)
  max.idx<-length(t.arr)-1
  start.idx<-c(1)
  end.idx<-numeric()
  for(i in 1:max.idx){
    curr.val<-t.arr[i]
    next.val<-t.arr[i+1]
    if(direction=="up"){
      if(curr.val>=next.val){
        segments.start.idx.push<-c(segments.start.idx.push, i)
        direction<-"down"
        segments.direction.push<-c(segments.direction.push, direction)
        end.idx<-c(end.idx, i)
        start.idx<-c(start.idx, i)
      }
    }else{
      if(curr.val<next.val){
        segments.start.idx.push<-c(segments.start.idx.push, i)
        direction<-"up"
        segments.direction.push<-c(segments.direction.push, direction)
        end.idx<-c(end.idx, i)
        start.idx<-c(start.idx, i)
      }
    }
  }
  end.idx<-c(end.idx, max.idx+1)
  pattern.list<-list(start.idx = start.idx, end.idx = end.idx, pattern = segments.direction.push)
  return(pattern.list)
}


return.center<-function(master.list, gene.arr){
  fit.count<-master.list$fitted.count
  fit.count$Symbol<-master.list$master.table$Symbol[match(fit.count$Gene, master.list$master.table$Gene)]
  fit.count<-fit.count %>% dplyr::filter(Symbol %in% gene.arr)
  df<-ddply(fit.count, .(Gene), function(df){
    base<-df$Fit.Count[1]
    df$Center.logFC<-log(df$Fit.Count/base, 2)
    return(df)
  })
  return(df)
}


#### Modify functions from package MetaLonDA R package functions

intervalArea.new = function(dd.1, dd.2){
  size = length(dd.1$Time)
  ar = numeric(size - 1)
  ## Calculate the absoulte and the sign of each interval area
  ar.abs = numeric(size - 1)
  ar.sign = numeric(size - 1)
  for(i in 1:(size - 1)){
    area.0 = trapz(dd.1$Time[i:(i+1)], dd.1$Count[i:(i+1)])
    area.1 = trapz(dd.2$Time[i:(i+1)], dd.2$Count[i:(i+1)])
    ar[i] = (abs(area.0) - abs(area.1)) / max(abs(area.0), abs(area.1))
    if(is.na(ar[i])){
      ar[i]=1
    }
    ar.abs[i] = abs(ar[i])
    ar.sign[i] = ar[i]/abs(ar[i])
  }
  area<-list(ar = ar, ar.abs = ar.abs, ar.sign = ar.sign)
  return(area)
}

#' CurveFitting function
#' 
#' To perform Local Polynomial Regression Fitting
#' 
#' @param perm.dat.1, data to fit the curve.
#' @param points, array of time values.
#' @return fitted count over time
#' @export
#'
curveFitting.new = function(perm.dat.1, points){
  group.1<-perm.dat.1 %>% dplyr::filter(Group == 1)
  group.2<-perm.dat.1 %>% dplyr::filter(Group == 2)
  mod.1 = loess(Count ~ Time, data = group.1)
  mod.2 = loess(Count ~ Time, data = group.2)
  est.1 = predict(mod.1, data.frame(Time = points), se = TRUE)
  est.2 = predict(mod.2, data.frame(Time = points), se = TRUE)
  dd.1 = data.frame(Time = points, Count = est.1$fit, Group = group.1.name)
  dd.2 = data.frame(Time = points, Count = est.2$fit, Group = group.2.name)
  output = list(dd.1, dd.2)
  return(output)
}


findSigInterval = function(adjusted.pvalue, threshold = 0.05, sign)
{
  sig = which(adjusted.pvalue < threshold)
  sign = sign[sig]
  padj = adjusted.pvalue[sig]
  start = numeric()
  end = numeric()
  p = numeric()
  dom = numeric()
  
  if(length(sig) == 0)
  {
    cat("No Significant Intevals Found \n")
  }
  else if(length(sig) == 1)
  {
    start = sig[1]
    end = sig [1]
    p = padj[1]
    dom = sign[1]
  }
  else
  {
    start = sig[1]
    
    if((sig[2] - sig[1]) != 1 | sign[2] != sign[1])
    {
      end = c(end, sig[1])
      dom = c(dom, sign[1])
      p = c(p, padj[1])
    }
    
    for(i in 2:length(sig))
    {
      if(i != length(sig))
      {
        if((sig[i] - sig[i-1]) > 1 | sign[i] != sign[i-1])
        {
          start= c(start, sig[i])
        }
        
        if((sig[i+1] - sig[i]) != 1 | sign[i+1] != sign[i])
        {
          end = c(end, sig[i])
          dom = c(dom, sign[i])
          p = c(p, mean(adjusted.pvalue[start[length(start)] : end[length(end)]]))
        }
      }
      else
      {
        if((sig[i]-sig[i-1]) > 1 | sign[i] != sign[i-1])
        {
          start= c(start, sig[i])
        }
        end= c(end, sig[i])
        dom = c(dom, sign[i])
        p = c(p, mean(adjusted.pvalue[start[length(start)] : end[length(end)]]))
      }
    }
  }
  
  return(list(start = start, end = end, pvalue = p, dominant = dom))
}


#' Permutation test
#' 
#' Function to perform permutation test.
#' 
#' @param perm.dat, temporal data frame to perform permutation.
#' @param n.perm, time to perm
#' @param points, array of time values
#' @param parall, run using multiple core
#' @return permutation test output data
#' @export
#'

permutation.new<-function(perm.dat, n.perm = 100, points, parall = FALSE){
  
  ## Start permutation
  cat("Start Permutation \n")
  pp = list() 
  perm = 0 # to be able to store the value
  n.subjects = length(unique(perm.dat$ID))
  #cat("# of Subjects = ", n.subjects, "\n")
  aggregate.df<-perm.dat
  
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
      cat("QQQQ")
      assign(paste("Model", j, sep = "_"), perm)
    }
  }, .parallel = parall, .progress = "text", .inform = TRUE,
  .paropts = list(.export=ls(.GlobalEnv),
                  .packages=.packages(all.available=T)))
  
  
  if(parall == TRUE) {
    stopCluster(cl)
  }
  
  pp[sapply(pp, is.null)] = NULL
  return(pp)
}


areaPermutation.new = function(perm)
{
  ar.list = list()
  list.len = length(perm)
  for (j in 1:list.len)
  {
    ar.list[[j]] = intervalArea.new(dd.1 = as.data.frame(perm[[j]][1]), dd.2 = as.data.frame(perm[[j]][2]))
  }
  
  return(ar.list)
}


