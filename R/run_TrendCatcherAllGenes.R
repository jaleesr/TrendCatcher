
run_TrendCatcherAllGenes<-function(raw.df, gene.dispersion, cores, baseline.t, show.verbose){
  if(show.verbose == FALSE){
    # Make sure cores doesn't exceed limit
    max.core<-parallel::detectCores()
    if(is.na(cores) | cores>=max.core){
      cores=parallel::detectCores()
      cl <- parallel::makeCluster(cores[1]-1)
    } else{
      cl <- parallel::makeCluster(cores[1])
    }
    doSNOW::registerDoSNOW(cl)
    
    # Define progress bar
    pb<-txtProgressBar(0,nrow(raw.df),style=3)
    progress<-function(n){
      setTxtProgressBar(pb,n)
    }
    opts<-list(progress=progress)
    
    
    # Run para
    system.time({
      finalMatrix <- foreach(i=1:nrow(raw.df), .combine=rbind, .packages=c('stringr','gss'), .verbose = show.verbose,
                             .export = c("run_TrendCatcherSingleGene", "get_time_array", "get_rep_array", "transform_single_gene_df", "fit_single_gene_spline",
                                         "ConstNB", "ConstNB_comp", "fit_single_gene_const", "cal_p", "cal_time_p_single_gene"),
                             .options.snow=opts
      ) %dopar% {
        tempMatrix = run_TrendCatcherSingleGene(raw.df = raw.df, i = i,
                                                gene.dispersion = gene.dispersion, baseline.t = baseline.t)
        tempMatrix
      }
    })
    
    #stop cluster
    parallel::stopCluster(cl)
    message("")
    return(finalMatrix)
  } else{
    message("You are setting verbose to TRUE, this function only use 1 CPU, so it may be slower than set it to FALSE.")
    finalMatrix.list<-list()
    for(i in 1:nrow(raw.df)){
      gene.name<-rownames(raw.df)[i]
      print(gene.name)
      finalMatrix.list[[i]] = run_TrendCatcherSingleGene(raw.df = raw.df, i = i,
                                                         gene.dispersion = gene.dispersion, baseline.t = baseline.t)
    }
    finalMatrix<-do.call(rbind, finalMatrix.list)
    return(finalMatrix)
  }
}
