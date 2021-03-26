run_TrendCatcherSingleGene<-function(raw.df, i, gene.dispersion, baseline.t){
  # Get gene count
  gene.row.info<-raw.df[i,]
  # Get gene name
  gene.name<-rownames(raw.df)[i]
  # Get gene disp
  gene.disp<-gene.dispersion[i]
  # Get time array
  time.arr<-get_time_array(raw.count.df = raw.df)
  # Get rep array
  rep.arr<-get_rep_array(raw.count.df = raw.df)

  # Transform gene count array to gene data frame
  data.trans<-transform_single_gene_df(gene.row.info = gene.row.info,
                                       gene.name = gene.name,
                                       time.arr = time.arr,
                                       rep.arr = rep.arr)
  # Get non-baseline count table, fit to spline model
  data.trans.sig<-data.trans[-which(data.trans$Time == baseline.t),]
  spline.output<-fit_single_gene_spline(data.trans = data.trans.sig)
  # Get baseline count table, fit to const NB model
  base.count<-data.trans$Count[which(data.trans$Time == baseline.t)]
  #const.output<-fit_single_gene_const(count.arr = base.count, disp.var = gene.disp)
  invisible(capture.output(const.output<-fit_single_gene_const(count.arr = base.count, disp.var = gene.disp)))
  # Combine the fitted data into fitted count table
  spline.output<-rbind(data.frame(Gene = gene.name, Time = baseline.t, Count = const.output$scaMu), spline.output)
  gene.p<-cal_time_p_single_gene(const.output = const.output, spline.output = spline.output)
  gene.i.fit.df<-data.frame(Gene = spline.output$Gene, Time = spline.output$Time, Fit.Count = spline.output$Count,
             mu = const.output$scaMu, disp = const.output$disp.varParam, t.p.val = gene.p$p.val)
  return(gene.i.fit.df)
}

