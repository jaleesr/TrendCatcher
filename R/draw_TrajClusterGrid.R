#' Draw Grouped DDEGs Trajectories in Grid Plot
#'
#' Group all DDEGs based on their sub-type trajectory patterns and plot their trajectories together, 
#' then layout all sub-type trajectory patterns which contains more than N genes in a grid plot. 
#' Each individual sub-grid plot is titled with sub-type trajectory pattern and number of genes included.
#' X-axis is the time, Y-aixs is log2 transformed fitted count trajectory.
#' 
#' @param master.list, a list object. The output from run_TrendCatcher function, contains master.table element.
#' @param min.traj.n, an integer variable. The minimum number of genes from the same sub-type trajectory. By default if 10.
#' @param save.as.PDF, a string. The absolute file path to save the figure as PDF file if needed. If set to NA, will print
#' the figure instead of saving it as PDF file.
#' @param pdf.width, a numeric variable. The PDF file width size. By default is 10.
#' @param pdf.height, a numeric variable. The PDF file height size. By default is 10.
#' @export
#'
draw_TrajClusterGrid<-function(master.list, min.traj.n = 10, save.as.PDF = NA, pdf.width =10, pdf.height = 10){
  
  if(is.na(master.list[[1]])){stop("master.list object is NA!!!")}
  
  ### 1. Remove the flat trajectory
  pattern.df<-master.list$master.table %>% dplyr::filter(pattern!="flat")

  ### 2. Group trajectory patterns, based on master type, and sub-type (time dependent)
  tab<-table(pattern.df$pattern_str)
  tab<-as.data.frame(tab)
  tab<-tab %>% dplyr::filter(Freq !=0 )
  tab<-tab %>% arrange(desc(Freq))
  colnames(tab)<-c("SubPattern", "Freq")

  # Draw the grid plot for trajectories
  tab.filter<-tab %>% dplyr::filter(Freq>min.traj.n)
  tab.filter<-tab.filter %>% dplyr::arrange(desc(Freq))
  nb.cols <- nrow(tab.filter)
  mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

  g<-list()
  for(i in 1:nrow(tab.filter)){
    pattern.i<-tab.filter[i,]$SubPattern
    sub<-pattern.df %>% dplyr::filter(pattern_str == pattern.i)
    gene.sel<-sub$Gene
    title.str<-as.character(pattern.i)
    title.str<-str_split(title.str, "_", simplify = T)
    title.str<-paste0(title.str[-length(title.str)], collapse = " ")
    n.traj.str<-paste0("\n ", tab.filter[i,]$Freq, " Genes", collapse = "")
    title.str<-paste0(title.str, n.traj.str, collapse = "")
    g[[i]]<-master.list$fitted.count %>% dplyr::filter(Gene %in% gene.sel) %>%
      ggplot(aes(x = Time, y = Fit.Count, group = Gene)) +
      geom_line(color = mycolors[i], alpha = 0.5) + scale_y_continuous(name = "log2 Count",
                                                                       trans = "log2", labels = scales::number_format(accuracy = 0.01)) +
      ggtitle(paste(title.str, collapse = " "))
  }

  n <- length(g)
  nCol <- floor(sqrt(n))
  if(!is.na(save.as.PDF)){
    pdf(save.as.PDF, height = pdf.height, width = pdf.width)
    p<-do.call("grid.arrange", c(g, ncol=nCol))
    dev.off()
  } else{
    p<-do.call("grid.arrange", c(g, ncol=nCol))
  }
}
