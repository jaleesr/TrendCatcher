my.pie<-function (x, labels = names(x), edges = 200, radius = 0.8, clockwise = FALSE,
                  init.angle = if (clockwise) 90 else 0, density = NULL, angle = 45,
                  col = NULL, border = NULL, lty = NULL, main = NULL, cex,...)
{
  if (!is.numeric(x) || any(is.na(x) | x < 0))
    stop("'x' values must be positive.")
  if (is.null(labels))
    labels <- as.character(seq_along(x))
  else labels <- as.graphicsAnnot(labels)
  x <- c(0, cumsum(x)/sum(x))
  dx <- diff(x)
  nx <- length(dx)
  plot.new()
  pin <- par("pin")
  xlim <- ylim <- c(-1, 1)
  if (pin[1L] > pin[2L])
    xlim <- (pin[1L]/pin[2L]) * xlim
  else ylim <- (pin[2L]/pin[1L]) * ylim
  dev.hold()
  on.exit(dev.flush())
  plot.window(xlim, ylim, "", asp = 1)
  if (is.null(col))
    col <- if (is.null(density))
      c("white", "lightblue", "mistyrose", "lightcyan",
        "lavender", "cornsilk")
  else par("fg")
  if (!is.null(col))
    col <- rep_len(col, nx)
  if (!is.null(border))
    border <- rep_len(border, nx)
  if (!is.null(lty))
    lty <- rep_len(lty, nx)
  angle <- rep(angle, nx)
  if (!is.null(density))
    density <- rep_len(density, nx)
  twopi <- if (clockwise)
    -2 * pi
  else 2 * pi
  t2xy <- function(t) {
    t2p <- twopi * t + init.angle * pi/180
    list(x = radius * cos(t2p), y = radius * sin(t2p), an=t2p)
  }
  for (i in 1L:nx) {
    n <- max(2, floor(edges * dx[i]))
    P <- t2xy(seq.int(x[i], x[i + 1], length.out = n))
    polygon(c(P$x, 0), c(P$y, 0), density = density[i], angle = angle[i],
            border = border[i], col = col[i], lty = lty[i])
    P <- t2xy(mean(x[i + 0:1]))
    lab <- as.character(labels[i])
    if (!is.na(lab) && nzchar(lab)) {
      lines(c(1, 1.05) * P$x, c(1, 1.05) * P$y)
      text(1.05 * P$x, 1.05 * P$y, labels[i], xpd = TRUE,
           #srt = ifelse(P$x < 0, P$an/pi*180+180, P$an/pi*180),
           adj = ifelse(P$x < 0, 1, 0), cex = cex, ...)
    }
  }
  title(main = main, ...)
  invisible(NULL)
}


#' Draw Grouped DDEGs Main-type and Sub-type Composition in Hierarchical Pie Chart
#'
#' Group all DDEGs based on their main-type and sub-type trajectory patterns and plot their composition in a 
#' hierachical pie chart. Inner pie chart represents the main-type trajectory pattern composition. The outer pie chart
#' represents sub-type trajectory pattern composition. 
#' 
#' 
#' @param master.list, a list object. The output from run_TrendCatcher function, contains master.table element.
#' @param fig.title, a string. The main title of the figure. By default if "". 
#' @param inner.radius, a numeric variable. The inner pie chart radius size. By default is 0.7.
#' @param cex.out, a numeric variable. The text size of label of outer pie chart. By default is 1.
#' @param cex.in, a numeric variable. The text size of the label of inner pie chart. By default is 1.
#' @export
#'

draw_TrajClusterPie<-function(master.list, fig.title = "", inner.radius = 0.7, cex.out = 1, cex.in = 1){
  
  if(is.na(master.list[[1]])){stop("master.list object is NA!!!")}
  
  ### 1. Remove the flat trajectory
  pattern.df<-master.list$master.table %>% dplyr::filter(pattern!="flat")

  ### 2. Group trajectory patterns, based on master type, and sub-type (time dependent)
  tab<-table(pattern.df$pattern, pattern.df$start.t)
  tab<-as.data.frame(tab)
  tab<-tab %>% dplyr::filter(Freq !=0 )
  tab<-tab %>% dplyr::arrange(desc(Freq))
  pattern.str<-rep("", nrow(tab))
  for(i in 1:nrow(tab)){
    tab.info<-tab[i,]
    pattern.i<-tab.info$Var1
    start.t.i<-tab.info$Var2
    title.str<-as.character(pattern.i)
    title.str<-str_split(title.str, "_", simplify = T)
    title.str<-title.str[-length(title.str)]
    act.str<-as.character(start.t.i)
    act.str<-str_split(act.str, "_", simplify = T)
    act.str<-act.str[-length(act.str)]
    pattern.str[i]<-paste(act.str, title.str, collapse = " ")
  }
  tab$SubPattern<-pattern.str
  var1.df<-str_split(tab$Var1, "_", simplify = T)
  tab$MasterPattern<-apply(var1.df, 1, function(x){x<-x[x!=""]; paste(x, collapse = " ")})
  ### 3.Prepare data for pie-chart
  data.out<-tab
  data.out$Perc<-data.out$Freq*100/sum(data.out$Freq)
  data.out<-data.out %>% arrange(MasterPattern)
  data.in<-aggregate(Perc~MasterPattern, data.out, sum)
  colnames(data.in)[2]<-"Share"
  data.out<-merge(data.in,data.out)
  data.out<-data.out %>% group_by(MasterPattern, SubPattern) %>% arrange(desc(Share))
  data.in<-data.in %>% arrange(desc(Share))
  data.out<-as.data.frame(data.out)
  data.in<-as.data.frame(data.in)

  ### 4.Create a new pie function to save overwriting original
  newpie <- pie
  newlbs <- quote(if (!is.na(lab) && nzchar(lab)) {
    #lines(c(1, 1.05) * P$x, c(1, 1.05) * P$y)
    text(0.4 * P$x, 0.6 * P$y, labels[i], xpd = TRUE,
         adj = ifelse(P$x < 0, 1, 0), ...)
  })
  # add in the new lines of code - trial and error found the right position
  body(newpie)[[22]][[4]][[7]] <- newlbs

  subcolors <- function(.dta,main,mainCol){
    tmp_dta = cbind(.dta,1,'col')
    tmp1 = unique(.dta[[main]])
    for (i in 1:length(tmp1)){
      tmp_dta$"col"[.dta[[main]] == tmp1[i]] = mainCol[i]
    }
    u <- unlist(by(tmp_dta$"1",tmp_dta[[main]],cumsum))
    n <- dim(.dta)[1]
    subcol=rep(rgb(0,0,0),n);
    for(i in 1:n){
      t1 = col2rgb(tmp_dta$col[i])/256
      subcol[i]=rgb(t1[1],t1[2],t1[3],alpha = 1/(1+u[i]))
    }
    return(subcol);
  }

  n <- length(unique(data.out$MasterPattern))
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  mainCol=col_vector[10:(n+10)]

  data.out$label<-ifelse(data.out$Perc<0.5, "", paste0(data.out$SubPattern,"(",round(data.out$Perc,1),"%", ")"))
  data.in$label<-ifelse(data.in$Share<5, "", paste0(as.character(data.in$MasterPattern), "\n", round(data.in$Share,1),"%"))
  # Draw the pie chart
  par(mai= c(1,2,1,2))
  my.pie(data.out$Perc, col=subcolors(data.out,"MasterPattern",mainCol),
         data.out$label, radius = 1, border = 1, main = fig.title, cex = cex.out)
  par(new=T)
  newpie(data.in$Share, col=mainCol, border = F,labels = data.in$label,radius = inner.radius, cex=cex.in)
}
