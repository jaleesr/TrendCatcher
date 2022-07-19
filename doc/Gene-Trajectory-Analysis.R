## ---- results='hide', message=FALSE-------------------------------------------
library(TrendCatcher)

## -----------------------------------------------------------------------------
demo.master.list.path<-system.file("extdata", "BrainMasterList.rda", package = "TrendCatcher")
load(demo.master.list.path)

## -----------------------------------------------------------------------------
### In case bioMart has connection issue, you can load from example data

#gene.symbol.df<-get_GeneEnsembl2Symbol(ensemble.arr = master.list$master.table$Gene)
#master.table.new<-cbind(master.list$master.table, gene.symbol.df[match(master.list$master.table$Gene, gene.symbol.df$Gene), c("Symbol", "description")])
#master.list$master.table<-master.table.new
#head(master.list$master.table)

demo.master.list.path<-system.file("extdata", "BrainMasterList_Symbol.rda", package = "TrendCatcher")
load(demo.master.list.path)
head(master.list$master.table)

## -----------------------------------------------------------------------------
### ONLY use this command if CSV file is using GENE SYMBOL as row name!!!!!!
#master.list$master.table$Symbol<-master.list$master.table$Gene 

## ---- fig.width=10, fig.height=6----------------------------------------------
gene.symbol.arr<-unique(master.list$master.table$Symbol)[1:6]
p<-draw_GeneTraj(master.list = master.list, gene.symbol.arr = gene.symbol.arr, ncol = 3, nrow = 2)
p

## ---- fig.width=10, fig.height=10---------------------------------------------
draw_TrajClusterGrid(master.list = master.list, min.traj.n = 10)

## ---- fig.height=7------------------------------------------------------------
#par(mar=c(1,1,1,1))
#draw_TrajClusterPie(master.list = master.list,inner.radius = 0.7, cex.out = 1, cex.in = 1, fig.title = "Hierarchical Pie Chart")

