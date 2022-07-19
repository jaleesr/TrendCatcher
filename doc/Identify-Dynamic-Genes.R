## ---- results='hide', message=FALSE, echo=FALSE-------------------------------
library(TrendCatcher)

## -----------------------------------------------------------------------------
example.file.path<-system.file("extdata", "Brain_DemoCountTable.csv", package = "TrendCatcher")
tb<-read.csv(example.file.path, row.names = 1)
head(tb)

## ---- results='hide', eval=FALSE----------------------------------------------
#  example.file.path<-system.file("extdata", "Brain_DemoCountTable.csv", package = "TrendCatcher")
#  
#  master.list<-run_TrendCatcher(count.table.path = example.file.path,
#  baseline.t = 0,
#  time.unit = "h",
#  min.low.count = 1,
#  para.core.n = NA,
#  dyn.p.thres = 0.05)

## -----------------------------------------------------------------------------
demo.master.list.path<-system.file("extdata", "BrainMasterList.rda", package = "TrendCatcher")
load(demo.master.list.path)

## -----------------------------------------------------------------------------
names(master.list)

## -----------------------------------------------------------------------------
print(c(master.list$time.unit, master.list$baseline.t))

## -----------------------------------------------------------------------------
master.list$t.arr

## -----------------------------------------------------------------------------
master.list$Project.name

## -----------------------------------------------------------------------------
head(master.list$raw.df)

## -----------------------------------------------------------------------------
head(master.list$fitted.count)

## -----------------------------------------------------------------------------
head(master.list$master.table)

