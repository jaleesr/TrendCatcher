## ---- results='hide', message=FALSE-------------------------------------------
library(TrendCatcher)

## -----------------------------------------------------------------------------
severe.path<-system.file("extdata", "MasterListSevere.rda", package = "TrendCatcher")
load(severe.path)
moderate.path<-system.file("extdata", "MasterListModerate.rda", package = "TrendCatcher")
load(moderate.path)
ht.path<-system.file("extdata", "htSevere.rda", package = "TrendCatcher")
load(ht.path)

## -----------------------------------------------------------------------------
#head(ht.severe$GO.df)
head(ht.severe$GO.df)

## -----------------------------------------------------------------------------
ht.severe$merge.df %>% filter(Description == "neutrophil activation")

## ---- fig.height=10, fig.width=10---------------------------------------------
g<-draw_CurveComp(master.list.1 = master.list.severe, master.list.2 = master.list.moderate, ht.1 = ht.severe, pathway = "neutrophil activation",group.1.name = "severe", group.2.name = "moderate")
print(g)


## ---- echo=TRUE, results='hide'-----------------------------------------------
perm_output<-draw_CurveComp_Perm(master.list.1 = master.list.severe, 
                                 master.list.2 = master.list.moderate, 
                                 ht.1 = ht.severe, 
                                 pathway = "neutrophil activation", 
                                 group.1.name = "severe", 
                                 group.2.name = "moderate", 
                                 n.perm = 100, 
                                 parall = FALSE, 
                                 pvalue.threshold = 0.05)

## -----------------------------------------------------------------------------
names(perm_output)

## ----fig.height=10, fig.width=10----------------------------------------------
perm_output$plot

