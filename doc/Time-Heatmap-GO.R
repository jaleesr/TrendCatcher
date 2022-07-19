## ---- results='hide', message=FALSE-------------------------------------------
library(TrendCatcher)

## -----------------------------------------------------------------------------
demo.master.list.path<-system.file("extdata", "BrainMasterList_Symbol.rda", package = "TrendCatcher")
load(demo.master.list.path)

## ---- echo=FALSE, fig.height=15, fig.width=15---------------------------------
#time_heatmap<-draw_TimeHeatmap_GO(master.list = master.list, logFC.thres = 0, top.n = 10, dyn.gene.p.thres = 0.05, keyType = "SYMBOL", OrgDb = "org.Mm.eg.db", ont = "BP", term.width = 80, GO.enrich.p = 0.05, figure.title = "TimeHeatmap")  


## -----------------------------------------------------------------------------
# To save time, directely load from extdata
demo.time.heatmap.path<-system.file("extdata", "Brain_TimeHeatmap.rda", package = "TrendCatcher")
load(demo.time.heatmap.path)
names(time_heatmap)

## ---- fig.width=9, fig.height=10----------------------------------------------
require("ComplexHeatmap")
print(time_heatmap$time.heatmap)

## -----------------------------------------------------------------------------
head(time_heatmap$merge.df[,1:5])

## -----------------------------------------------------------------------------
head(time_heatmap$GO.df[,1:5])

## ---- fig.width=8, fig.height=5-----------------------------------------------
go.terms<-unique(time_heatmap$GO.df$Description)[1:5]
time_heatmap_selGO<-draw_TimeHeatmap_selGO(time_heatmap = time_heatmap, sel.go = go.terms, master.list = master.list, GO.perc.thres = 0, nDDEG.thres = 0, save.tiff.path = NA)

## ---- fig.height=10, fig.width=7----------------------------------------------
go.terms<-c("response to lipopolysaccharide",
            "response to interferon-beta",
            "cytokine-mediated signaling pathway",
            "response to interferon-gamma",
            "response to virus",
            "leukocyte migration",
            "mitotic nuclear division",
            "regulation of vasculature development",
            "extracellular structure organization",
            "regulation of epithelial cell proliferation")
gene.GO.df<-draw_GOHeatmap(master.list = master.list, time.window = "0h-6h", 
                           go.terms = go.terms,  merge.df = time_heatmap$merge.df, 
                           logFC.thres = 5)

## -----------------------------------------------------------------------------
head(gene.GO.df$GOheatmapDat)

