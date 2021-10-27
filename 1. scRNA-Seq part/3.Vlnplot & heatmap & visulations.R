
# 3.Vlnplot & Heatmap & Featureplot visulations

# 3.1 library packages
library(dplyr)
library(Seurat)
library(patchwork)
library(clustree)
library(ggsci)
library(ggplot2) 
options(future.globals.maxSize = 20000 * 1024^2)


scc_integrated <-readRDS("0.Data/SCT_combine.sample.anno.rds")
DimPlot(scc_integrated)

# 3.2 Draw heatamap plot

pdf("Heatmap.pdf",width = 30,height = 15)

DoHeatmap(scc_integrated, slot = "data",features = c(
  "Slc27a2","Lrp2","Slc22a8","Slc5a2","Slc5a12","Fxyd2","Slc17a3", # proximal tubule(PT)
  "Atp11a", "Slc13a3","Slc34a1","Gpx3"# Proximal convoluted tubule cell(PTC)
  "Aqp1","Bst1" # Descending loop of Henle (DLH)
  "Slc12a1","Umod","Cldn8","Krt18","Krt8",  # Ascending loop of Henle(ALH)
  "Slc12a3",  # Distal convoluted tubule(DCT)
  "Atp6v0d2","Atp6v1g3","Slc4a1","Aqp6","Slc26a4","Hmx2" # Collecting duct intercalated cell(CD-IC)
  "Aqp2","Hsd11b2", #  Collecting duct principal / epithelial cell (CD-PC)
  "Rhbg","Insrr","Stmn1", # Collecting duct transitional cell (CD-TC)
  "Cdca3","Mki67" # Novel cell
  "Nrp1","Kdr","Ehd3","Plat","Vim","S100a4","Aqp1","Bst1", # Endothelial(Endo)
  "Nphs1", "Nphs2", # Podocyte (Podo)
  "Vim","S100a4" # Pericytes and vascular smooth muscle (Peri)
  "Plac8", # Fibroblast (Fibro)
  "C1qa","C1qb", # Macrophage (Macro)
  "Cd79a", "Cd79b", # B lymphocyte (B lymph)
  "Cxcr6","Ltb","Il7r","Cd3d","Cd3e","Ifng",  # T lymphocyte (T lymph)
  "Gzma","Nkg7","Gnly", # Natural killer cell (NK)
  "Lyz","Cd14", # Monocytes (Mono)
  "S100a8","S100a9" # Neutrophil (Neutro)
)) + NoLegend()

dev.off()

# 3.3 Draw  Vlnplot

modify_vlnplot <- function(obj,
                          feature,
                          pt.size = 0,
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "mm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  +
    xlab("") + ylab(feature) + ggtitle("") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(),
          axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
          plot.margin = plot.margin )
  return(p)
}

StackedVlnPlot<- function(obj, features,
                          pt.size = 0,
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

my36colors <-c(
  '#E95C59', '#53A85F','#F1BB72','#F3B1A0','#D6E7A3','#57C3F3',
  '#E63863','#E4C755','#E59CC4','#AB3282','#23452F','#BD956A',
  '#8C549C','#58A4C3',"#00BFC4",'#E39A35')

feature = c(
  "Slc27a2","Lrp2", # PT
  "Slc34a1","Gpx3", # PT-C
  "Aqp1","Bst1", # DLH
  "Slc12a1","Umod", # ALH
  "Slc12a3", # DCT
  "Atp6v0d2","Atp6v1g3", # CD-IC
  "Aqp2","Hsd11b2", # CD-PC
  "Mki67","Stmn1", # Novel
  "Nrp1","Kdr",  # Endo
  "Lamb2","Wt1", # Podo
  "Vim","S100a4", # Peri
  "S100a4","Plac8", #  Fibro
  "C1qa","C1qb", # Macro
  "S100a8","S100a9", # Neutro
  "Cd79a","Cd79b", # B lymph
  "Ltb","Il7r", # T lymph
  "Gzma","Nkg7" # NK
)

StackedVlnPlot(scc_integrated, c(feature), pt.size=0, cols=my36colors)

# 3.3 Draw  Featureplot

# 3.3.1 PT
p1 <- FeaturePlot(scc_integrated, features = c("Slc27a2"),cols = c("white",'#E95C59'),label=F)+ 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)
p1

p2 <- FeaturePlot(scc_integrated, features = c("Lrp2"),cols = c("white",'#E95C59'),label=F)+ 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)
p2

p1 + p2

# 3.3.2 DLH_feature_plot

p1 <- FeaturePlot(scc_integrated, features = c("Aqp1"),cols = c("white",'#53A85F'),label=F)+ 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)
p1

p2 <- FeaturePlot(scc_integrated, features = c("Bst1"),cols = c("white",'#53A85F'),label=F,slot = "scale.data")+ 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)
p2

p1 + p2

# 3.3.3 ALH_feature_plot

p1 <- FeaturePlot(scc_integrated, features = c("Slc12a1"),cols = c("white",'#F1BB72'),label=F)+ 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)
p1

p2 <- FeaturePlot(scc_integrated, features = c("Umod"),cols = c("white",'#F1BB72'),label=F,slot = "scale.data")+ 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)
p2

p1 + p2

# 3.3.4 CD-IC_feature_plot

p1 <- FeaturePlot(scc_integrated, features = c("Atp6v0d2"),cols = c("white",'#D6E7A3'),label=F)+ 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)
p1

p2 <- FeaturePlot(scc_integrated, features = c("Atp6v1g3"),cols = c("white",'#D6E7A3'),label=F,slot = "scale.data")+ 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)
p2

p1 + p2

# 3.3.5 CD-PC_feature_plot

p1 <- FeaturePlot(scc_integrated, features = c("Aqp2"),cols = c("white",'#57C3F3'),label=F)+ 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)
p1

p2 <- FeaturePlot(scc_integrated, features = c("Hsd11b2"),cols = c("white",'#57C3F3'),label=F,slot = "scale.data")+ 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)
p2

p1 + p2

# 3.3.6 Novel_feature_plot

p1 <- FeaturePlot(scc_integrated, features = c("Mki67"),cols = c("white",'#E63863'),label=F)+ 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)
p1

p2 <- FeaturePlot(scc_integrated, features = c("Stmn1"),cols = c("white",'#E63863'),label=F,slot = "scale.data")+ 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)
p2

p1 + p2

# 3.3.7 Endo_feature_plot

p1 <- FeaturePlot(scc_integrated, features = c("Nrp1"),cols = c("white",'#E4C755'),label=F)+ 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)

p2 <- FeaturePlot(scc_integrated, features = c("Kdr"),cols = c("white",'#E4C755'),label=F,slot = "scale.data")+ 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)

p1 + p2

# 3.3.8 Podo_feature_plot

p1 <- FeaturePlot(scc_integrated, features = c("Lamb2"),cols = c("white",'#E59CC4'),label=F)+ 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)

p2 <- FeaturePlot(scc_integrated, features = c("Wt1"),cols = c("white",'#E59CC4'),label=F,slot = "scale.data")+ 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)

p1 + p2

# 3.3.9 Peri_feature_plot

p1 <- FeaturePlot(scc_integrated, features = c("Vim"),cols = c("white",'#AB3282'),label=F)+ 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)

p2 <- FeaturePlot(scc_integrated, features = c("S100a4"),cols = c("white",'#AB3282'),label=F,slot = "scale.data")+ 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)

p1 + p2

# 3.3.10 Fibro_feature_plot

p1 <- FeaturePlot(scc_integrated, features = c("S100a4"),cols = c("white",'#23452F'),label=F)+ 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)

p2 <- FeaturePlot(scc_integrated, features = c("Plac8"),cols = c("white",'#23452F'),label=F,slot = "scale.data")+ 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)

p1 + p2

# 3.3.11 Macro_feature_plot

p1 <- FeaturePlot(scc_integrated, features = c("C1qa"),cols = c("white",'#BD956A'),label=F)+ 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)

p2 <- FeaturePlot(scc_integrated, features = c("C1qb"),cols = c("white",'#BD956A'),label=F,slot = "scale.data")+ 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)

p1 + p2

# 3.3.12 Neutro_feature_plot

p1 <- FeaturePlot(scc_integrated, features = c("S100a8"),cols = c("white",'#8C549C'),label=F)+ 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)

p2 <- FeaturePlot(scc_integrated, features = c("S100a9"),cols = c("white",'#8C549C'),label=F,slot = "scale.data")+ 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)

p1 + p2

# 3.3.13 B lymph_feature_plot

p1 <- FeaturePlot(scc_integrated, features = c("Cd79a"),cols = c("white",'#58A4C3'),label=F)+ 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)

p2 <- FeaturePlot(scc_integrated, features = c("Cd79b"),cols = c("white",'#58A4C3'),label=F,slot = "scale.data")+ 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)

p1 + p2

# 3.3.14 T lymph & NK_feature_plot

p1 <- FeaturePlot(scc_integrated, features = c("Il7r"),cols = c("white",'#00BFC4'),label=F,slot = "scale.data")+ 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)

p2 <- FeaturePlot(scc_integrated, features = c("Nkg7"),cols = c("white",'#E39A35'),label=F,slot = "scale.data")+ 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)

p1 + p2

# 3.3.15 DCT_feature_plot

p1 <- FeaturePlot(scc_integrated, features = c("Slc12a3"),cols = c("white",'#F3B1A0'),label=F)+ 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)
p1




