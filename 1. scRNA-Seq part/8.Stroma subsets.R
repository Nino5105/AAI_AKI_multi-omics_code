# 8.Stroma subsets cell

# 8.1 library packages

library(dplyr)
library(Seurat)
library(patchwork)
library(clustree)
library(ggsci)
library(ggplot2) 

Stroma <- readRDS("1.Data/5.3 SCT_Stroma.rds")
DimPlot(Stroma)
Stroma$cell_type <- Stroma@active.ident
table(Stroma$cell_type)
# Endo  Peri Fibro 
# 1566  1151  1056

DimPlot(Stroma,label = T,pt.size = 1,label.size = 5,group.by = "integrated_snn_res.0.8")

# 8.2 datasets dimension reduction and re-cluster

Stroma <- RunPCA(object = Stroma,verbose = FALSE)
ElbowPlot(Stroma,ndims = 20) # 选择主成分为15

Stroma <- FindNeighbors(Stroma,dim=1:15)
Stroma <- FindClusters(Stroma,resolution = 0.8)
Stroma <- RunUMAP(Stroma,reduction="pca", dims = 1:15)

DimPlot(Stroma,label = T,pt.size = 1,label.size = 5)
DimPlot(Stroma,label = T,pt.size = 1,label.size = 5,group.by = "cell_type")

# 8.3 marker gene expression

sub <- subset(Stroma,idents = c(0,1,7,8),invert =T)
DoHeatmap(sub, slot = "scale.data",features = c(
  "Ehd3","Plat","Slc14a1", # Endo
  "Mmp2","Emilin1","Sfrp2",
  "Plac8", # Fboroblasts
  "Acta2","Pdgfrb", # Myofibroblasts
  "Vim","S100a4", # Peri
  "Pecam1","Ptprb","Aqp1","Sema3g","Cldn5","" # GE
  "Pecam1","Ptprb","Kdr", # Descending vasa recta
  "Pecam1","Ptprb","Kdr","Vcam1","Plvap" # Ascending vasa recta
  "Scg2","Chgb","Stmn2" # Ganglia
)) + NoLegend()

# C9_marker <- FindMarkers(subset,ident.1 = 6,only.pos = T)

# 8.4 cell type annotation

# Stroma@active.ident <- Stroma$integrated_snn_res.0.8
head(Stroma@active.ident)
DimPlot(Stroma,label = T,group.by = "integrated_snn_res.0.8",label.size = 5,pt.size = 1)

new.cluster.ids <- c("GE",
                     "GE",
                     "Peri",
                     "Peri",
                     "Peri",
                     "Fibro",
                     "Peri",
                     "Endo",
                     "MyoFibro",
                     "Fibro",
                     "Fibro",
                     "Peri",
                     "Peri")

names(new.cluster.ids) <- levels(Stroma)
Stroma <- RenameIdents(Stroma, new.cluster.ids)

levels(Stroma) <- c("GE","Endo","Fibro","MyoFibro","Peri")

Stroma$cell_type_sub <- Stroma@active.ident

p1 <- DimPlot(Stroma,label = TRUE,reduction = "umap",label.size = 5,pt.size = 1) + scale_color_tron() + theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title = "Stroma cell")
p1

p2 <- DimPlot(Stroma,label = F,group.by = "type",label.size = 5,pt.size = 1,cols = c('Control' = '#00BFC4', 'Treatment' = '#F8766D')) + theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title =  paste0("Dataset type"))
p2

p3 <- DimPlot(Stroma,label = TRUE,reduction = "umap",label.size = 5,pt.size = 1,split.by = "type") + scale_color_tron()
p3

# DimPlot(Stroma,label = F,group.by = "sample",label.size = 5)+ theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title =  paste0("Sample id")) + scale_color_nejm()
# DimPlot(Stroma,label = T,group.by = "id",label.size = 5) + theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title =  paste0("Cluster id"))

# 8.5 draw heatmap plot

AverageExpression_value <-  AverageExpression(Stroma, assays = "SCT", features = c(
  # "Pecam1","Ptprb","Sema3g",
  # "S100a4","Plac8","Vim"
  "Pecam1","Ptprb","Kdr","Vcam1","Plvap", "Sema3g","Cldn5", # GE
  "Ehd3","Plat","Slc14a1", # Endo
  "Aqp1","Mmp2",
  "S100a4","Plac8", # Fboroblasts
  "Acta2","Pdgfrb", # Myofibroblasts
  "Vim","S100a4" # Peri
  ), return.seurat = FALSE, group.by = "ident")

pheatmap::pheatmap(AverageExpression_value$SCT,cluster_rows = F,cluster_cols = F,
                   border_color = "black",
                   angle_col = 45,
                   # gaps_row = c(4,7,19,21,28,30,32),#gaps_col = c(4,7,8),
                   cellwidth = 30,cellheight = 15,
                   treeheight_col = 4,scale = "row")


# 8.6 save metadata and Stroma.rds

Stroma$cell_type_sub <- Stroma@active.ident
write.csv(Stroma@meta.data,"1.Data/Stroma_meta_data.csv")
saveRDS(Stroma,"1.Data/5.3 Stroma_anno.rds")

# 8.7 count the proportions of each group and each sample in different subclusters

table(Stroma$type,Stroma@active.ident)

dfsam1 <- as.data.frame(table(Stroma$type,Stroma@active.ident))
dfsam1$Var1 <- factor(dfsam1$Var1,levels = rev(c("Control","Treatment")))
dfsam1$Var2 <- factor(dfsam1$Var2,levels = rev(levels(dfsam1$Var2)))

p4 <- ggplot(data = dfsam1,aes(x = Var2,y=Freq,fill = Var1)) +
  geom_bar(stat="identity",position = "fill",width = 0.7)+
  scale_y_continuous(labels = scales::percent)+
  scale_fill_manual(values = rev(c("#00BFC4","#F8766D")), name = "" )+
  theme_classic() + 
  labs(y = 'Fraction of different type',x="") +
  coord_flip()  +
  NoLegend()

dfsam2 <- as.data.frame(table(Stroma$sample,Stroma$cell_type_sub))
dfsam2$Var1 <- factor(dfsam2$Var1,levels = rev(c("KC1","KC2","KC3","KA1","KA2","KA3")))
dfsam2$Var2 <- factor(dfsam2$Var2,levels = rev(levels(dfsam2$Var2)))

p5 <- ggplot(data = dfsam2,aes(x = Var2,y=Freq,fill = Var1)) +
  geom_bar(stat="identity",position = "fill",width = 0.7)+
  scale_y_continuous(labels = scales::percent)+
  theme_classic() + 
  labs(y = 'Fraction of different sample',x="") +
  coord_flip() +
  scale_fill_nejm() + 
  NoLegend()


p4/p5 #3X8

# 8.8 Draw pipe chart

dfsam1 <- as.data.frame(table(Stroma$type,Stroma$cell_type_sub))

control <- dfsam1[dfsam1$Var1 == "Control",]
myLabel = as.vector(control$Var2)
myLabel = paste(myLabel, "(", round(control$Freq / sum(control$Freq) * 100, 2), "%)"
                , sep = "")  
p6 <- ggplot(control, aes(x = "", y = Freq,fill = Var2)) + 
  geom_bar(stat = "identity")+
  coord_polar(theta = "y") + 
  theme_bw() + 
  labs(x = "", y = "", title = "")+
  theme(axis.ticks = element_blank())+
  theme(axis.text.x = element_blank()) + 
  theme(panel.grid=element_blank()) +    
  theme(panel.border=element_blank())+
  scale_fill_d3(labels = myLabel) + 
  theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title =  paste0("Control group"))# 6x6
# 6x6

treatment <- dfsam1[dfsam1$Var1 == "Treatment",]
myLabel = as.vector(treatment$Var2)
myLabel = paste(myLabel, "(", round(treatment$Freq / sum(treatment$Freq) * 100, 2), "%)"
                , sep = "")  

p7 <- ggplot(treatment, aes(x = "", y = Freq,fill = Var2)) + 
  geom_bar(stat = "identity")+
  coord_polar(theta = "y") + 
  theme_bw() + 
  labs(x = "", y = "", title = "")+
  theme(axis.ticks = element_blank())+
  theme(axis.text.x = element_blank()) + 
  theme(panel.grid=element_blank()) +    
  theme(panel.border=element_blank())+
  scale_fill_d3(labels = myLabel)+
  theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title =  paste0("Treatment group"))# 6x6
# 6x6

p6/p7


# 8.9 Draw feature plot & vlnplot 

p1 <- FeaturePlot(Stroma, features = c("Pecam1","Ptprb","Sema3g",
                                        "S100a4","Plac8","Vim"),ncol = 3,slot = "scale.data",
                  cols = c("white","red"),label=F,pt.size = 0.2)+ 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)
p1

Myeloid <- readRDS("1.Data/5.3 Myeloid_anno.rds")
DimPlot(Myeloid)
Stroma <- subset(Myeloid,idents = c("Stroma M1","Stroma M2","Stroma Pro"))

p1 <- VlnPlot(Stroma, features = c("Il1b","Tnf","Tgfb1","Il10",""),
              ncol = 1,
              cols = c("#00BFC4","#F8766D"),
              # group.by = "type",
              split.by = "type",
              y.max = 4,
              # stack = T,
              assay = "SCT",
              slot = "scale.data",
              # flip = 90,
              pt.size = 0)
p1

p2 <- VlnPlot(Stroma, features = c("Arg1"),
              ncol = 1,
              cols = c("#00BFC4","#F8766D"),
              # group.by = "type",
              split.by = "type",
              y.max = 0.2,
              # stack = T,
              assay = "SCT",
              slot = "scale.data",
              # flip = 90,
              pt.size = 0)

p1+p2

# 8.10  Examine the marker gene expression

Stroma <- readRDS("1.Data/5.3 Stroma_anno.rds")
Stroma

p1 <- VlnPlot(Stroma, features = c("Icam1"),
              ncol = 1,
              cols = c("#00BFC4","#F8766D"),
              # group.by = "type",
              split.by = "type",
              y.max = 4,
              # stack = T,
              assay = "SCT",
              slot = "data",
              # flip = 90,
              pt.size = 0) + NoLegend()
p1

p2 <- VlnPlot(Stroma, features = c("Vcam1"),
              ncol = 1,
              cols = c("#00BFC4","#F8766D"),
              # group.by = "type",
              split.by = "type",
              y.max = 4,
              # stack = T,
              assay = "SCT",
              slot = "data",
              # flip = 90,
              pt.size = 0)+ NoLegend()
p2

p3 <- VlnPlot(Stroma, features = c("Acta2"),
              ncol = 1,
              cols = c("#00BFC4","#F8766D"),
              # group.by = "type",
              split.by = "type",
              y.max = 6,
              # stack = T,
              assay = "SCT",
              slot = "data",
              # flip = 90,
              pt.size = 0)+ NoLegend()
p3

p4 <- VlnPlot(subset(Stroma,Col1a1>0), features = c("Col1a1"),
              ncol = 1,
              cols = c("#00BFC4","#F8766D"),
              # group.by = "type",
              split.by = "type",
              # y.max = 0.5,
              # stack = T,
              assay = "SCT",
              slot = "data",
              # flip = 90,
              pt.size = 0)
p4


p5 <- VlnPlot(Stroma, features = c("Ccn2"),
              ncol = 1,
              cols = c("#00BFC4","#F8766D"),
              # group.by = "type",
              split.by = "type",
              # y.max = 0.5,
              # stack = T,
              assay = "SCT",
              slot = "data",
              # flip = 90,
              pt.size = 0)+ NoLegend()
p5


p6 <- VlnPlot(Stroma, features = c("Fn1"),
              ncol = 1,
              cols = c("#00BFC4","#F8766D"),
              # group.by = "type",
              split.by = "type",
              # y.max = 0.5,
              # stack = T,
              assay = "SCT",
              slot = "data",
              # flip = 90,
              pt.size = 0)+ NoLegend()
p6


p1 | p2 | p6 | p3 | p5 | p4

 # FeaturePlot(Stroma,features = c("Acta2"),pt.size = 1)
