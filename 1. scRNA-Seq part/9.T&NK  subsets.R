# 9.T lymphocyte & NK subsets cell

# 9.1 library packages#

library(dplyr)
library(Seurat)
library(patchwork)
library(clustree)
library(ggsci)
library(ggplot2) 
library(clusterProfiler)
library(org.Mm.eg.db)
library(monocle)

T_NK_lymph <- readRDS("1.Data/5.2 SCT_T_NK_lymph.rds")
DimPlot(T_NK_lymph)

# 9.2 datasets dimension reduction and re-cluster

T_NK_lymph <- RunPCA(object = T_NK_lymph,verbose = FALSE)
ElbowPlot(T_NK_lymph,ndims = 20) # 选择主成分为15

T_NK_lymph <- FindNeighbors(T_NK_lymph,dim=1:20)
T_NK_lymph <- FindClusters(T_NK_lymph,resolution = 0.8)
T_NK_lymph <- RunUMAP (T_NK_lymph,reduction="pca", dims = 1:20)
DimPlot(T_NK_lymph,label = T,pt.size = 1)

T_NK_lymph@assays

# 9.3 marker gene expression

DimPlot(T_lymph,label = T,pt.size = 1,label.size = 5)

DoHeatmap(T_NK_lymph, slot = "scale.data",features = c(
  "Cd3d","Il7r","Cd4","Cd8a", # T cell markers
  "Ccr7","Sell","Lef1","Tcf7", # T Naive markers
  "Il2","Il4","Il6","Il17a", # T Effector (,"Il9","Il10","Il22")
  "Ifng","Fasl","Nkg7","Gzma","Gzmk",# T Cytotoxic  markers
  "S100a4","Cxcr3", "Cd40lg","Cd69", # T Memory markers ("Cxcr6", )
  "Il2ra","Tnfrsf18","Ctla4","Lag3","Pdcd1",  # T Regular markers
  "Mki67","Stmn1", # proliferation
  "Ncr1","Tyrobp" # NK markers
)) + NoLegend()
  
# 9.4 cell type annotation

T_NK_lymph@active.ident <- T_NK_lymph$integrated_snn_res.0.8
head(T_NK_lymph@active.ident)
DimPlot(T_NK_lymph,label = T,group.by = "integrated_snn_res.0.8",label.size = 5,pt.size = 1)

T_NK_lymph$id <- T_NK_lymph@active.ident 

head(T_NK_lymph@meta.data)
new.cluster.ids <- c("CD4+ Tem",
                     "CD4+ Te",
                     "CD8+ CTL",
                     "CD4+ Tem",
                     "CD4+ Tem",
                     "CD8+ Tn",
                     "CD8+ CTL",
                     "CD4+ Tn",
                     "CD4+ Treg",
                     "CD4+ Te",
                     "CD8+ CTL",
                     "T pro",
                     "CD4+ Te",
                     "CD4+ Treg",
                     "NK",
                     "CD8+ Tem")

names(new.cluster.ids) <- levels(T_NK_lymph)
T_NK_lymph <- RenameIdents(T_NK_lymph, new.cluster.ids)

levels(T_NK_lymph) <- c("CD4+ Tn",
                        "CD4+ Te",
                        "CD4+ Tem",
                        "CD4+ Treg",
                        "CD8+ Tn",
                        "CD8+ CTL",
                        "CD8+ Tem",
                        "T pro",
                        "NK")

p1 <- DimPlot(T_NK_lymph,label = TRUE,reduction = "umap",label.size = 5,pt.size = 1) + scale_color_lancet() + theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title = "T_NK_lymph cells")
p2 <- DimPlot(T_NK_lymph,label = F,group.by = "type",label.size = 5,pt.size = 1,cols = c('Control' = '#00BFC4', 'Treatment' = '#F8766D')) + theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title =  paste0("Dataset type"))
p1 + p2

p3 <- DimPlot(T_NK_lymph,label = TRUE,reduction = "umap",label.size = 5,,pt.size = 1,split.by = "type") + scale_color_lancet()
p3

# DimPlot(T_NK_lymph,label = F,group.by = "sample",label.size = 5)+ theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title =  paste0("Sample id")) + scale_color_nejm()
# DimPlot(T_NK_lymph,label = T,group.by = "id",label.size = 5) + theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title =  paste0("Cluster id"))

# 9.5 draw heatmap plot

AverageExpression_value <-  AverageExpression(T_NK_lymph, assays = "SCT", features = c(
  "Cd3d","Il7r","Cd4","Cd8a", # T cell markers
  "Ccr7","Sell","Lef1","Tcf7", # T Naive markers
  "Il2","Il4","Il6","Il17a", # T Effector (,"Il9","Il10","Il22")
  "Ifng","Fasl","Nkg7","Gzma","Gzmk",# T Cytotoxic  markers
  "S100a4","Cxcr3", "Cd40lg","Cd69", # T Memory markers ("Cxcr6", )
  "Il2ra","Tnfrsf18","Ctla4","Lag3","Pdcd1",  # T Regular markers
  "Mki67","Stmn1", # proliferation
  "Ncr1","Tyrobp" # NK markers
), return.seurat = FALSE, group.by = "ident")

pheatmap::pheatmap(AverageExpression_value$SCT,cluster_rows = F,border_color = "black",
                   angle_col = 45,gaps_col = c(4,7,8),gaps_row = c(4,8,12,17,21,26,28,30),
                   cellwidth = 20,cellheight = 10,cluster_cols = F,scale = "row")


# 9.6 save metadata andT_NK_lymph.rds

T_NK_lymph$cell_type_sub <- T_NK_lymph@active.ident
write.csv(T_NK_lymph@meta.data,"T_NK_lymph_meta_data.csv")
saveRDS(T_NK_lymph,"3.Result/6.2 T&NK/T_NK_lymph_anno.rds")

# 9.7 count the proportions of each group and each sample in different subclusters

dfsam1 <- as.data.frame(table(T_NK_lymph$type,T_NK_lymph$cell_type_sub))
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

dfsam2 <- as.data.frame(table(T_NK_lymph$sample,T_NK_lymph$cell_type_sub))
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

# 9.8 draw pipe chart

dfsam1 <- as.data.frame(table(T_NK$type,T_NK$cell_type_sub))

control <- dfsam1[dfsam1$Var1 == "Control",]
ggplot(control, aes(x = "", y = Freq,fill = Var2)) + 
  geom_bar(stat = "identity")+
  coord_polar(theta = "y") + 
  theme_bw() + 
  labs(x = "", y = "", title = "")+
  theme(axis.ticks = element_blank())+
  theme(axis.text.x = element_blank()) + 
  theme(panel.grid=element_blank()) +   
  theme(panel.border=element_blank())+
  scale_fill_lancet(labels = myLabel) # 6x6

myLabel = as.vector(control$Var2)
myLabel = paste(myLabel, "(", round(control$Freq / sum(control$Freq) * 100, 2), "%)"
                , sep = "")  

treatment <- dfsam1[dfsam1$Var1 == "Treatment",]
ggplot(treatment, aes(x = "", y = Freq,fill = Var2)) + 
  geom_bar(stat = "identity")+
  coord_polar(theta = "y") + 
  theme_bw() + 
  labs(x = "", y = "", title = "")+
  theme(axis.ticks = element_blank())+
  theme(axis.text.x = element_blank()) + 
  theme(panel.grid=element_blank()) +    
  theme(panel.border=element_blank())+
  scale_fill_lancet(labels = myLabel) # 6x6

myLabel = as.vector(treatment$Var2)
myLabel = paste(myLabel, "(", round(treatment$Freq / sum(treatment$Freq) * 100, 2), "%)"
                , sep = "")  

# 9.9 T lymph & NK cell different expressed gene (DEGs) and GO enrichment analysis

T_NK <- readRDS("1.Data/5.2 SCT_T_NK_lymph_anno.rds")
T_NK@active.assay <- "SCT"
T_NK@active.ident <- as.factor(T_NK$type)
DimPlot(T_NK)
table(T_NK@active.ident)

head(T_NK@meta.data)

CD4_Tn  <- subset(T_NK,cell_type_sub == "CD4+ Tn")
CD4_Te  <- subset(T_NK,cell_type_sub == "CD4+ Te")
CD4_Tem <- subset(T_NK,cell_type_sub == "CD4+ Tem")
CD4_Treg   <- subset(T_NK,cell_type_sub == "CD4+ Treg")
CD8_Tn  <- subset(T_NK,cell_type_sub == "CD8+ Tn")
CD8_CTL  <- subset(T_NK,cell_type_sub == "CD8+ CTL")
CD8_Tem <- subset(T_NK,cell_type_sub == "CD8+ Tem")
T_pro <- subset(T_NK,cell_type_sub == "T pro")
NK <- subset(T_NK,cell_type_sub == "NK")

CD4_Tn_DEG <- FindMarkers(CD4_Tn,ident.1 = "Treatment",ident.2 = "Control",
                          min.pct = 0.1, logfc.threshold = 0.25) # n = 818
CD4_Tn_DEG$change = ifelse(CD4_Tn_DEG$p_val_adj < 0.05 & abs(CD4_Tn_DEG$avg_log2FC) >= 0.25, 
                         ifelse(CD4_Tn_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(CD4_Tn_DEG$change)
# Down Stable     Up 
# 271    458     89

CD4_Te_DEG <- FindMarkers(CD4_Te,ident.1 = "Treatment",ident.2 = "Control",
                        min.pct = 0.1, logfc.threshold = 0.25) # n = 2390
CD4_Te_DEG$change = ifelse(CD4_Te_DEG$p_val_adj < 0.05 & abs(CD4_Te_DEG$avg_log2FC) >= 0.25, 
                           ifelse(CD4_Te_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(CD4_Te_DEG$change)
# Down Stable     Up 
# 1408     77    905 

CD4_Tem_DEG <- FindMarkers(CD4_Tem,ident.1 = "Treatment",ident.2 = "Control",
                          min.pct = 0.1, logfc.threshold = 0.25) # n = 922
CD4_Tem_DEG$change = ifelse(CD4_Tem_DEG$p_val_adj < 0.05 & abs(CD4_Tem_DEG$avg_log2FC) >= 0.25, 
                           ifelse(CD4_Tem_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(CD4_Tem_DEG$change)
# Down Stable     Up 
# 509    206    207

CD4_Treg_DEG <- FindMarkers(CD4_Treg,ident.1 = "Treatment",ident.2 = "Control",
                            min.pct = 0.1, logfc.threshold = 0.25) # n = 1406
CD4_Treg_DEG$change = ifelse(CD4_Treg_DEG$p_val_adj < 0.05 & abs(CD4_Treg_DEG$avg_log2FC) >= 0.25, 
                           ifelse(CD4_Treg_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(CD4_Treg_DEG$change)
# Down Stable     Up 
# 434    952     20

CD8_Tn_DEG <- FindMarkers(CD8_Tn,ident.1 = "Treatment",ident.2 = "Control",
                            min.pct = 0.1, logfc.threshold = 0.25) # n = 817
CD8_Tn_DEG$change = ifelse(CD8_Tn_DEG$p_val_adj < 0.05 & abs(CD8_Tn_DEG$avg_log2FC) >= 0.25, 
                             ifelse(CD8_Tn_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(CD8_Tn_DEG$change)
# Down Stable     Up 
# 351    332    134 

CD8_CTL_DEG <- FindMarkers(CD8_CTL,ident.1 = "Treatment",ident.2 = "Control",
                            min.pct = 0.1, logfc.threshold = 0.25) # n = 838
CD8_CTL_DEG$change = ifelse(CD8_CTL_DEG$p_val_adj < 0.05 & abs(CD8_CTL_DEG$avg_log2FC) >= 0.25, 
                           ifelse(CD8_CTL_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(CD8_CTL_DEG$change)
# Down Stable     Up 
# 454    217    167

CD8_Tem_DEG <- FindMarkers(CD8_Tem,ident.1 = "Treatment",ident.2 = "Control",
                          min.pct = 0.1, logfc.threshold = 0.25) # n = 3113
CD8_Tem_DEG$change = ifelse(CD8_Tem_DEG$p_val_adj < 0.05 & abs(CD8_Tem_DEG$avg_log2FC) >= 0.25, 
                            ifelse(CD8_Tem_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(CD8_Tem_DEG$change)
# Down Stable     Up 
# 1253   1626    234

T_pro_DEG <- FindMarkers(T_pro,ident.1 = "Treatment",ident.2 = "Control",
                            min.pct = 0.1, logfc.threshold = 0.25) # n = 3031
T_pro_DEG$change = ifelse(T_pro_DEG$p_val_adj < 0.05 & abs(T_pro_DEG$avg_log2FC) >= 0.25, 
                            ifelse(T_pro_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(T_pro_DEG$change)
# Down Stable     Up 
# 1317   1467    247

NK_DEG <- FindMarkers(NK,ident.1 = "Treatment",ident.2 = "Control",
                            min.pct = 0.1, logfc.threshold = 0.25) # n = 1187
NK_DEG$change = ifelse(NK_DEG$p_val_adj < 0.05 & abs(NK_DEG$avg_log2FC) >= 0.25, 
                          ifelse(NK_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(NK_DEG$change)
# Down Stable     Up 
# 375    775     37

## extract the up-regularate gene of each subtypes 

CD4_Tn_DEG_up <- row.names(CD4_Tn_DEG[CD4_Tn_DEG$change == "Up",]) # n = 89
CD4_Te_DEG_up <- row.names(CD4_Te_DEG[CD4_Te_DEG$change == "Up",]) # n = 905
CD4_Tem_DEG_up <- row.names(CD4_Tem_DEG[CD4_Tem_DEG$change == "Up",])# n = 207
CD4_Treg_DEG_up <- row.names(CD4_Treg_DEG[CD4_Treg_DEG$change == "Up",])# n = 20
CD8_Tn_DEG_up <- row.names(CD8_Tn_DEG[CD8_Tn_DEG$change == "Up",])# n = 134
CD8_CTL_DEG_up <- row.names(CD8_CTL_DEG[CD8_CTL_DEG$change == "Up",])# n = 167
CD8_Tem_DEG_up <- row.names(CD8_Tem_DEG[CD8_Tem_DEG$change == "Up",])# n = 134
T_pro_DEG_up <- row.names(T_pro_DEG[T_pro_DEG$change == "Up",])# n = 247
NK_DEG_up <- row.names(NK_DEG[NK_DEG$change == "Up",])# n = 37

ALL_DEG_up <- union(CD4_Tn_DEG_up,c(CD4_Te_DEG_up,CD4_Tem_DEG_up,CD4_Treg_DEG_up,
                      CD8_Tn_DEG_up,CD8_CTL_DEG_up,CD8_Tem_DEG_up,
                      T_pro_DEG_up,NK_DEG_up))
ALL_DEG_up <- as.data.frame(ALL_DEG_up)
row.names(ALL_DEG_up) <- ALL_DEG_up$ALL_DEG_up
colnames(ALL_DEG_up) <- c("ID")
head(ALL_DEG_up)

T_NK_DEG_up <- as.data.frame(T_NK_DEG_up)
CD4_Tn_DEG_up <- as.data.frame(CD4_Tn_DEG_up)
CD4_Te_DEG_up <- as.data.frame(CD4_Te_DEG_up)
CD4_Tem_DEG_up <- as.data.frame(CD4_Tem_DEG_up)
CD4_Treg_DEG_up <- as.data.frame(CD4_Treg_DEG_up)
CD8_Tn_DEG_up <- as.data.frame(CD8_Tn_DEG_up)
CD8_CTL_DEG_up <- as.data.frame(CD8_CTL_DEG_up)
CD8_Tem_DEG_up <- as.data.frame(CD8_Tem_DEG_up)
T_pro_DEG_up <- as.data.frame(T_pro_DEG_up)
NK_DEG_up <- as.data.frame(NK_DEG_up)

rownames(T_NK_DEG_up) <- T_NK_DEG_up$T_NK_DEG_up
rownames(CD4_Tn_DEG_up) <- CD4_Tn_DEG_up$CD4_Tn_DEG_up
rownames(CD4_Te_DEG_up) <- CD4_Te_DEG_up$CD4_Te_DEG_up
rownames(CD4_Tem_DEG_up) <- CD4_Tem_DEG_up$CD4_Tem_DEG_up
rownames(CD4_Treg_DEG_up) <- CD4_Treg_DEG_up$CD4_Treg_DEG_up
rownames(CD8_Tn_DEG_up) <- CD8_Tn_DEG_up$CD8_Tn_DEG_up
rownames(CD8_CTL_DEG_up) <- CD8_CTL_DEG_up$CD8_CTL_DEG_up
rownames(CD8_Tem_DEG_up) <- CD8_Tem_DEG_up$CD8_Tem_DEG_up
rownames(T_pro_DEG_up) <- T_pro_DEG_up$T_pro_DEG_up
rownames(NK_DEG_up) <- NK_DEG_up$NK_DEG_up

CD4_Tn_DEG_up$ID <- rownames(CD4_Tn_DEG_up)
head(CD4_Tn_DEG_up)
CD4_Te_DEG_up$ID <- rownames(CD4_Te_DEG_up)
head(CD4_Te_DEG_up)
CD4_Tem_DEG_up$ID <- rownames(CD4_Tem_DEG_up)
head(CD4_Tem_DEG_up)
CD4_Treg_DEG_up$ID <- rownames(CD4_Treg_DEG_up)
head(CD4_Treg_DEG_up)
CD8_Tn_DEG_up$ID <- rownames(CD8_Tn_DEG_up)
head(CD8_Tn_DEG_up)
CD8_CTL_DEG_up$ID <- rownames(CD8_CTL_DEG_up)
head(CD8_CTL_DEG_up)
CD8_Tem_DEG_up$ID <- rownames(CD8_Tem_DEG_up)
head(CD8_Tem_DEG_up)
T_pro_DEG_up$ID <- rownames(T_pro_DEG_up)
head(T_pro_DEG_up)
NK_DEG_up$ID <- rownames(NK_DEG_up)
head(NK_DEG_up)

merge_DEG_up <- left_join(ALL_DEG_up,CD4_Tn_DEG_up,by="ID") 
merge_DEG_up <- left_join(ALL_DEG_up,T_NK_DEG_up,by="ID") 
merge_DEG_up <- left_join(merge_DEG_up,CD4_Te_DEG_up,by="ID")
merge_DEG_up <- left_join(merge_DEG_up,CD4_Tem_DEG_up,by="ID")
merge_DEG_up <- left_join(merge_DEG_up,CD4_Treg_DEG_up,by="ID")
merge_DEG_up <- left_join(merge_DEG_up,CD8_Tn_DEG_up,by="ID")
merge_DEG_up <- left_join(merge_DEG_up,CD8_CTL_DEG_up,by="ID")
merge_DEG_up <- left_join(merge_DEG_up,CD8_Tem_DEG_up,by="ID")
merge_DEG_up <- left_join(merge_DEG_up,T_pro_DEG_up,by="ID")
merge_DEG_up <- left_join(merge_DEG_up,NK_DEG_up,by="ID")

head(merge_DEG_up)
write.csv(merge_DEG_up,"T&NK_merge_DEG_up.csv")

merge_DEG_up_t <- merge_DEG_up
rownames(merge_DEG_up_t) <- merge_DEG_up_t$ID

head(merge_DEG_up_t)

merge_DEG_up_t[which(!is.na(merge_DEG_up_t),arr.ind = T)]<-1
merge_DEG_up_t[which(is.na(merge_DEG_up_t),arr.ind = T)]<-0
head(merge_DEG_up_t)
merge_DEG_up_t <- as.data.frame(lapply(merge_DEG_up_t,as.numeric))
merge_DEG_up_t <- merge_DEG_up_t[,-1]
str(merge_DEG_up_t)

# merge_DEG_up_t$ID <- row.names(merge_DEG_up_t)
# row.names(merge_DEG_up_t) <- row.names(merge_DEG_up)
# merge_DEG_up_t <- as.data.frame(merge_DEG_up_t)
# str(merge_DEG_up_t)

colnames(merge_DEG_up_t) <- c("CD4+ Tn","CD4+ Te","CD4+ Tem",
                              "CD4+ Treg","CD8+ Tn","CD8+ CTL",
                              "CD8+ Tem","T pro","NK")

upset(merge_DEG_up_t, nsets = 10,nintersects = 35,
      mb.ratio = c(0.6, 0.4),
      order.by = c("freq"),
      decreasing = c(TRUE,T), 
      sets.bar.color = c("#ED0000","#ADB6B6","#AD002A","#42B540",
                         "#FDAF91","#925E9F", "#00468B","#1B1919","#0099B4"),
      mainbar.y.label = "Intersection size of DEG (Up)", sets.x.label = "DEG number",
      text.scale = 1.4)

ALL_DEG_up_t <- union(CD4_Tn_DEG_up,c(CD4_Te_DEG_up,CD4_Tem_DEG_up,CD4_Treg_DEG_up,
                                    CD8_Tn_DEG_up,CD8_CTL_DEG_up,CD8_Tem_DEG_up,
                                    T_pro_DEG_up,NK_DEG_up))
T_NK_DEG <- read.csv("T_NK_DEG_all.csv",row.names = 1,check.names = F)

length(intersect(ALL_DEG_up_t,T_NK_DEG_up)) # n =734

venn.plot <- 
  venn.diagram(
    x = list(
      'Union_subgroup_DEG' = ALL_DEG_up_t,
      'Whole_T&NK_DEG' = T_NK_DEG_up),
    filename = NULL,
    col = "black",
    fill = c("red", "skyblue"),
    alpha = 0.6,
    cex = 0.8,
    cat.col = 'black',
    cat.cex = 0.8,
    cat.fontface = "bold",
    margin = 0.05,
    # main = "Overlap protein detected in different methods",
    main.cex = 1.2)
pdf(file="Overlap DEG.pdf")
grid.draw(venn.plot)
dev.off()

## extract the down-regularate gene of each subtypes 

CD4_Tn_DEG_down <- row.names(CD4_Tn_DEG[CD4_Tn_DEG$change == "Down",]) # n = 89
CD4_Te_DEG_down <- row.names(CD4_Te_DEG[CD4_Te_DEG$change == "Down",]) # n = 905
CD4_Tem_DEG_down <- row.names(CD4_Tem_DEG[CD4_Tem_DEG$change == "Down",])# n = 207
CD4_Treg_DEG_down <- row.names(CD4_Treg_DEG[CD4_Treg_DEG$change == "Down",])# n = 20
CD8_Tn_DEG_down <- row.names(CD8_Tn_DEG[CD8_Tn_DEG$change == "Down",])# n = 134
CD8_CTL_DEG_down <- row.names(CD8_CTL_DEG[CD8_CTL_DEG$change == "Down",])# n = 167
CD8_Tem_DEG_down <- row.names(CD8_Tem_DEG[CD8_Tem_DEG$change == "Down",])# n = 134
T_pro_DEG_down <- row.names(T_pro_DEG[T_pro_DEG$change == "Down",])# n = 247
NK_DEG_down <- row.names(NK_DEG[NK_DEG$change == "Down",])# n = 37

ALL_DEG_down <- union(CD4_Tn_DEG_down,c(CD4_Te_DEG_down,CD4_Tem_DEG_down,CD4_Treg_DEG_down,
                                        CD8_Tn_DEG_down,CD8_CTL_DEG_down,CD8_Tem_DEG_down,
                                        T_pro_DEG_down,NK_DEG_down))
ALL_DEG_down <- as.data.frame(ALL_DEG_down)
row.names(ALL_DEG_down) <- ALL_DEG_down$ALL_DEG_down
colnames(ALL_DEG_down) <- c("ID")
head(ALL_DEG_down)

CD4_Tn_DEG_down <- as.data.frame(CD4_Tn_DEG_down)
CD4_Te_DEG_down <- as.data.frame(CD4_Te_DEG_down)
CD4_Tem_DEG_down <- as.data.frame(CD4_Tem_DEG_down)
CD4_Treg_DEG_down <- as.data.frame(CD4_Treg_DEG_down)
CD8_Tn_DEG_down <- as.data.frame(CD8_Tn_DEG_down)
CD8_CTL_DEG_down <- as.data.frame(CD8_CTL_DEG_down)
CD8_Tem_DEG_down <- as.data.frame(CD8_Tem_DEG_down)
T_pro_DEG_down <- as.data.frame(T_pro_DEG_down)
NK_DEG_down <- as.data.frame(NK_DEG_down)

rownames(CD4_Tn_DEG_down) <- CD4_Tn_DEG_down$CD4_Tn_DEG_down
rownames(CD4_Te_DEG_down) <- CD4_Te_DEG_down$CD4_Te_DEG_down
rownames(CD4_Tem_DEG_down) <- CD4_Tem_DEG_down$CD4_Tem_DEG_down
rownames(CD4_Treg_DEG_down) <- CD4_Treg_DEG_down$CD4_Treg_DEG_down
rownames(CD8_Tn_DEG_down) <- CD8_Tn_DEG_down$CD8_Tn_DEG_down
rownames(CD8_CTL_DEG_down) <- CD8_CTL_DEG_down$CD8_CTL_DEG_down
rownames(CD8_Tem_DEG_down) <- CD8_Tem_DEG_down$CD8_Tem_DEG_down
rownames(T_pro_DEG_down) <- T_pro_DEG_down$T_pro_DEG_down
rownames(NK_DEG_down) <- NK_DEG_down$NK_DEG_down

CD4_Tn_DEG_down$ID <- rownames(CD4_Tn_DEG_down)
head(CD4_Tn_DEG_down)
CD4_Te_DEG_down$ID <- rownames(CD4_Te_DEG_down)
head(CD4_Te_DEG_down)
CD4_Tem_DEG_down$ID <- rownames(CD4_Tem_DEG_down)
head(CD4_Tem_DEG_down)
CD4_Treg_DEG_down$ID <- rownames(CD4_Treg_DEG_down)
head(CD4_Treg_DEG_down)
CD8_Tn_DEG_down$ID <- rownames(CD8_Tn_DEG_down)
head(CD8_Tn_DEG_down)
CD8_CTL_DEG_down$ID <- rownames(CD8_CTL_DEG_down)
head(CD8_CTL_DEG_down)
CD8_Tem_DEG_down$ID <- rownames(CD8_Tem_DEG_down)
head(CD8_Tem_DEG_down)
T_pro_DEG_down$ID <- rownames(T_pro_DEG_down)
head(T_pro_DEG_down)
NK_DEG_down$ID <- rownames(NK_DEG_down)
head(NK_DEG_down)

merge_DEG_down <- left_join(ALL_DEG_down,CD4_Tn_DEG_down,by="ID") 
merge_DEG_down <- left_join(ALL_DEG_down,T_NK_DEG_down,by="ID") 
merge_DEG_down <- left_join(merge_DEG_down,CD4_Te_DEG_down,by="ID")
merge_DEG_down <- left_join(merge_DEG_down,CD4_Tem_DEG_down,by="ID")
merge_DEG_down <- left_join(merge_DEG_down,CD4_Treg_DEG_down,by="ID")
merge_DEG_down <- left_join(merge_DEG_down,CD8_Tn_DEG_down,by="ID")
merge_DEG_down <- left_join(merge_DEG_down,CD8_CTL_DEG_down,by="ID")
merge_DEG_down <- left_join(merge_DEG_down,CD8_Tem_DEG_down,by="ID")
merge_DEG_down <- left_join(merge_DEG_down,T_pro_DEG_down,by="ID")
merge_DEG_down <- left_join(merge_DEG_down,NK_DEG_down,by="ID")

head(merge_DEG_down)
write.csv(merge_DEG_down,"T&NK_merge_DEG_down.csv")

merge_DEG_down_t <- merge_DEG_down
rownames(merge_DEG_down_t) <- merge_DEG_down_t$ID

head(merge_DEG_down_t)

merge_DEG_down_t[which(!is.na(merge_DEG_down_t),arr.ind = T)]<-1
merge_DEG_down_t[which(is.na(merge_DEG_down_t),arr.ind = T)]<-0
head(merge_DEG_down_t)
merge_DEG_down_t <- as.data.frame(lapply(merge_DEG_down_t,as.numeric))
merge_DEG_down_t <- merge_DEG_down_t[,-1]
str(merge_DEG_down_t)

# merge_DEG_down_t$ID <- row.names(merge_DEG_down_t)
# row.names(merge_DEG_down_t) <- row.names(merge_DEG_down)
# merge_DEG_down_t <- as.data.frame(merge_DEG_down_t)
# str(merge_DEG_down_t)

colnames(merge_DEG_down_t) <- c("CD4+ Tn","CD4+ Te","CD4+ Tem",
                                "CD4+ Treg","CD8+ Tn","CD8+ CTL",
                                "CD8+ Tem","T pro","NK")

upset(merge_DEG_down_t, nsets = 10,nintersects = 35,
      mb.ratio = c(0.6, 0.4),
      order.by = c("freq"),
      decreasing = c(TRUE,T), 
      sets.bar.color = c("#ED0000",	"#ADB6B6","#AD002A","#42B540",
                         "#FDAF91","#0099B4","#1B1919","#925E9F","#00468B"),
      mainbar.y.label = "Intersection size of DEG (down)", sets.x.label = "DEG number",
      text.scale = 1.4)

ALL_DEG_down_t <- union(CD4_Tn_DEG_down,c(CD4_Te_DEG_down,CD4_Tem_DEG_down,CD4_Treg_DEG_down,
                                          CD8_Tn_DEG_down,CD8_CTL_DEG_down,CD8_Tem_DEG_down,
                                          T_pro_DEG_down,NK_DEG_down))
T_NK_DEG <- read.csv("T_NK_DEG_all.csv",row.names = 1,check.names = F)
T_NK_DEG_down <- row.names(T_NK_DEG[T_NK_DEG$change == "Down",])

length(intersect(ALL_DEG_down_t,T_NK_DEG_down)) # n = 986

venn.plot <- 
  venn.diagram(
    x = list(
      'Union_subgroup_DEG' = ALL_DEG_down_t,
      'Whole_T&NK_DEG' = T_NK_DEG_down),
    filename = NULL,
    col = "black",
    fill = c("red", "skyblue"),
    alpha = 0.6,
    cex = 0.8,
    cat.col = 'black',
    cat.cex = 0.8,
    cat.fontface = "bold",
    margin = 0.05,
    # main = "Overlap protein detected in different methods",
    main.cex = 1.2)
pdf(file="Overlap DEG(down).pdf")
grid.draw(venn.plot)
dev.off()

# 9.10 T lymph & NK cell different expressed gene (DEGs) and GO enrichment analysis

T_NK <- readRDS("1.Data/5.2 SCT_T_NK_lymph_anno.rds")
T_NK@active.assay <- "SCT"
T_NK@active.ident <- as.factor(T_NK$type)
DimPlot(T_NK)

T_NK_DEG <- FindMarkers(T_NK,ident.1 = "Treatment",ident.2 = "Control",
                        min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
head(T_NK_DEG)

T_NK_DEG$change = ifelse(T_NK_DEG$p_val_adj < 0.05 & abs(T_NK_DEG$avg_log2FC) >= 0.25, 
                         ifelse(T_NK_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(T_NK_DEG$change)
# Down Stable     Up 
# 990      2    784 

write.csv(T_NK_DEG,"T_NK_DEG_all.csv")

T_NK_DEG <- read.csv("T_NK_DEG_all.csv",row.names = 1,check.names = F)
T_NK_DEG_up <- row.names(T_NK_DEG[T_NK_DEG$change == "Up",]) # n = 784
head(T_NK_DEG_up)

T_NK_DEG_down <- row.names(T_NK_DEG[T_NK_DEG$change == "Down",]) # n = 990
head(T_NK_DEG_down)

T_NK@active.ident <- as.factor(T_NK$sample)
DimPlot(T_NK)
T_NK_averge <- AverageExpression(T_NK)
head(T_NK_averge$SCT)
T_NK_exper <- T_NK_averge$SCT[c(T_NK_DEG_up,T_NK_DEG_down),]

pdf("pheatmap in T_NK.pdf",width = 5,height = 10)

annotation_col <- data.frame(Sample = factor(c("KA1","KA2","KA3","KC1","KC2","KC3")),Type = factor(c(rep("Treatment",3),rep("Control",3))))
row.names(annotation_col) <- c("KA1","KA2","KA3","KC1","KC2","KC3")
ann_colors = list(Sample = c(KA1="#BC3C29",KA2="#0072B5",KA3="#E18727",KC1="#20854E",KC2="#7876B1",KC3="#6F99AD" ), CellType = c(Control = "#00BFC4",Treatment = "#F8766D"))

pheatmap::pheatmap(T_NK_exper,show_colnames = T,show_rownames = F,
                   cellwidth = 15,cellheight = 0.2,scale = "row",
                   treeheight_row = 0,treeheight_col = 0,
                   cluster_rows = F,cluster_cols = F,border_color = "black",
                   annotation_col = annotation_col,annotation_colors = ann_colors)
dev.off()

Go_BP_up <- enrichGO(gene = T_NK_DEG_up, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                     ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)

Go_BP_down <- enrichGO(gene = T_NK_DEG_down, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                       ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)

p1 <- dotplot(Go_BP_up, showCategory = 10,title = "The GO enrichment(BP) analysis of AA vs Con up DEG in T_NK")
p2 <- dotplot(Go_BP_down, showCategory = 10,title = "The GO enrichment(BP) analysis of AA vs Con down DEG in T_NK")

pdf("GO enrichment analysis of Treatment vs Control DEG in T_NK.pdf",width = 10,height = 10)

p1/p2

write.csv(Go_BP_up@result,"The GO enrichment(BP) analysis of AA vs Con up DEG in T_NK.csv" )
write.csv(Go_BP_down@result,"The GO enrichment(BP) analysis of AA vs Con down DEG in T_NK.csv" )

# 9.11 T lymph cell trajectory analysis

T_NK <- readRDS(" SCT_T_NK_lymph_anno.rds")
DefaultAssay(T_NK) <- "RNA"
DimPlot(T_NK)

expr_matrix <- as(as.matrix(T_NK@assays$RNA@counts), 'sparseMatrix')

p_data <- T_NK@meta.data

f_data <- data.frame(gene_short_name = row.names(T_NK),row.names = row.names(T_NK))

pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)

cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size()) 

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

cds <- detectGenes(cds, min_expr = 0.1)
print(head(fData(cds)))

expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >= 10))
length(expressed_genes) # m = 15878

diff <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr="~cell_type_sub",cores=10)

head(diff)

deg <- subset(diff, qval < 0.001) 
deg <- deg[order(deg$qval,decreasing=F),]
head(deg)

write.csv(deg,file="monocle_DEG.csv",quote=F)

ordergene <- rownames(deg) 
cds <- setOrderingFilter(cds, ordergene)  
table(cds@featureData@data[["use_for_ordering"]]) 
plot_ordering_genes(cds) 

cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')

cds <- orderCells(cds)
cds2 <- orderCells(cds,root_state = 3) 

# p0 <- plot_cell_trajectory(cds,color_by="State",show_tree = T,cell_size = 1,cell_link_size = 2,show_backbone=TRUE) 
p1 <- plot_cell_trajectory(cds,color_by="Pseudotime",show_tree = T,cell_size = 1,cell_link_size = 1.5,show_backbone=TRUE) 
p2 <- plot_cell_trajectory(cds,color_by="type",show_tree = T,cell_size = 1,cell_link_size = 1.5,show_backbone=TRUE) + scale_color_manual(breaks = c("Control", "Treatment"), values=c("#00BFC4", "#F8766D"))
p3 <- plot_cell_trajectory(cds,color_by="cell_type_sub",show_tree = T,cell_size = 1,cell_link_size = 1.5,show_backbone=TRUE) + scale_color_lancet()
p4 <- plot_cell_trajectory(cds,color_by="sample",show_tree = T,cell_size = 1,cell_link_size = 2,show_backbone=TRUE)+ scale_color_nejm()

p1 + p2 + p3 +p4

plot_cell_trajectory(cds, color_by = "cell_type_sub", cell_size = 1,cell_link_size = 2,show_backbone=TRUE) + facet_wrap("~type", nrow = 1)+ scale_color_lancet()


table(T_NK$cell_type_sub)


Time_diff <- differentialGeneTest(cds[ordergene,], cores = 5, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
head(Time_diff)
Time_diff <- Time_diff[,c(5,2,3,4,1,6,7)]

Time_genes <- top_n(Time_diff, n = 1, desc(qval)) %>% pull(gene_short_name) %>% as.character()
head(Time_genes)

p <- plot_pseudotime_heatmap(cds[Time_genes,], num_clusters = 4, show_rownames=T, return_heatmap=T)
p$tree_row

clusters <- cutree(p$tree_row, k = 4)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
table(clustering)
# clustering
# 1   2   3   4 
# 160  22 161  11

clustering$gene <- row.names(clustering)
C1 <- row.names(clustering[clustering$Gene_Clusters == "1",])
C1_Go_BP_up <- enrichGO(gene = C1, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                     ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
dotplot(C1_Go_BP_up, showCategory = 10,title = "The GO enrichment(BP) analysis of C1")


C2 <- row.names(clustering[clustering$Gene_Clusters == "2",])
C2_Go_BP_up <- enrichGO(gene = C2, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                        ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
dotplot(C2_Go_BP_up, showCategory = 10,title = "The GO enrichment(BP) analysis of C2")

C3 <- row.names(clustering[clustering$Gene_Clusters == "3",])
C3_Go_BP_up <- enrichGO(gene = C3, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                        ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
dotplot(C3_Go_BP_up, showCategory = 10,title = "The GO enrichment(BP) analysis of C3")


C4 <- row.names(clustering[clustering$Gene_Clusters == "4",])
C4_Go_BP_up <- enrichGO(gene = C4, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                        ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
dotplot(C4_Go_BP_up, showCategory = 10,title = "The GO enrichment(BP) analysis of C3")

T_NK_DEG <- read.csv("../4.DEG_analysis/T_NK_DEG.csv",row.names = 1)
T_NK_DEG_up <- row.names(T_NK_DEG[T_NK_DEG$change == "Up",])
T_NK_DEG_down <- row.names(T_NK_DEG[T_NK_DEG$change == "Down",])

length(intersect(T_NK_DEG_up,c(C1,C2))) # n=150
length(intersect(T_NK_DEG_down,c(C3,C4))) # n=161

genes <- row.names(subset(fData(cds),gene_short_name %in% c())) 
genes <- c("Cd4","Cd8a","Cd8b1","Ccr7","Lef1","Tcf7","S100a4","Cxcr3")

genes <- c("Il2","Il4","Il6","Il17a",
"Ifng", "Fasl","Nkg7","Gzma",
"Il2ra","Tnfrsf18","Ctla4","Pdcd1")

plot_genes_in_pseudotime(cds[genes,], color_by = "type",ncol = 4,nrow = 3,label_by_short_name = F,relative_expr = T)+ 
  scale_color_manual(breaks = c("Control", "Treatment"), values=c("#00BFC4", "#F8766D")) # 8x5