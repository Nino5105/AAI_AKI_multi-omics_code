# 7.Other Epithelial subsets cell

# 7.1 library packages

library(dplyr)
library(Seurat)
library(patchwork)
library(clustree)
library(ggsci)
library(ggplot2) 
library(clusterProfiler)
library(org.Mm.eg.db)

Other_segment <- readRDS("Other_segment.rds")
DimPlot(Other_segment,)

# 7.2 datasets dimension reduction and re-cluster

Other_segment <- RunPCA(object = Other_segment,verbose = FALSE)
ElbowPlot(Other_segment,ndims = 20) # 选择主成分为15

Other_segment <- FindNeighbors(Other_segment,dim=1:25)
Other_segment <- FindClusters(Other_segment,resolution = 0.8)
Other_segment <- RunUMAP (Other_segment,reduction="pca", dims = 1:20)
Other_segment@active.ident <- Other_segment$integrated_snn_res.0.8
DimPlot(Other_segment,label = T,pt.size = 1,label.size = 5)

# Other_segment@assays

# 7.3 marker gene expression

 DoHeatmap(Other_segment, slot = "scale.data",features = c(
  "Aqp1","Bst1", # DLH markers
  "Slc12a1","Umod","Cldn16", # ALH markers
  "Slc12a3", # DCT
  "Atp6v0d2","Atp6v1g3", # CD-IC
  "Aqp2","Hsd11b2",# CD-PC
  "Wt1","Lamb2","Podxl","Ptpro"# Podo
 )) + NoLegend()

table(Idents(Other_segment),Other_segment$orig.ident)

# 7.4 cell type annotation

Other_segment@active.ident <- Other_segment$integrated_snn_res.0.8

new.cluster.ids=c("DCT","ALH","ALH","DLH","Podo","DLH","DCT","ALH","CD-PC","ALH","ALH","Podo","CD-PC",
                  "ALH","Podo","Podo","CD-IC","DLH","CD-IC","ALH","DLH","Podo","Podo")
names(new.cluster.ids) <- levels(Other_segment)
Other_segment <- RenameIdents(Other_segment, new.cluster.ids)
levels(Other_segment) <- c("DLH","ALH","DCT","CD-IC","CD-PC","Podo")
Other_segment$cell_type <- Other_segment@active.ident
DimPlot(Other_segment)

head(Other_segment@active.ident)
DimPlot(Other_segment,label = T,group.by = "cell_type",label.size = 5,pt.size = 1)

p1 <- DimPlot(Other_segment,label = TRUE,reduction = "umap",label.size = 5,pt.size = 1,cols = c("DLH" = '#53A85F', "ALH" = '#F1BB72',"DCT" = '#F3B1A0',
                                                                                               "CD-IC" = '#D6E7A3',"CD-PC" = '#57C3F3',"Podo" = "#E59CC4")) + theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title = "Other_segment cells")
p2 <- DimPlot(Other_segment,label = F,group.by = "type",label.size = 5,pt.size = 1,cols = c('Control' = '#00BFC4', 'Treatment' = '#F8766D')) + theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title =  paste0("Dataset type"))
p1 + p2

p3 <- DimPlot(Other_segment,label = TRUE,reduction = "umap",label.size = 5,,pt.size = 1,split.by = "type",cols = c("DLH" = '#53A85F', "ALH" = '#F1BB72',"DCT" = '#F3B1A0',
                                                                                                                   "CD-IC" = '#D6E7A3',"CD-PC" = '#57C3F3',"Podo" = "#E59CC4"))
p3

# DimPlot(Other_segment,label = F,group.by = "sample",label.size = 5)+ theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title =  paste0("Sample id")) + scale_color_nejm()
# DimPlot(Other_segment,label = T,group.by = "id",label.size = 5) + theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title =  paste0("Cluster id"))

# 7.5 draw heatmap plot

AverageExpression_value <-  AverageExpression(Other_segment, assays = "SCT", features = c(
  "Aqp1","Bst1", # DLH markers
  "Slc12a1","Umod","Cldn16", # ALH markers
  "Slc12a3", # DCT
  "Atp6v0d2","Atp6v1g3", # CD-IC
  "Aqp2","Hsd11b2",# CD-PC
  "Wt1","Lamb2","Podxl","Ptpro"# Podo
), return.seurat = FALSE, group.by = "ident")

pheatmap::pheatmap(AverageExpression_value$SCT,cluster_rows = F,border_color = "black",
                   angle_col = 45,#gaps_col = c(4,7,8),gaps_row = c(4,8,12,17,21,26,28,30),
                   cellwidth = 25,cellheight = 20,cluster_cols = F,scale = "row")


# 7.6 save metadata and Other epithelial.rds

Other_segment$cell_type <- Other_segment@active.ident
write.csv(Other_segment@meta.data,"3.Result/6.4 Other segment/Other_segment_meta_data.csv")
saveRDS(Other_segment,"1.Data/5.5 Other_segment_anno.rds")

# 7.7 count the proportions of each group and each sample in different subclusters

dfsam1 <- as.data.frame(table(Other_segment$type,Other_segment$cell_type))[c(3:12,17:18),]
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

dfsam2 <- as.data.frame(table(Other_segment$sample,Other_segment$cell_type))[c(7:36,49:54),]
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


p4/p5 

# 6.8 Other epithelial different expressed gene (DEGs) and GO enrichment analysis

Other_epi <- readRDS("1.Data/5.5 Other_epi_anno.rds")

Other_epi_count <- Other_epi@assays$RNA@counts
Other_epi_count <- as.matrix(Other_epi_count)
head(Other_epi_count)[1:5,1:5]
write.csv(Other_epi_count,"Other_epi_count.csv")


DimPlot(Other_epi,group.by = "cell_type")
table(Other_epi$cell_type)
# DLH   ALH   DCT CD-IC CD-PC  Podo 
# 1201  2416  1239   272   481  1166
 
# BPA.total.normalize.DF.Singlet.Acrosomal <- subset(BPA.total.normalize.DF.Singlet,idents = "Acrosomal")
# Idents(BPA.total.normalize.DF.Singlet.Acrosomal) <- 'group'

DimPlot(Other_epi)

DLH <- subset(Other_epi,idents = "DLH")
ALH <- subset(Other_epi,idents = "ALH")
DCT <- subset(Other_epi,idents = "DCT")
CD_IC <- subset(Other_epi,idents = "CD-IC")
CD_PC <- subset(Other_epi,idents = "CD-PC")
Podo <- subset(Other_epi,idents = "Podo")

Idents(Other_Epi_DEG) <- as.factor(Other_Epi_DEG$type)
Idents(DLH) <- as.factor(DLH$type)
Idents(ALH) <- as.factor(ALH$type)
Idents(DCT) <- as.factor(DCT$type)
Idents(CD_IC) <- as.factor(CD_IC$type)
Idents(CD_PC) <- as.factor(CD_PC$type)
Idents(Podo) <- as.factor(Podo$type)

Other_Epi_DEG 
Other_Epi_DEG <- FindMarkers(Other_Epi_DEG,ident.1 = 'Treatment', ident.2 = 'Control', verbose = FALSE,min.pct = 0.1, logfc.threshold = 0.25)
Other_Epi_DEG$Type <-  ifelse(Other_Epi_DEG$p_val < 0.05 & abs(Other_Epi_DEG$avg_log2FC) >= 0.25, 
                    ifelse(Other_Epi_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(Other_Epi_DEG$Type)
# Down Stable     Up 
# 794    372    668 
# Other_Epi_DEG$Cell_type <- "All"
# table(Other_Epi_DEG$Type)

DLH_DEG <- FindMarkers(DLH,ident.1 = 'Treatment', ident.2 = 'Control', verbose = FALSE,min.pct = 0.1, logfc.threshold = 0.25)
DLH_DEG$Type <-  ifelse(DLH_DEG$p_val < 0.05 & abs(DLH_DEG$avg_log2FC) >= 0.25, 
                        ifelse(DLH_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(DLH_DEG$Type)
# Down Stable     Up 
# 324    392    888
DLH_DEG$Cell_type <- "DLH"

ALH_DEG <- FindMarkers(ALH,ident.1 = 'Treatment', ident.2 = 'Control', verbose = FALSE,min.pct = 0.1, logfc.threshold = 0.25)
ALH_DEG$Type <-  ifelse(ALH_DEG$p_val < 0.05 & abs(ALH_DEG$avg_log2FC) >= 0.25, 
                        ifelse(ALH_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(ALH_DEG$Type)
# Down Stable     Up 
# 396    296    710
ALH_DEG$Cell_type <- "ALH"

DCT_DEG <- FindMarkers(DCT,ident.1 = 'Treatment', ident.2 = 'Control', verbose = FALSE,min.pct = 0.1, logfc.threshold = 0.25)
DCT_DEG$Type <-  ifelse(DCT_DEG$p_val < 0.05 & abs(DCT_DEG$avg_log2FC) >= 0.25, 
                        ifelse(DCT_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(DCT_DEG$Type)
# Down Stable     Up 
# 319    375    568
DCT_DEG$Cell_type <- "DCT"

CD_IC_DEG <- FindMarkers(CD_IC,ident.1 = 'Treatment', ident.2 = 'Control', verbose = FALSE,min.pct = 0.1, logfc.threshold = 0.25)
CD_IC_DEG$Type <-  ifelse(CD_IC_DEG$p_val < 0.05 & abs(CD_IC_DEG$avg_log2FC) >= 0.25, 
                          ifelse(CD_IC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(CD_IC_DEG$Type)
# Down Stable     Up 
# 228    582    376 
CD_IC_DEG$Cell_type <- "CD-IC"

CD_PC_DEG <- FindMarkers(CD_PC,ident.1 = 'Treatment', ident.2 = 'Control', verbose = FALSE,min.pct = 0.1, logfc.threshold = 0.25)
CD_PC_DEG$Type <-  ifelse(CD_PC_DEG$p_val < 0.05 & abs(CD_PC_DEG$avg_log2FC) >= 0.25, 
                          ifelse(CD_PC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(CD_PC_DEG$Type)
# Down Stable     Up 
# 381    543    328
CD_PC_DEG$Cell_type <- "CD-PC"

Podo_DEG <- FindMarkers(Podo,ident.1 = 'Treatment', ident.2 = 'Control', verbose = FALSE,min.pct = 0.1, logfc.threshold = 0.25)
Podo_DEG$Type <-  ifelse(Podo_DEG$p_val < 0.05 & abs(Podo_DEG$avg_log2FC) >= 0.25, 
                         ifelse(Podo_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(Podo_DEG$Type)
# Down Stable     Up 
# 793    645    549 
Podo_DEG$Cell_type <- "Podo"

# differential genes visualization
ALL_DEG <- rbind(DLH_DEG,ALH_DEG,DCT_DEG,CD_IC_DEG,CD_PC_DEG,Podo_DEG)
ALL_DEG$Cell_type  <- factor(ALL_DEG$Cell_type,levels = rev(c("DLH","ALH","DCT","CD-IC","CD-PC","Podo")))

ALL_DEG_sub <- ALL_DEG[ALL_DEG$Type != "Stable",]
table(ALL_DEG_sub$Type)

ggplot() + geom_point() +
  geom_jitter(ALL_DEG_sub, mapping=aes(x=Cell_type, y= avg_log2FC, color=Type),width=0.4,size=1) +
  theme_test()    +
  theme(axis.text.x=element_text(size=10,angle=90,face ='bold'),
        axis.text.y=element_text(size=10,face ='bold'),
        axis.title.x=element_text(size = 10),
        axis.title.y=element_text(size = 14)) + 
  xlab(NULL) + ylab('Log2FC(Treatment vs Control)')  + 
  geom_hline(yintercept = 0,lty=1,lwd=2,alpha=0.5)+
  scale_colour_manual(values = c("#6B9AC7","#E24D36")) +
  guides(color=guide_legend(override.aes = list(size=6))) + coord_flip() # 6X10

# extract the up-regulated genes for enrichment analysis
DLH_DEG_up <- rownames(DLH_DEG[DLH_DEG$Type =="Up",])
ALH_DEG_up <- rownames(ALH_DEG[ALH_DEG$Type =="Up",])
DCT_DEG_up <- rownames(DCT_DEG[DCT_DEG$Type =="Up",])
CD_IC_DEG_up <- rownames(CD_IC_DEG[CD_IC_DEG$Type =="Up",])
CD_PC_DEG_up <- rownames(CD_PC_DEG[CD_PC_DEG$Type =="Up",])
Podo_DEG_up <- rownames(Podo_DEG[Podo_DEG$Type =="Up",])

DLH_Go_BP_up <- enrichGO(gene = DLH_DEG_up, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                         ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
ALH_Go_BP_up <- enrichGO(gene = ALH_DEG_up, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                         ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
DCT_Go_BP_up <- enrichGO(gene = DCT_DEG_up, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                         ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
CD_IC_Go_BP_up <- enrichGO(gene = CD_IC_DEG_up, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                           ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
CD_PC_Go_BP_up <- enrichGO(gene = CD_PC_DEG_up, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                           ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
Podo_Go_BP_up <- enrichGO(gene = Podo_DEG_up, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                          ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)

DLH_Go_BP_up_5 <- DLH_Go_BP_up@result[c(1:5),c("Description","pvalue")]
DLH_Go_BP_up_5$log10Pvalue <- -log10(DLH_Go_BP_up_5$pvalue)
DLH_Go_BP_up_5$subtype <- "DLH"

ALH_Go_BP_up_5 <- ALH_Go_BP_up@result[c(1:5),c("Description","pvalue")]
ALH_Go_BP_up_5$log10Pvalue <- -log10(ALH_Go_BP_up_5$pvalue)
ALH_Go_BP_up_5$subtype <- "ALH"

DCT_Go_BP_up_5 <- DCT_Go_BP_up@result[c(1:5),c("Description","pvalue")]
DCT_Go_BP_up_5$log10Pvalue <- -log10(DCT_Go_BP_up_5$pvalue)
DCT_Go_BP_up_5$subtype <- "DCT"

CD_IC_Go_BP_up_5 <- CD_IC_Go_BP_up@result[c(1:5),c("Description","pvalue")]
CD_IC_Go_BP_up_5$log10Pvalue <- -log10(CD_IC_Go_BP_up_5$pvalue)
CD_IC_Go_BP_up_5$subtype <- "CD-IC"

CD_PC_Go_BP_up_5 <- CD_PC_Go_BP_up@result[c(1:5),c("Description","pvalue")]
CD_PC_Go_BP_up_5$log10Pvalue <- -log10(CD_PC_Go_BP_up_5$pvalue)
CD_PC_Go_BP_up_5$subtype <- "CD-PC"

Podo_Go_BP_up_5 <- Podo_Go_BP_up@result[c(1:5),c("Description","pvalue")]
Podo_Go_BP_up_5$log10Pvalue <- -log10(Podo_Go_BP_up_5$pvalue)
Podo_Go_BP_up_5$subtype <- "Podo"

Go_Epi_DEG_up <- rbind(DLH_Go_BP_up_5,ALH_Go_BP_up_5,DCT_Go_BP_up_5,CD_IC_Go_BP_up_5,CD_PC_Go_BP_up_5,Podo_Go_BP_up_5)
Go_Epi_DEG_up$item <- row.names(Go_Epi_DEG_up)

Go_Epi_DEG_up$subtype <- factor(Go_Epi_DEG_up$subtype,levels = rev(c("DLH","ALH","DCT","CD-IC","CD-PC","Podo")))

ggbarplot(Go_Epi_DEG_up, x="item", y="log10Pvalue", fill = "subtype", color = "white",
          palette =  c("DLH" = '#53A85F', "ALH" = '#F1BB72',"DCT" = '#F3B1A0',
                       "CD-IC" = '#D6E7A3',"CD-PC" = '#57C3F3',"Podo" = "#E59CC4"),
          sort.val = "asc", 
          sort.by.groups=TRUE,
          x.text.angle=0, 
          xlab = NULL) + coord_flip() 

table(Go_Epi_DEG_up$Description)

write.csv(Go_Epi_DEG_up,"Go_Epi_DEG_up.csv")

## extract the down-regulated genes for enrichment analysis
DLH_DEG_down <- rownames(DLH_DEG[DLH_DEG$Type =="Down",])
ALH_DEG_down <- rownames(ALH_DEG[ALH_DEG$Type =="Down",])
DCT_DEG_down <- rownames(DCT_DEG[DCT_DEG$Type =="Down",])
CD_IC_DEG_down <- rownames(CD_IC_DEG[CD_IC_DEG$Type =="Down",])
CD_PC_DEG_down <- rownames(CD_PC_DEG[CD_PC_DEG$Type =="Down",])
Podo_DEG_down <- rownames(Podo_DEG[Podo_DEG$Type =="Down",])

DLH_Go_BP_down <- enrichGO(gene = DLH_DEG_down, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                         ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
ALH_Go_BP_down <- enrichGO(gene = ALH_DEG_down, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                         ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
DCT_Go_BP_down <- enrichGO(gene = DCT_DEG_down, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                         ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
CD_IC_Go_BP_down <- enrichGO(gene = CD_IC_DEG_down, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                           ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
CD_PC_Go_BP_down <- enrichGO(gene = CD_PC_DEG_down, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                           ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
Podo_Go_BP_down <- enrichGO(gene = Podo_DEG_down, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                          ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)

DLH_Go_BP_down_5 <- DLH_Go_BP_down@result[c(2,3,4,8,9),c("Description","pvalue")]
DLH_Go_BP_down_5$log10Pvalue <- -log10(DLH_Go_BP_down_5$pvalue)
DLH_Go_BP_down_5$subtype <- "DLH"

ALH_Go_BP_down_5 <- ALH_Go_BP_down@result[c(2,3,4,5,10),c("Description","pvalue")]
ALH_Go_BP_down_5$log10Pvalue <- -log10(ALH_Go_BP_down_5$pvalue)
ALH_Go_BP_down_5$subtype <- "ALH"

DCT_Go_BP_down_5 <- DCT_Go_BP_down@result[c(2,3,4,5,9),c("Description","pvalue")]
DCT_Go_BP_down_5$log10Pvalue <- -log10(DCT_Go_BP_down_5$pvalue)
DCT_Go_BP_down_5$subtype <- "DCT"

CD_IC_Go_BP_down_5 <- CD_IC_Go_BP_down@result[c(1:3,5,10),c("Description","pvalue")]
CD_IC_Go_BP_down_5$log10Pvalue <- -log10(CD_IC_Go_BP_down_5$pvalue)
CD_IC_Go_BP_down_5$subtype <- "CD-IC"

CD_PC_Go_BP_down_5 <- CD_PC_Go_BP_down@result[c(1:5),c("Description","pvalue")]
CD_PC_Go_BP_down_5$log10Pvalue <- -log10(CD_PC_Go_BP_down_5$pvalue)
CD_PC_Go_BP_down_5$subtype <- "CD-PC"

Podo_Go_BP_down_5 <- Podo_Go_BP_down@result[c(1:5),c("Description","pvalue")]
Podo_Go_BP_down_5$log10Pvalue <- -log10(Podo_Go_BP_down_5$pvalue)
Podo_Go_BP_down_5$subtype <- "Podo"

Go_Epi_DEG_down <- rbind(DLH_Go_BP_down_5,ALH_Go_BP_down_5,DCT_Go_BP_down_5,CD_IC_Go_BP_down_5,CD_PC_Go_BP_down_5,Podo_Go_BP_down_5)
Go_Epi_DEG_down$item <- row.names(Go_Epi_DEG_down)

Go_Epi_DEG_down$subtype <- factor(Go_Epi_DEG_down$subtype,levels = rev(c("DLH","ALH","DCT","CD-IC","CD-PC","Podo")))

ggbarplot(Go_Epi_DEG_down, x="item", y="log10Pvalue", fill = "subtype", color = "white",
          palette =  c("DLH" = '#53A85F', "ALH" = '#F1BB72',"DCT" = '#F3B1A0',
                       "CD-IC" = '#D6E7A3',"CD-PC" = '#57C3F3',"Podo" = "#E59CC4"),
          sort.val = "asc", 
          sort.by.grodowns=TRUE,
          x.text.angle=0, 
          xlab = NULL) + coord_flip() 

table(Go_Epi_DEG_down$Description)

write.csv(Go_Epi_DEG_down,"Go_Epi_DEG_down.csv")

# 6.8 Other epithelial gene set variance analysis (GSVA) 

## Species gene ID conversion
library(tidyverse)
library(biomaRt) # BiocManager::install("biomaRt")
library(GSVA) # BiocManager::install("GSVA")

listMarts() # available datasets

#        biomart                version
# 1 ENSEMBL_MART_ENSEMBL    Ensembl Genes 104
# 2 ENSEMBL_MART_MOUSE      Mouse strains 104

listDatasets(human) #
listDatasets(mouse) 
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"), verbose = FALSE) 
human <- useMart('ensembl',dataset = "hsapiens_gene_ensembl") 
mouse <- useMart('ensembl',dataset = "mmusculus_gene_ensembl") 

mouse.gene <- row.names(expr) # 

m2h.g <- getLDS(attributes = c("mgi_symbol"),
                filters = "mgi_symbol", 
                values = mouse.gene,  
                mart = mouse, 
                attributesL = c("hgnc_symbol","chromosome_name","start_position"),
                martL = human, 
                uniqueRows = T)

m.g <- m2h.g$MGI.symbol 
h.g <- m2h.g$HGNC.symbol

expr <- as.data.frame(expr)
expr_sub <- expr[m.g,]
expr_sub <- expr_sub[nchar(h.g) != 0,]
expr_sub <- expr_sub[!duplicated(h.g[nchar(h.g) != 0]),]
row.names(expr_sub) <- h.g[nchar(h.g) != 0][!duplicated(h.g[nchar(h.g) != 0])]
expr_sub <- as.matrix(expr_sub)

saveRDS(expr_sub,file="Other epithelial_expr.rds")

## perform GAVA for PT cells

Other_epi = readRDS("SCT_Other_epi.rds")
# table(ST@meta.data$orig.ident)
# PT@meta.data$cell_type = ST@active.ident

meta <- Other_epi@meta.data[,c("seurat_clusters","cell_type")]
meta

expr <- as.matrix(PT@assays$SCT@counts)

gmt2list <- function(gmtfile){ 
  sets <- as.list(read_lines(gmtfile)) 
  for(i in 1:length(sets)){ 
    tmp = str_split(sets[[i]], '\t') 
    n = length(tmp[[1]]) 
    names(sets)[i] = tmp[[1]][1] 
    sets[[i]] = tmp[[1]][3:n] 
    rm(tmp, n) 
  } 
  return(sets) 
}

 s.sets  <-  gmt2list("mouse_HMouse.gmt")

es.matrix = gsva(expr, s.sets, kcdf = "Poisson",method = "gsva", parallel.sz = 50)

write.csv(es.matrix,"Other_epi_gsva.csv")


# 6.9 Other_epi cell cycle phase analysis

Other_EPI <- readRDS("SCT_PT_anno.rds")
DimPlot(Other_EPI)

View(cc.genes)

g2m_genes <- cc.genes$g2m.genes 
s_genes <- cc.genes$s.genes 

## Transform the cell cycle gene set
convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

m.g2m.genes <- convertHumanGeneList(cc.genes$g2m.genes)
m.s.genes <- convertHumanGeneList(cc.genes$s.genes)

## Use CellCycleScoring function to perform cell cycle scoring
Other_EPI
DefaultAssay(Other_EPI) <- "SCT"
Other_EPI <- CellCycleScoring(Other_EPI, g2m.features=m.g2m.genes, s.features=m.s.genes)

colnames(Other_EPI@meta.data)
table(Other_EPI$Phase)
# G1  G2M    S 
# 7151 3948 6285

DimPlot(Other_EPI,label = F,group.by = 'Phase') + theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title =  paste0("Cell cycle phase")) + scale_color_npg()
DimPlot(Other_EPI,label = F,group.by = 'Phase',split.by = "type") + theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title =  paste0("Cell cycle phase"))

saveRDSOther_EPI,"SCT_Other_EPI_anno+CC.rds")


## draw bar plot
dfsam2 <- as.data.frame(table(PT$Phase,PT$cell_type_sub,PT$type))
dfsam2

ggplot(data = dfsam2,aes(x = Var2,y=Freq,fill = Var1)) +
  geom_bar(stat="identity",position = "fill",width = 0.8,color='black')+
  scale_y_continuous(labels = scales::percent)+
  theme_classic() + 
  labs(y = 'Fraction of different PT subtype (%)',x="") +
  # coord_flip() + 
  scale_fill_jama() +
  facet_grid(Var3 ~ .) 
# NoLegend() # 6x5





