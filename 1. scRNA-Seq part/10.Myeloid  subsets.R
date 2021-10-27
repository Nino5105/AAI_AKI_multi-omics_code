# 9.Myelodi subsets cell

# 10.1 library packages

library(dplyr)
library(Seurat)
library(patchwork)
library(clustree)
library(ggsci)
library(ggplot2) 

Myeloid <- subset(scc_integrated,idents = "Myeloid")
Myeloid$cell_type_sub <- "Myeloid"
DimPlot(Myeloid,label = T,pt.size = 1,label.size = 5,group.by = "integrated_snn_res.0.8")

# 10.2 datasets dimension reduction and re-cluster

Myeloid <- RunPCA(object = Myeloid,verbose = FALSE)
ElbowPlot(Myeloid,ndims = 20) # 选择主成分为15

Myeloid <- FindNeighbors(Myeloid,dim=1:20)
Myeloid <- FindClusters(Myeloid,resolution = 0.8)
Myeloid <- RunUMAP(Myeloid,reduction="pca", dims = 1:20)
DimPlot(Myeloid,label = T,pt.size = 1,label.size = 5)

# 10.3 marker gene expression

DoHeatmap(Myeloid, slot = "scale.data",features = c(
  "C1qa","C1qb","Apoe","Cd74",# Macrophage
  "Cd80","Cd86","Cd68","Tlr2","Tnf","Socs3", # Macrophage M1 markers
  "Il1b","!Il2","!Il6","Il12b","!Il23a","!Csf1",
  "Ccl2","Ccl3","Ccl5","Cxcl9","Cxcl10","Cxcl16", # Macrophage M1 cytokines
  "Cd163","Msr1","Mrc1","!Cd200r1","!Tgm2","Arg1","Il10", # Macrophage M2 markers
  "!Tgfb1","Ccl5","Ccl8","!Ccl17","!Ccl22","!Ccl24","Cd36", # Macrophage M2 cytokines
  "Cx3cr1","Cd209d","Clec10a", # M2-like
  "Mki67","Cdca3", # Macrophage prliferation
  "!Cst3","!Cpa3","!Cma1","!Tpsb2","!Cd34","Enpp3","!Il2ra","Cd2","Kit" # MAST cell
)) + NoLegend()


# 10.4 cell type annotation

Myeloid@active.ident <- Myeloid$integrated_snn_res.0.8
head(Myeloid@active.ident)
DimPlot(Myeloid,label = T,group.by = "integrated_snn_res.0.8",label.size = 5,pt.size = 1)

new.cluster.ids <- c("Macro M2",
                     "Macro M1",
                     "Macro M2",
                     "Monocytes",
                     "Macro M1",
                     "Macro M2",
                     "Mast cell",
                     "Macro M1",
                     "Macro M1",
                     "Macro M1",
                     "Macro M1",
                     "Macro Pro",
                     "Macro Pro",
                     "Macro Pro")

names(new.cluster.ids) <- levels(Myeloid)
Myeloid <- RenameIdents(Myeloid, new.cluster.ids)

levels(Myeloid) <- c("Macro M1",
                   "Macro M2",
                   "Macro Pro",
                   "Monocytes",
                   "Mast cell")

Myeloid$cell_type_sub <- Myeloid@active.ident

p1 <- DimPlot(Myeloid,label = TRUE,reduction = "umap",label.size = 5,pt.size = 1) + scale_color_d3() + theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title = "Myeloid immune cell")
p1

p2 <- DimPlot(Myeloid,label = F,group.by = "type",label.size = 5,pt.size = 1,cols = c('Control' = '#00BFC4', 'Treatment' = '#F8766D')) + theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title =  paste0("Dataset type"))
p2

p3 <- DimPlot(Myeloid,label = TRUE,reduction = "umap",label.size = 5,,pt.size = 1,split.by = "type") + scale_color_d3()
p3

# DimPlot(Other_immune,label = F,group.by = "sample",label.size = 5)+ theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title =  paste0("Sample id")) + scale_color_nejm()
# DimPlot(Other_immune,label = T,group.by = "id",label.size = 5) + theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title =  paste0("Cluster id"))

# 10.5 draw heatmap plot

AverageExpression_value <-  AverageExpression(Myeloid, assays = "SCT", features = c(
  "C1qa","C1qb","Apoe","Cd74",# Macrophage
  "Cd80","Cd86","Cd68","Tlr2","Tnf","Socs3", # Macrophage M1 markers
  "Il1b","!Il2","!Il6","Il12b","!Il23a","!Csf1","Nos2",
  "Ccl2","Ccl3","Ccl5","Cxcl9","Cxcl10","Cxcl16", # Macrophage M1 cytokines
  
  "Cd163","Mrc1","!Msr1","!Cd200r1","!Tgm2","!Arg1","Il10", # Macrophage M2 markers
  "!Tgfb1","Ccl5","Ccl8","Ccl17","!Ccl18","!Ccl22","!Ccl24","!Cd36", # Macrophage M2 cytokines
  "!H2-Ab1","Adgre1","Cx3cr1","Cd209d","Clec10a", # M2-like
  # "Itgam",,"Csf1r","Ly6c1",
  "Mki67","Cdca3", # Macrophage prliferation
   "Lyz1","Cd14",   # Monocytes
  "!Cst3","!Cpa3","!Cma1","!Tpsb2","!Cd34","Enpp3","!Il2ra","Cd2","Kit", # MAST cell
"!Ppp2ca"), return.seurat = FALSE, group.by = "ident")

pheatmap::pheatmap(AverageExpression_value$SCT,cluster_rows = F,cluster_cols = F,
                   border_color = "black",
                   angle_col = 45,gaps_row = c(4,7,19,21,28,30,32),#gaps_col = c(4,7,8),
                   cellwidth = 30,cellheight = 15,
                   treeheight_col = 4,scale = "row")


# 10.6 save metadata and myeloid.rds

Myeloid$cell_type_sub <- Macro@active.ident
write.csv(Myeloid@meta.data,"3.Result/6.3 Myeloid/meta_data.csv")
saveRDS(Myeloid,"1.Data/5.3 Myeloid_anno.rds")

# 10.7 count the proportions of each group and each sample in different subclusters

table(Myeloid$type,Macro@active.ident)

dfsam1 <- as.data.frame(table(Myeloid$type,Macro@active.ident))
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

dfsam2 <- as.data.frame(table(Myeloid$sample,Macro$cell_type_sub))
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

# 10.8 draw feature plot & vlnplot

p1 <- FeaturePlot(Myeloid, features = c("Cd86","Il1b","Mrc1","Ccl8","Mki67","Lyz1","Kit"),ncol = 4,slot = "scale.data",cols = c("white","red"),label=F,pt.size = 0.2)+ 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)
p1

Myeloid <- readRDS("Myeloid_anno.rds")
DimPlot(Myeloid)
Macro <- subset(Myeloid,idents = c("Macro M1","Macro M2","Macro Pro"))

p1 <- VlnPlot(Macro, features = c("Il1b","Tnf","Tgfb1","Il10"),
              ncol = 1,
              cols = c("#00BFC4","#F8766D"),
              # group.by = "type",
              split.by = "type",
              y.max = 5,
              # stack = T,
              assay = "SCT",
              slot = "scale.data",
              # flip = 90,
              pt.size = 0)
p1

p2 <- VlnPlot(Macro, features = c("Arg1"),
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

# 10.9 macrophage cell trajectory analysis and GO enrichment analysis

Myeloid <- readRDS("1.Data/5.3 Myeloid_anno.rds")
table(Myeloid@active.ident)
Macro <- subset(Myeloid,idents = c("Macro M1","Macro M2","Macro Pro"))

DefaultAssay(Macro) <- "RNA"
DimPlot(Macro)

expr_matrix <- as(as.matrix(Macro@assays$RNA@counts), 'sparseMatrix')

colnames(Macro@meta.data)
p_data <- Macro@meta.data[,c(4,5,18)]

f_data <- data.frame(gene_short_name = row.names(Macro),row.names = row.names(Macro))

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

expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 10))
length(expressed_genes)

diff <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr="~cell_type_sub",cores=10)

head(diff)

deg <- subset(diff, qval < 0.001)
table(diff$qval <= 0.001)

deg <- deg[order(deg$qval,decreasing=F),]
head(deg)

write.csv(deg,file="monocle_Myeloid_DEG.csv",quote=F)

ordergene <- rownames(deg) 
cds <- setOrderingFilter(cds, ordergene)  
table(cds@featureData@data[["use_for_ordering"]]) 
plot_ordering_genes(cds) 

cds <- reduceDimension(cds, max_components = 2,method = 'DDRTree')

cds <- orderCells(cds)

p0 <- plot_cell_trajectory(cds,color_by="Pseudotime",show_tree = T,cell_size = 1,cell_link_size = 2,show_backbone=TRUE) 
p1 <- plot_cell_trajectory(cds,color_by="State",show_tree = T,cell_size = 1,cell_link_size = 2,show_backbone=TRUE) +
  scale_color_manual(breaks = c("1", "2", "3"), values=c("#979797","#EF4F5B","#7990C8"))
p2 <- plot_cell_trajectory(cds,color_by="type",show_tree = T,cell_size = 1,cell_link_size = 2,show_backbone=TRUE) + scale_color_manual(breaks = c("Control", "Treatment"), values=c("#00BFC4", "#F8766D"))
p3 <- plot_cell_trajectory(cds,color_by="cell_type_sub",show_tree = T,cell_size = 1,cell_link_size = 2,show_backbone=F) + 
  scale_color_manual(breaks = c("Macro M1", "Macro M2", "Macro Pro"), values=c("#1F77B4","#FF7F0E","#2CA02C"))

(p0 + p1)/(p2+p3) 

plot_cell_trajectory(cds, color_by = "State")

BEAM_res <- BEAM(cds[ordergene,], branch_point = 1, cores = 10)
head(BEAM_res)

BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

table(BEAM_res$qval < 0.001)
BEAM_genes <- row.names(subset(BEAM_res,qval < 0.001)) # 选择qval< 0.001,n =3361

plot <- plot_genes_branched_heatmap(cds[BEAM_genes,],
                                    branch_point = 1,num_clusters = 5,
                                    cores = 10, 
                                    use_gene_short_name = T,
                                    show_rownames = T,
                                    return_heatmap = T) #

# ggsave(plot$ph_res, "BEAM_heatmap.pdf", width = 4, height = 5)

plot <- plot_genes_branched_heatmap(cds[BEAM_genes,],
                                    branch_point = 1,num_clusters = 5,
                                    cores = 10, 
                                    use_gene_short_name = T,
                                    show_rownames = T,
                                    return_heatmap = T) #

plot$ph_res

gene_group <- plot$annotation_row
gene_group$gene <- rownames(gene_group)
allcluster_go=data.frame()

C1_gene <- gene_group[gene_group$Cluster == 1,]$gene #n=945
C1_BP <- enrichGO(gene = C1_gene,OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
                  ont = "BP",pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
C1_BP@result[c(1:20),c("Description","pvalue")]

C2_gene <- gene_group[gene_group$Cluster == 2,]$gene #n=791
C2_BP <- enrichGO(gene = C2_gene,OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
                  ont = "BP",pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
C2_BP@result[c(1:10),c("Description","pvalue")]

C3_gene <- gene_group[gene_group$Cluster == 3,]$gene #n=464
C3_BP <- enrichGO(gene = C3_gene,OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
                  ont = "BP",pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
C3_BP@result[c(1:10),c("Description","pvalue")]

C4_gene <- gene_group[gene_group$Cluster == 4,]$gene #n=635
C4_BP <- enrichGO(gene = C4_gene,OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
                  ont = "BP",pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
C4_BP@result[c(1:10),c("Description","pvalue")]

C5_gene <- gene_group[gene_group$Cluster == 5,]$gene #n=526
C5_BP <- enrichGO(gene = C5_gene,OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
                  ont = "BP",pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
C5_BP@result[c(1:10),c("Description","pvalue")]


head(allcluster_go[,c("ID","Description","qvalue","cluster")])