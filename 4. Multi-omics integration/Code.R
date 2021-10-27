# 1.labrary packages

library(VennDiagram) 
library(ggplot2)
library(ggpubr)
library(patchwork)
library(clusterProfiler)
library(org.Mm.eg.db)

# 2.loading datasets

DEG_scRNA <- read.csv("DEG in scRNA-Seq.csv",check.names = F,row.names = 1)
DEG_RNA <- read.csv("DEG in Bulk-RNA Seq.csv",check.names = F,row.names = 1)
DEP <- read.csv("DEP in label-free.csv",check.names = F,row.names = 1)

table(DEG_scRNA$type)
DEG_scRNA_up <- rownames(DEG_scRNA[DEG_scRNA$type == "Up(3514)",])
DEG_scRNA_down <- rownames(DEG_scRNA[DEG_scRNA$type == "Down(3118)",])
DEG_scRNA_all <- rownames(DEG_scRNA[DEG_scRNA$type != "Stable(8118)",])

table(DEG_RNA$type)
DEG_RNA_up <- rownames(DEG_RNA[DEG_RNA$type == "Up(2871)",])
DEG_RNA_down <- rownames(DEG_RNA[DEG_RNA$type == "Down(1794)",])
DEG_RNA_all <- rownames(DEG_RNA[DEG_RNA$type != "Stable(9936)",])

table(DEP$change)
DEP_up <- rownames(DEP[DEP$change == "Up(1903)",])
DEP_down <- rownames(DEP[DEP$change == "Down(667)",])
DEP_all <- rownames(DEP[DEP$change != "Stable(2247)",])

overlap_up_list <- intersect(DEG_scRNA_up,intersect(DEG_RNA_up,DEP_up))

# Draw a Venn diagram of up-regulated genes

venn.plot <- 
  venn.diagram(
    x = list(
      'DEG_scRNA_up \n (n = 3514)' = DEG_scRNA_up,
      'DEG_RNA_up\n (n = 2871)' = DEG_RNA_up,
      "DEP_up\n (n = 1903)" = DEP_up),
    filename = NULL,
    col = "black",
    fill = c("#FC4E07","skyblue","#2CA02C"),
    alpha = 0.8,
    cex = 0.8,
    cat.col = 'black',
    cat.cex = 0.8,
    cat.fontface = "bold",
    margin = 0.05,
    main = "Overlap up DEP/DEG in different methods",
    main.cex = 1.2)

pdf(file="Overlap up DEP&DEG in different methods.pdf")
grid.draw(venn.plot)
dev.off()

# Draw a Venn diagram of down-regulated genes

venn.plot <- 
  venn.diagram(
    x = list(
      'DEG_scRNA_down \n (n = 3117)' = DEG_scRNA_down,
      'DEG_RNA_down\n (n = 1794)' = DEG_RNA_down,
      "DEP_down\n (n = 667)" = DEP_down),
    filename = NULL,
    col = "black",
    fill = c("#FC4E07","skyblue","#2CA02C"),
    alpha = 0.8,
    cex = 0.8,
    cat.col = 'black',
    cat.cex = 0.8,
    cat.fontface = "bold",
    margin = 0.05,
    main = "Overlap down DEP/DEG in different methods",
    # rotation.degree = 180, 
    main.cex = 1.2)

pdf(file="Overlap down DEP&DEG in different methods.pdf")
grid.draw(venn.plot)
dev.off()

# Draw correlation diagrams of different genes

## scRNA & Bulk RNA

scRNA_RNA_list <- intersect(DEG_scRNA_all,DEG_RNA_all) # n = 3238

scRNA_RNA_list_data <- merge(data.frame(DEG_scRNA[scRNA_RNA_list,c("logFC")],rep("scRNA",3238),,scRNA_RNA_list),
      data.frame(DEG_RNA[scRNA_RNA_list,c("logFC")],rep("RNA",3238),scRNA_RNA_list),by = "scRNA_RNA_list")
row.names(scRNA_RNA_list_data) <- scRNA_RNA_list_data$scRNA_RNA_list
colnames(scRNA_RNA_list_data) <- c("scRNA_RNA_list_data","scRNA_Seq","type1","Bulk_RNA_Seq","type2")
head(scRNA_RNA_list_data)

p1 <- ggscatter(data=scRNA_RNA_list_data, x = "scRNA_Seq", y = "Bulk_RNA_Seq",size = 1,
          title = "Log2(Fold change) [Treatment vs Control]",
          conf.int = T,) +
  geom_smooth(method = "lm",fill = "lightblue",size = 1.5)+
  stat_cor(method = "spearman")+ 
  xlab("scRNA-Seq") + ylab("Bulk RNA-Seq") + 
  geom_hline(yintercept=0, linetype=2, colour="black") +
  geom_vline(xintercept=0, linetype=2, colour="black") +
  xlim(-10,10) +ylim(-10,10) +theme(plot.title = element_text(hjust = 0.5))
p1

## scRNA & Protein

scRNA_Protein_list <- intersect(DEG_scRNA_all,DEP_all) # n = 1086

scRNA_Protein_list_data <- merge(data.frame(DEG_scRNA[scRNA_Protein_list,c("logFC")],rep("scRNA",1086),scRNA_Protein_list),
                             data.frame(DEP[scRNA_Protein_list,c("logFC")],rep("Protein",1086),scRNA_Protein_list),by = "scRNA_Protein_list")
row.names(scRNA_Protein_list_data) <- scRNA_Protein_list_data$scRNA_Protein_list
colnames(scRNA_Protein_list_data) <- c("scRNA_RNA_list_data","scRNA_Seq","type1","Mass_spec","type2")
head(scRNA_Protein_list_data)

p2 <- ggscatter(data=scRNA_Protein_list_data, x = "scRNA_Seq", y = "Mass_spec",size = 1,
          title = "Log2(Fold change) [Treatment vs Control]",
          conf.int = T,) +
  geom_smooth(method = "lm",fill = "lightblue",size = 1.5)+
  stat_cor(method = "spearman")+ 
  xlab("scRNA-Seq") + ylab("Mass spec (Label free)") + 
  geom_hline(yintercept=0, linetype=2, colour="black") +
  geom_vline(xintercept=0, linetype=2, colour="black") +
  xlim(-10,10) +ylim(-10,10) +theme(plot.title = element_text(hjust = 0.5))
p2

## RNA & Protein
RNA_Protein_list <- intersect(DEG_RNA_all,DEP_all) # n = 862

RNA_Protein_list_data <- merge(data.frame(DEG_RNA[RNA_Protein_list,c("logFC")],rep("RNA",862),RNA_Protein_list),
                                 data.frame(DEP[RNA_Protein_list,c("logFC")],rep("Protein",862),RNA_Protein_list),by = "RNA_Protein_list")
row.names(RNA_Protein_list_data) <- RNA_Protein_list_data$RNA_Protein_list
colnames(RNA_Protein_list_data) <- c("scRNA_RNA_list_data","RNA_Seq","type1","Mass_spec","type2")
head(RNA_Protein_list_data)

p3 <- ggscatter(data=RNA_Protein_list_data, x = "RNA_Seq", y = "Mass_spec",size = 1,
                title = "Log2(Fold change) [Treatment vs Control]",
                conf.int = T,) +
  geom_smooth(method = "lm",fill = "lightblue",size = 1.5)+
  stat_cor(method = "spearman")+ 
  xlab("Bulk RNA-Seq") + ylab("Mass spec (Label free)") + 
  geom_hline(yintercept=0, linetype=2, colour="black") +
  geom_vline(xintercept=0, linetype=2, colour="black") +
  xlim(-10,10) +ylim(-10,10) +theme(plot.title = element_text(hjust = 0.5))
p3

p1 + p2 + p3 

# Draw heatmap plot of marker genes

scRNA <- read.csv("scRNA-Seq counts.csv",row.names = 1,check.names = F)
RNA <- read.csv("Bulk RNA-Seq counts.csv",row.names = 1,check.names = F)
Mass <- read.csv("Protein abundance.csv",row.names = 1,check.names = F)

marker_list <- c(""

condition = factor(c(rep("Control",3),rep("Treatment",3)),levels = c("Control","Treatment"))
condition

ann_colors = list(type = c(Treatment = "#F8766D",Control = "#00BFC4"))

scRNA_genelist = DGEList(counts = scRNA, group = condition)
# scRNA_keep <- rowSums(cpm(scRNA_genelist)>1)>=2 
# scRNA_genelist.filted <- scRNA_genelist[scRNA_genelist,,keep.lib.sizes=FALSE] 
scRNA_genelist.norm <- calcNormFactors(scRNA_genelist)  
scRNA_logCPM <- cpm(scRNA_genelist.norm, log=TRUE, prior.count=3)
scRNA_sub <- scRNA_logCPM[marker_list,]

scRNA_sub_heatmap <- pheatmap(scRNA_sub, 
                              show_rownames = T, show_colnames = F,
                              cluster_cols = F, cluster_rows = F,
                              scale = "row",cellwidth = 20,cellheight = 20,
                              annotation_col = meta,
                              annotation_colors = ann_colors,
                              treeheight_row = 0,
                              treeheight_col = 0,
                              border_color = "black",viridis(100))


RNA_genelist = DGEList(counts = RNA, group = condition)
RNA_genelist.norm <- calcNormFactors(RNA_genelist)  
RNA_logCPM <- cpm(RNA_genelist.norm, log=TRUE, prior.count=3)
RNA_sub <- RNA_logCPM[marker_list,]
RNA_sub_heatmap <- pheatmap(RNA_sub, 
                              show_rownames = T, show_colnames = F,
                              cluster_cols = F, cluster_rows = F,
                              scale = "row",cellwidth = 20,cellheight = 20,
                              annotation_col = meta,
                              annotation_colors = ann_colors,
                              treeheight_row = 0,
                              treeheight_col = 0,
                              border_color = "black",viridis(100))

Mass_sub <- log2(Mass+1)
Mass_sub <- Mass[marker_list,]
Mass_sub_heatmap <- pheatmap(Mass_sub[c(1,3,4,8),], 
                            show_rownames = T, show_colnames = F,
                            cluster_cols = F, cluster_rows = F,
                            scale = "row",cellwidth = 20,cellheight = 20,
                            annotation_col = meta,
                            annotation_colors = ann_colors,
                            treeheight_row = 0,
                            treeheight_col = 0,
                            border_color = "black",viridis(100))

# Go enrichment analysis for overlap genes in three omics-datasets

Go_BP_up <- enrichGO(gene = overlap_up_list, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                         ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)

Go_BP_down <- enrichGO(gene = overlap_down_list, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                     ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)

dotplot(Go_BP_up, showCategory = 10,title = "The GO enrichment(BP) analysis of overlap up DEGs")

dotplot(Go_BP_down, showCategory = 10,title = "The GO enrichment(BP) analysis of overlap down DEGs")


Go_BP_up_list <- Go_BP_up@result[c(1,5,6,7,10),c("Description","pvalue")]
Go_BP_up_list$type <- "Up-regulate"

Go_BP_down_list <-Go_BP_down@result[c(1,2,3,4,5),c("Description","pvalue")]
Go_BP_down_list$type <- "Down-regulate"

Go_BP_list <- rbind(Go_BP_up_list,Go_BP_down_list)
Go_BP_list$log10Pvalue <- -log10(Go_BP_list$pvalue)


ggbarplot(Go_BP_list, x="Description", y="log10Pvalue", fill = "type", color = "white",
          palette =  c("Up-regulate" = "#ff4757", "Down-regulate" = "#546de5"),
          sort.val = "asc",#上升排序,区别于desc，具体看图演示 
          sort.by.grodowns=TRUE,#按组排序 
          x.text.angle=0, 
          xlab = NULL) + coord_flip()