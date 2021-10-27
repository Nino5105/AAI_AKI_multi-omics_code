
# 6.Proximal tubules cell

# 6.1 library packages
library(Seurat)
library(tidyverse)
library(patchwork)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggsci)

PT <- readRDS("SCT_PT.rds")
DimPlot(PT)

# 6.2 datasets dimension reduction and re-cluster

PT <- RunPCA(object = PT,verbose = FALSE)
ElbowPlot(PT,ndims = 50) 

PT <- FindNeighbors(PT,dim=1:15) # set PCA number = 15
PT <- FindClusters(PT,resolution = 0.8)
PT <- RunUMAP (PT,reduction="pca", dims = 1:15)
DimPlot(PT)

PT@assays

# 6.3 marker gene expression
DimPlot(PT,label = T,pt.size = 1,label.size = 5)

DoHeatmap(PT, slot = "scale.data",features = c(
  "Slc17a3","Slc22a8","Fxyd2", # proximal tubule(PT) 
  "Slc5a2","Slc5a12", # Proximal convoluted tubule cell(PTC)
  "Atp11a", "Slc13a3", # proximal straight tubules(PST)
  "Krt8","Krt18","Cd24","Vcam1")) + NoLegend()

# 6.4 cell type annotation

head(PT@active.ident)
new.cluster.ids <- c("PT-S","PCT","PST","PST","PCT","PT-S","PCT","PT-S",
                     "PCT","PT-S","PT-S","PST","PCT","PST","PCT","PT-S","PT-S","PCT")

names(new.cluster.ids) <- levels(PT)
PT <- RenameIdents(PT, new.cluster.ids)
levels(PT) <- c("PT-S","PCT","PST")

p1 <- DimPlot(PT,label = TRUE,reduction = "umap",label.size = 5) + scale_color_npg() + theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title = "PT cells")
p2 <- DimPlot(PT,label = F,group.by = "type",label.size = 5,cols = c('Control' = '#00BFC4', 'Treatment' = '#F8766D')) + theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title =  paste0("Dataset type"))
p3 <- DimPlot(PT,label = TRUE,reduction = "umap",label.size = 5,split.by = "type") + scale_color_npg()

p1 + p2
p3

# 6.4 save metadata and PT.rds

PT$cell_type_sub <- PT@active.ident
write.csv(PT@meta.data,"PT_meta_data.csv")
saveRDS(PT,"3.Result/6.1 PT/PT_anno.rds")

# 6.5 count the proportions of each group and each sample in different subclusters

dfsam1 <- as.data.frame(table(PT$type,PT$cell_type_sub))
dfsam1$Var1 <- factor(dfsam1$Var1,levels = rev(c("Control","Treatment")))
dfsam1$Var2 <- factor(dfsam1$Var2,levels = levels(dfsam1$Var2))

p4 <- ggplot(data = dfsam1,aes(x = Var2,y=Freq,fill = Var1)) +
  geom_bar(stat="identity",position = "fill",width = 0.7)+
  scale_y_continuous(labels = scales::percent)+
  scale_fill_manual(values = rev(c("#00BFC4","#F8766D")), name = "" )+
  theme_classic() + 
  labs(y = 'Fraction of different type',x="") + 
  # coord_flip()  + 
  NoLegend()

dfsam2 <- as.data.frame(table(PT$sample,PT$cell_type_sub))
dfsam2$Var1 <- factor(dfsam2$Var1,levels = rev(c("KC1","KC2","KC3","KA1","KA2","KA3")))
dfsam2$Var2 <- factor(dfsam2$Var2,levels = levels(dfsam2$Var2))

p5 <- ggplot(data = dfsam2,aes(x = Var2,y=Freq,fill = Var1)) +
  geom_bar(stat="identity",position = "fill",width = 0.7)+
  scale_y_continuous(labels = scales::percent)+
  theme_classic() + 
  labs(y = 'Fraction of different sample',x="") +
  # coord_flip() + 
  scale_fill_nejm() + 
  NoLegend()


p4/p5

# 6.6 draw violin plot

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
  '#E64933', '#4DBBD5','#009C81')

feature = c(
  # "Slc27a2","Lrp2",
  "Fxyd2","Gpx3",
  "Slc5a2","Slc5a12",
  "Atp11a", "Slc13a3"
)

StackedVlnPlot(PT, c(feature), pt.size=0, cols=my36colors)

# 6.7 PT (including PT-S, PCT and PST) different expressed gene (DEGs)

## DEGs analysis（AA treatment vs Control）

PT_DEG <- FindMarkers(PT,ident.1 = "Treatment",ident.2 = "Control",
                      min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
head(PT_DEG)

PT_DEG$change = ifelse(PT_DEG$p_val_adj < 0.05 & abs(PT_DEG$avg_log2FC) >= 0.25, 
                       ifelse(PT_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(PT_DEG$change)
write.csv(PT_DEG,"PT_DEG.csv")

PT_DEG_up <- row.names(PT_DEG[PT_DEG$avg_log2FC>0 & PT_DEG$p_val_adj < 0.05,]) 
head(PT_DEG_up)

PT_DEG_down <- row.names(PT_DEG[PT_DEG$avg_log2FC<0 & PT_DEG$p_val_adj < 0.05,]) 
head(PT_DEG_down)

## Extract the average value of the expression of each sample and draw heatmap

PT@active.ident <- as.factor(PT$sample)
DimPlot(PT)
PT_averge <- AverageExpression(PT)
head(PT_averge$SCT)
PT_exper <- PT_averge$SCT[c(PT_DEG_up,PT_DEG_down),]

pdf("pheatmap in PT.pdf",width = 5,height = 10)

annotation_col <- data.frame(Sample = factor(c("KA1","KA2","KA3","KC1","KC2","KC3")),Type = factor(c(rep("Treatment",3),rep("Control",3))))
row.names(annotation_col) <- c("KA1","KA2","KA3","KC1","KC2","KC3")
ann_colors = list(Sample = c(KA1="#BC3C29",KA2="#0072B5",KA3="#E18727",KC1="#20854E",KC2="#7876B1",KC3="#6F99AD" ), CellType = c(Control = "#00BFC4",Treatment = "#F8766D"))

pheatmap::pheatmap(PT_exper,show_colnames = T,show_rownames = F,
                   cellwidth = 10,cellheight = 0.2,scale = "row",
                   treeheight_row = 0,treeheight_col = 0,
                   cluster_rows = F,cluster_cols = F,
                   annotation_col = annotation_col,annotation_colors = ann_colors)
dev.off()


## GO enrichment analysis of up- and down-regulated genes

Go_BP_up <- enrichGO(gene = PT_DEG_up, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                     ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)

Go_BP_down <- enrichGO(gene = PT_DEG_down, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                       ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)

p1 <- dotplot(Go_BP_up, showCategory = 10,title = "The GO enrichment(BP) analysis of Treatment vs Control up DEG in PT")
p2 <- dotplot(Go_BP_down, showCategory = 10,title = "The GO enrichment(BP) analysis of Treatment vs Control down DEG in PT")

pdf("GO enrichment analysis of Treatment vs Control DEG in PT.pdf",width = 12,height = 10)

p1/p2

dev.off()

write.csv(Go_BP_up@result,"The GO enrichment(BP) analysis of Treatment vs Control up DEG in PT.csv" )
write.csv(Go_BP_down@result,"The GO enrichment(BP) analysis of Treatment vs Control down DEG in PT.csv" )

# 6.8 PT (including PT-S, PCT and PST) gene set variance analysis (GSVA) 

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

saveRDS(expr_sub,file="PT_expr.rds")

## perform GAVA for PT cells

PT = readRDS("SCT_PT.rds")
# table(ST@meta.data$orig.ident)
# PT@meta.data$cell_type = ST@active.ident

meta <- PT@meta.data[,c("seurat_clusters","cell_type")]
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

write.csv(es.matrix,"PT_gsva.csv")


gsva_value <- read.csv("1.Data/gsva_Hallmarks_mouse/PT_gsva.csv",check.names = F,row.names = 1)
row.names(gsva_value) <- substring(row.names(gsva_value),10,)

GSVA_overlap_pathway <- read.csv("GSVA_overlap_pathway.csv",row.names = 1)
gsva_value <- gsva_value[row.names(GSVA_overlap_pathway),]

gsva_value_t <- as.data.frame(t(gsva_value))

meta_data <- read.csv("1.Data/PT_meta_data.csv",check.names = F,row.names = 1)[,c("sample","type","cell_type_sub")]
meta_data$id <- row.names(meta_data)

## draw vlnplot for gava 

 i = 0
 for (i in 1:50){
   gsva_sig <- as.data.frame(gsva_value_t[,colnames(gsva_value_t)[i]],row.names = row.names(gsva_value_t))
   gsva_sig$id <- row.names(gsva_sig)
   gsva_vln_data <- merge(meta_data,gsva_sig,by = "id")
   colnames(gsva_vln_data)[5] <- "value"
   
   gsva_vln_data_summary <- summarySE(gsva_vln_data, measurevar="value", groupvars=c("type","cell_type_sub"))
   
   pdf(paste0(i,"_PT_",colnames(gsva_value_t)[i],"_GSVA_all.pdf"),width = 5,height = 4)
   
   vln_plot <- ggplot(data=gsva_vln_data, aes(x=cell_type_sub, y=value,fill=type)) + 
     geom_split_violin(trim=T,color="black") +
     geom_errorbar(data = gsva_vln_data_summary,aes(ymin = value-sd, ymax=value+sd), 
                width=0.1, 
                position=position_dodge(0.4), 
                color="black",
                alpha = 0.7,
                size=0.5) + 
     scale_fill_manual(values = c("#00BFC4","#F8766D")) +
     theme_bw()+
     theme(axis.text.y=element_text(size=10,face="bold",colour="black"), 
        axis.title.y=element_text(size = 12,face="bold",colour="black"), 
        axis.text.x=element_text(size=12,face="bold",colour="black"), 
        axis.title.x=element_text(size = 12,face="bold",colour="black"),
        legend.title=element_text(size=12,face="bold",colour="black"), 
        legend.text=element_text(size=12,face="bold",colour="black"), 
        title = element_text(size = 13,face="bold",colour="black"),
        # panel.border = element_blank(),axis.line = element_line(colour = "black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+  
     ylab("Enrichment score (value ± sd)")+ xlab("") + 
     theme(plot.title = element_text(vjust = 1,hjust = 0.5)) + 
     labs(title =  paste0(colnames(gsva_value_t)[i]))
   print(vln_plot)
   
   dev.off()
   }
 

# 6.8 PT (including PT-S, PCT and PST) SECNIC analysis

library(SCENIC)
library(foreach)
library(pheatmap)

rm(list=ls())

## Establish an analysis environment
## dir.create("SCENIC")
## dir.create("SCENIC/int")

## load PT scRNA-Seq expression matrix and meta information
scRNA = read.csv("PT_count.csv",header = T,check.names=F,row.names = 1)
head(scRNA)[1:5,1:5]

scRNA = t(scRNA)
exprMat = as.matrix(scRNA)

## Set SCENIC configuration information
mydbDIR <- "mydbs/"
mydbs <- c("mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather",
           "mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")

names(mydbs) <- c("500bp", "10kb")

scenicOptions <- initializeScenic(org="mgi", 
                                  nCores=1,
                                  dbDir=mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "AA")


## Gene filter (the sum of gene expression>cell number*3%, and it is expressed in 1% of the cells)
genesKept <- geneFiltering(exprMat, scenicOptions, 
                           minCountsPerGene = 3 * 0.01 * ncol(exprMat), 
                           minSamples = ncol(exprMat) * 0.1)
exprMat_filtered <- exprMat[genesKept, ]

## Calculate the correlation matrix and log the data
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1)
 
## Run GENIE3 regression analysis to calculate the correlation of TF-Targets
runGenie3(exprMat_filtered_log, scenicOptions, nParts = 10)
save(exprMat_filtered,scenicOptions,file = "input_GENIE3_data.Rdata")

## infer the co-gene expression modules
runSCENIC_1_coexNetwork2modules(scenicOptions)

## infer the transcript regulon
runSCENIC_2_createRegulons(scenicOptions,minGenes = 1)

## Regulon activity score
runSCENIC_3_scoreCells(scenicOptions, exprMat=exprMat_filtered_log)

## Regulon active binary conversion and visualization
runSCENIC_4_aucell_binarize(scenicOptions, exprMat=exprMat_filtered_log)

## visualation of results
logMat <- exprMat_filtered_log
aucellApp <- plotTsne_AUCellApp(scenicOptions, logMat)
savedSelections <- shiny::runApp(aucellApp)

save(exprMat_filtered_log,scenicOptions,file = "input_GENIE4_data.Rdata")

## use pheatmap for visualization 
library(pheatmap)

## load the raw regulonAUC matrix
AUCmatrix <- readRDS("int/3.4_regulonAUC.Rds")
AUCmatrix <- AUCmatrix@assays@data@listData$AUC
AUCmatrix <- data.frame(t(AUCmatrix), check.names=F)
AUCmatrix <- t(AUCmatrix)
pheatmap::pheatmap(AUCmatrix,show_rownames = T,show_colnames = F)

write.csv(AUCmatrix,"AUCmatrix.csv",row.names = T)

myAUCmatrix = read.csv("v1/AUCmatrix.csv",header = T,row.names = 1,check.names = F)

## select significant regulons
my.regulons <- c(
 "Rela_extended (262g)",
 "Klf13_extended (143g)",
 "Cebpd_extended (19g)",
 "Trp53 (17g)",
 "Stat3 (77g)",
 "Fosl2 (70g)",
 "Ets2_extended (1435g)",
 "Tgif1_extended (22g)",
 "Stat1 (47g)",
 "Nfkb1 (54g)")
myAUCmatrix_sub = myAUCmatrix[my.regulons,]
myAUCmatrix_sub


ann_colors = list(sample = c(KA1="#BC3C29",KA2="#0072B5",KA3="#E18727",KC1="#20854E",KC2="#7876B1",KC3="#6F99AD" ), 
                  type = c(Treatment = "#F8766D",Control = "#00BFC4"),
                  cell_type_sub = c(PT="#E5462F",PCT="#009F86",PST="#A0DAE8"))


PT_meta <- read.csv("PT_meta_data.csv",row.names = 1,check.names = F)[,c("sample","type","cell_type_sub")]
head(PT_meta)
table(colnames(myAUCmatrix_sub) %in% row.names(PT_meta))

PT_scenic_heatmap <- pheatmap(myAUCmatrix_sub, scale = "row",
                              show_rownames = T, show_colnames = F,
                              cluster_cols = F,cluster_rows = T,
                              annotation_col = PT_meta, gaps_col = 2846,
                              annotation_colors = ann_colors,cellwidth = 0.02,cellheight = 20,
                              treeheight_row = 1,border_color = "black",color = colorRampPalette(colors = c("blue","white","red"))(100))
                              # color = colorRampPalette(rev(brewer.pal(n=5,name='PRGn')))(round(10)))
table(PT_meta$type)

PT_scenic_heatmap

write.csv(myAUCmatrix_sub,"myAUCmatrix_sub.csv") 

# 6.9 PT (including PT-S, PCT and PST) RNA velocity analysis

library(velocyto.R)
library(SeuratWrappers)
library(dplyr)
library(hash)
source("/public/group/chenjunhui/scRNA/CELL/pipe/utils/CellFunction_v3.R")

## load the merged loom files

ldat <- ReadVelocity(file = "merged.loom")

PT  <- readRDS("SCT_PT_anno.rds")
bm <- PT
bm_count_data <- as.matrix(bm@assays$RNA@counts)

## read sample infomation
file=read.csv("sample.csv",header=F)
samH <- hash(file$V1,file$V2)

## Create a matrix of the same size
tmp_zero <- matrix(0, nrow=dim(bm_count_data)[1], ncol=dim(bm_count_data)[2])
colnames(tmp_zero) <- colnames(bm_count_data)
rownames(tmp_zero) <- rownames(bm_count_data)

DATA_MK <- function(data) {
    ldat_spliced_data <- dmatrix(data)
	ldat_spliced_data <-do.call(cbind,ldat_spliced_data)
	colnames(ldat_spliced_data) <- gsub("x$", "",colnames(ldat_spliced_data))
	sample<-as.data.frame(limma::strsplit2(split=":",colnames(ldat_spliced_data)))
	sample$V3 <-"SAM"
	sample$V3[which(sample$V1 == "KA1")] <- paste(sample$V2[which(sample$V1 == "KA1")],"-1_",samH[["KA1"]],sep="")
	sample$V3[which(sample$V1 == "KA2")] <- paste(sample$V2[which(sample$V1 == "KA2")],"-1_",samH[["KA2"]],sep="")
	sample$V3[which(sample$V1 == "KA3")] <- paste(sample$V2[which(sample$V1 == "KA3")],"-1_",samH[["KA3"]],sep="")
	sample$V3[which(sample$V1 == "KC1")] <- paste(sample$V2[which(sample$V1 == "KC1")],"-1_",samH[["KC1"]],sep="")
	sample$V3[which(sample$V1 == "KC2")] <- paste(sample$V2[which(sample$V1 == "KC2")],"-1_",samH[["KC2"]],sep="")
	sample$V3[which(sample$V1 == "KC3")] <- paste(sample$V2[which(sample$V1 == "KC3")],"-1_",samH[["KC3"]],sep="")
	colnames(ldat_spliced_data) <- sample$V3
        index <- duplicated(rownames(ldat_spliced_data))
        ldat_spliced_data_dup <- ldat_spliced_data[!index,]
	tmp<-dmatrix(ldat_spliced_data_dup)
	ldat_spliced_data_dup<-do.call(cbind,tmp)
	return(ldat_spliced_data_dup)
}

spliced_data <- DATA_MK(ldat$spliced)
unspliced_data <- DATA_MK(ldat$unspliced)

MAERGE_MATIX <- function(tmp_zero, ldat_spliced_data_dup) {
    ldat_spliced_data_dup <- ldat_spliced_data_dup[,colnames(ldat_spliced_data_dup) %in% colnames(tmp_zero)]
    ldat_spliced_data_dup <- ldat_spliced_data_dup[rownames(ldat_spliced_data_dup) %in% rownames(tmp_zero),]
    TURE_FAlES <- !colnames(tmp_zero) %in% colnames(ldat_spliced_data_dup)
    tmp_col <- tmp_zero[,TURE_FAlES]
    colnames(tmp_col) <- colnames(tmp_zero)[TURE_FAlES]
    if (length(colnames(tmp_zero)[TURE_FAlES])>0) {
        for (i in length(colnames(tmp_col))) {
            add_col <- matrix(0, nrow=dim(ldat_spliced_data_dup)[1], ncol=1)
            colnames(add_col) <- colnames(tmp_col)[i]
            rownames(add_col) <- rownames(ldat_spliced_data_dup)
            ldat_spliced_data_dup <- cbind(ldat_spliced_data_dup, add_col)
        }
    }

    TURE_FAlES <- !rownames(tmp_zero) %in% rownames(ldat_spliced_data_dup)
    tmp_row <- tmp_zero[TURE_FAlES,]
    rownames(tmp_row) <- rownames(tmp_zero)[TURE_FAlES]
    if (length(rownames(tmp_zero)[TURE_FAlES])>0) {
        for (j in 1:length(rownames(tmp_row))) {
            add_row <- matrix(0, nrow=1, ncol=dim(ldat_spliced_data_dup)[2])
            colnames(add_row) <- colnames(ldat_spliced_data_dup)
            rownames(add_row) <- rownames(tmp_row)[j]
            ldat_spliced_data_dup <- rbind(ldat_spliced_data_dup, add_row)
        }
    }
    return(ldat_spliced_data_dup)
}

tmp_zero_spliced <- MAERGE_MATIX(tmp_zero,spliced_data)
tmp_zero_unspliced <- MAERGE_MATIX(tmp_zero,unspliced_data)

bm[['spliced']] <- CreateAssayObject(counts = tmp_zero_spliced)
bm[['unspliced']] <- CreateAssayObject(counts = tmp_zero_unspliced)

bm <- RunVelocity(object = bm, deltaT = 1, kCells = 25, fit.quantile = 0.02)

saveRDS(bm, file="./RunVelocity.rds")

ident.colors <- c("#E5462F","#A0DAE8","#009F86")
names(x = ident.colors) <- levels(x = bm)
cell.colors <- ident.colors[Idents(object = bm)]
names(x = cell.colors) <- colnames(x = bm)

pdf_name = paste("./embedding.cor.pdf", sep="")
convert_line = paste("/usr/bin/convert ", pdf_name, png_name, sep="\t")
pdf(pdf_name)
show.velocity.on.embedding.cor(emb = Embeddings(object = bm, reduction = "umap"), vel = Tool(object = bm, slot = "RunVelocity"), 
n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, 
min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, do.par = FALSE, cell.border.alpha = 0.1)
dev.off()


# 6.10 PT (including PT-S, PCT and PST) cell cycle phase analysis

PT <- readRDS("SCT_PT_anno.rds")
DimPlot(PT)

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
PT
DefaultAssay(PT) <- "SCT"
PT <- CellCycleScoring(PT, g2m.features=m.g2m.genes, s.features=m.s.genes)

colnames(PT@meta.data)
table(PT$Phase)
# G1  G2M    S 
# 7151 3948 6285

DimPlot(PT,label = F,group.by = 'Phase') + theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title =  paste0("Cell cycle phase")) + scale_color_npg()
DimPlot(PT,label = F,group.by = 'Phase',split.by = "type") + theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title =  paste0("Cell cycle phase"))

saveRDS(PT,"1.Data/5.1 SCT_PT_anno+CC.rds")


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