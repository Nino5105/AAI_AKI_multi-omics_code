
# 1.Rawdata_process

# 1.1 library packages
library(dplyr)
library(Seurat)
library(patchwork)
library(clustree)
library(ggsci)
library(ggplot2) 
options(future.globals.maxSize = 20000 * 1024^2)

# 1.2 load 10X genomics datasets
KA1 <- Read10X(data.dir = "KA1_filtered_feature_bc_matrix/")
KA2 <- Read10X(data.dir = "KA2_filtered_feature_bc_matrix/")
KA3 <- Read10X(data.dir = "KA3_filtered_feature_bc_matrix/")
KC1 <- Read10X(data.dir = "KC1_filtered_feature_bc_matrix/")
KC2 <- Read10X(data.dir = "KC2_filtered_feature_bc_matrix/")
KC3 <- Read10X(data.dir = "KC3_filtered_feature_bc_matrix/")

KA1 <- CreateSeuratObject(counts = KA1, project = "KA1", min.cells = 3, min.features = 200)
KA2 <- CreateSeuratObject(counts = KA2, project = "KA2", min.cells = 3, min.features = 200)
KA3 <- CreateSeuratObject(counts = KA3, project = "KA3", min.cells = 3, min.features = 200)
KC1 <- CreateSeuratObject(counts = KC1, project = "KC1", min.cells = 3, min.features = 200)
KC2 <- CreateSeuratObject(counts = KC2, project = "KC2", min.cells = 3, min.features = 200)
KC3 <- CreateSeuratObject(counts = KC3, project = "KC3", min.cells = 3, min.features = 200)

KA1 # 20198 features across 11708 samples
KA1$sample <- "KA1"
KA1$type <- "Treatment"

KA2 # 19991 features across 11046 samples
KA2$sample <- "KA2"
KA2$type <- "Treatment"

KA3 # 20144 features across 11527 samples
KA3$sample <- "KA3"
KA3$type <- "Treatment"

KC1 # 19037 features across 11121 samples
KC1$sample <- "KC1"
KC1$type <- "Control"

KC2 # 18914 features across 9400 samples
KC2$sample <- "KC2"
KC2$type <- "Control"

KC3 # 19308 features across 12437 samples
KC3$sample <- "KC3"
KC3$type <- "Control"

# 1.3 Evaluate mitochondrial genes percentage
KA1[["percent.mt"]] <- PercentageFeatureSet(KA1, pattern = "mt-") 
KA2[["percent.mt"]] <- PercentageFeatureSet(KA2, pattern = "mt-") 
KA3[["percent.mt"]] <- PercentageFeatureSet(KA3, pattern = "mt-") 
KC1[["percent.mt"]] <- PercentageFeatureSet(KC1, pattern = "mt-") 
KC2[["percent.mt"]] <- PercentageFeatureSet(KC2, pattern = "mt-") 
KC3[["percent.mt"]] <- PercentageFeatureSet(KC3, pattern = "mt-")

# 1.4 Draw violin plot: detect the number of genes, reads, and mitochondrial  genes percentage
VlnPlot(KA1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # 12708 samples
VlnPlot(KA2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # 11046 samples
VlnPlot(KA3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # 11527 samples
VlnPlot(KC1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # 11121 samples
VlnPlot(KC2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # 9400 samples
VlnPlot(KC3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # 12437 samples

# 1.5 Filter the low quality cell 
KA1 <- subset(KA1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA > 500 & nCount_RNA < 50000 & percent.mt < 25)
KA2 <- subset(KA2, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA > 500 & nCount_RNA < 50000 & percent.mt < 25)
KA3 <- subset(KA3, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA > 500 & nCount_RNA < 50000 & percent.mt < 25)
KC1 <- subset(KC1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA > 500 & nCount_RNA < 50000 & percent.mt < 30)
KC2 <- subset(KC2, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA > 500 & nCount_RNA < 50000 & percent.mt < 30)
KC3 <- subset(KC3, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA > 500 & nCount_RNA < 50000 & percent.mt < 30)

# 1.6 Draw a violin plot: detect the number of genes, reads, and mitochondrial  genes percentage
VlnPlot(KA1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # 10068 samples
VlnPlot(KA2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # 9326 samples
VlnPlot(KA3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # 9561 samples
VlnPlot(KC1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # 8281 samples
VlnPlot(KC2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # 6272 samples
VlnPlot(KC3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # 8703 samples

all <- merge(KA1,c(KA2,KA3,KC1,KC2,KC3))
VlnPlot(all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0) # 8703 samples

# 1.7 Use SCT and CCA to integrate multiple samples
KA1 <- SCTransform(KA1, return.only.var.genes = FALSE,assay = "RNA", verbose = FALSE)
KA2 <- SCTransform(KA2, return.only.var.genes = FALSE,assay = "RNA", verbose = FALSE)
KA3 <- SCTransform(KA3, return.only.var.genes = FALSE,assay = "RNA", verbose = FALSE)
KC1 <- SCTransform(KC1, return.only.var.genes = FALSE,assay = "RNA", verbose = FALSE)
KC2 <- SCTransform(KC2, return.only.var.genes = FALSE,assay = "RNA", verbose = FALSE)
KC3 <- SCTransform(KC3, return.only.var.genes = FALSE,assay = "RNA", verbose = FALSE)

object_list = list(KA1,KA2,KA3,KC1,KC2,KC3)
selfeatures <- SelectIntegrationFeatures(object.list = object_list, nfeatures = 3000)
scc.list <- PrepSCTIntegration(object.list = object_list, anchor.features = selfeatures, verbose = FALSE)
scc.anchors <- FindIntegrationAnchors(object.list = scc.list, normalization.method = "SCT",anchor.features = selfeatures, verbose = FALSE)
scc_integrated <- IntegrateData(anchorset = scc.anchors, normalization.method = "SCT",verbose = FALSE)

remove(object_list,scc.anchors,scc.list)
scc_integrated # 45152 features across 52211 samples within 3 assays

head(scc_integrated@meta.data) 
scc_integrated@assays$RNA@data 

VlnPlot(scc_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
saveRDS(scc_integrated,file="SCT_combine.sample.rds")
saveRDS(c(KA1,KA2,KA3,KC1,KC2,KC3),file="SCT_each.sample.rds")

