
# 2.Reduction and visulations

# 2.1 library packages
library(dplyr)
library(Seurat)
library(patchwork)
library(clustree)
library(ggsci)
library(ggplot2) 
options(future.globals.maxSize = 20000 * 1024^2)

# 2.2 data reduction

# scc_integrated <-readRDS("SCT_combine.sample.rds")
scc_integrated <- RunPCA(object = scc_integrated,verbose = FALSE)
ElbowPlot(scc_integrated,ndims = 50)

# 2.3 set the PCA number

pdf("cluster.pdf")
for (i in 20:40){
  scc_integrated <- FindNeighbors(scc_integrated,dim=1:i)
  scc_integrated <- FindClusters(scc_integrated,resolution = 0.8)
  scc_integrated <- RunUMAP (scc_integrated,reduction="pca", dims = 1:i)
  p1 <- DimPlot(scc_integrated,label = TRUE,reduction = "umap") + theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 20)) + labs(title =  paste0("PCA n = ",i))
  scc_integrated <- RunTSNE(scc_integrated,dims = 1:i)
  p2 <- TSNEPlot (scc_integrated,label = TRUE) + theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 20)) + labs(title =  paste0("PCA n = ",i))
  print(p1/p2)
}
dev.off()

# 2.4 set the resolution number

scc_integrated <- FindClusters(scc_integrated,resolution = 0.1)
scc_integrated <- FindClusters(scc_integrated,resolution = 0.2)
scc_integrated <- FindClusters(scc_integrated,resolution = 0.4)
scc_integrated <- FindClusters(scc_integrated,resolution = 0.6)
scc_integrated <- FindClusters(scc_integrated,resolution = 0.8)
scc_integrated <- FindClusters(scc_integrated,resolution = 1)
clustree(scc_integrated)

scc_integrated$integrated_snn_res.0.8 

# 2.5 Draw Dimplot

DimPlot(scc_integrated,label = T,group.by = "id") + theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title =  paste0("Cluster id"))
DimPlot(scc_integrated,label = F,group.by = "type",cols = c('Control' = '#00BFC4', 'Treatment' = '#F8766D')) + theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title =  paste0("Dataset type"))
DimPlot(scc_integrated,label = F,group.by = "sample")+ theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title =  paste0("Sample id")) + scale_color_nejm()


# 2.6 Count the proportions of each group and each sample in different subclusters

dfsam <- as.data.frame(table(scc_integrated$sample,scc_integrated$id))
dfsam$Var1 <- factor(dfsam$Var1,levels = rev(c("KA1","KA2","KA3","KC1","KC2","KC3")))
dfsam$Var2 <- factor(dfsam$Var2,levels = levels(dfsam$Var2))

ggplot(data = dfsam,aes(x = Var2,y=Freq,fill = Var1)) +
  geom_bar(stat="identity",position = "fill",width = 0.5)+
  scale_y_continuous(labels = scales::percent)+
  theme_classic() + 
  labs(y = 'Fraction of different sample',x="") +
  coord_flip() + 
  scale_fill_nejm() + 
  NoLegend()

dfsam2 <- as.data.frame(table(scc_integrated$type,scc_integrated$id))
dfsam2$Var2 <- factor(dfsam2$Var2,levels = rev(levels(dfsam2$Var2)))

ggplot(data = dfsam2,aes(x = Var2,y=Freq,fill = Var1)) +
  geom_bar(stat="identity",position = "fill",width = 0.5)+
  scale_y_continuous(labels = scales::percent)+
  scale_fill_manual(values = rev(c("#F8766D","#00BFC4")), name = "" )+
  theme_classic() + 
  labs(y = 'Fraction of different group',x="") + 
  coord_flip()  + 
  NoLegend()

saveRDS(scc_integrated,file="SCT_combine.sample.umap.rds")