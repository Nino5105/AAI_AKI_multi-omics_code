library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Seurat)

data<-readRDS("data.rds")
meta<-data@meta.data

color<-read.table("color.use",row.names=1,sep="\t",check.names=F)
color$V2<-paste("#",color$V2,sep="")
color<-color[sort(row.names(color)),,drop=FALSE]

Control_meta<-subset(meta,type=="Control")
Treatment_meta<-subset(meta,type=="Treatment")

Control_data<-as.matrix(data@assays$SCT@data[,as.character(row.names(Control_meta))])
Treatment_data<-as.matrix(data@assays$SCT@data[,as.character(row.names(Treatment_meta))])

#########Control#########
cellchat <- createCellChat(object = Control_data, meta = Control_meta, group.by = "cell_type")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

#Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

#Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
pdf("Control-CellChat.pdf")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions"",color.use=color$V2)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength"",color.use=color$V2)
dev.off()
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
saveRDS(cellchat, "Control_cellchat.RDS")

########################merge##############################
control<-readRDS("Control_cellchat.RDS")
treatment<-readRDS("Treatment_cellchat.RDS")
object.list <- list(Control=control,Treatment=treatment)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
pdf("Diff_Num.pdf")
compareInteractions(cellchat, show.legend = F, group = c(1,2))
dev.off()
pdf("Diff_Heatmap.pdf")
netVisual_heatmap(cellchat,color.use=color$V2)
dev.off()
pdf("Diff_circ.pdf")
netVisual_diffInteraction(cellchat, weight.scale = T,color.use=color$V2)
dev.off()

#####select######
c1<-subsetCellChat(cellchat,idents.use=c("CD4+ Te","CD4+ Tem","CD4+ Tn","CD4+ Treg","CD8+ CTL","CD8+ Tem","CD8+ Tn","Endo","Fibro","GE","Macro M1","Macro M2","Macro Pro","Mast cell","Monocytes","MyoFibro","NK","PCT","PST","PT-S","T pro"))
groupSize <- as.numeric(table(c1@idents))
netVisual_circle(c1@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_bubble(c1, remove.isolate = FALSE,sources.use=c("PCT","PST","PT-S"),targets.use=c("Macro M1","Macro M2","CD8+ CTL","CD8+ Tem","GE","Fibro","MyoFibro"))+coord_flip()

#pairLR.use <- extractEnrichedLR(cellchat, signaling =cellchat@netP$Control$pathways)

#selected ligand-receptor pair
LR_list <- c("H2-D1_CD8A",
           "H2-D1_CD8B1",
           "H2-K1_CD8A",
           "H2-K1_CD8B1",
           "MIF_CD74_CXCR4",
           "MIF_CD74_CD44",
           "CCL5_CCR1",
           "CCL5_CCR5",
           "PTPRC_MRC1",
           "SPP1_CD44",
           "SPP1_ITGA4_ITGB1",
           "SPP1_ITGA9_ITGB1",
           "SPP1_ITGA8_ITGB1",
           "CCL6_CCR2",
           "COL1A2_ITGA1_ITGB1",
           "COL1A2_ITGA9_ITGB1",
           "SELL_CD34",
           "SELL_PODXL")
pair<-data.frame(LR_list)
colnames(pair)<-"interaction_name"

#selected cell-cell
cell<-c("PT-S","PCT","PST","MyoFibro","Macro M1","Macro M2","GE","Fibro","CD8+ CTL","CD8+ Tem")
out<-netVisual_bubble(cellchat,remove.isolate = FALSE,sources.use=cell,targets.use=cell,pairLR.use=pair,comparison=c(2,1),angle.x=45,color.text=c("blue","red"),return.data=TRUE)
source.target<-levels(out$communication$source.target)
col<-c(141,142,#PCT-CTL
161,162,#PST-CTL
181,182,#PT-S-CTL
143,144,#PCT-CD8+Tem
163,164,#PST-CD8+Tem
183,184
)
select<-(source.target[col])


source("netVisual_bubble2.r")#select cell-pair
netVisual_bubble2(cellchat,remove.isolate = FALSE,sources.use=cell,targets.use=cell,pairLR.use=pair,comparison=c(2,1),angle.x=45,color.text=c("blue","red"),select.source.target=select)+theme(axis.text.x=element_text(size=6))
