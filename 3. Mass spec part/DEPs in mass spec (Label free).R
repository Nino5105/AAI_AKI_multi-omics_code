# 1.载入R包

library(edgeR) 
library(limma)
library(ggplot2)
library(ggrepel)
library(DMwR2)
library(mice)
library(clusterProfiler)
library(org.Mm.eg.db)
library(patchwork)

#	2.数据预处理

raw_data <- read.csv("1.raw_data/1.raw_data.csv",header = T,row.names = 1,check.names = F)
head(raw_data)
md.pattern(raw_data) # 查看 NA 值的分布情况
colnames(raw_data)
trim_data <- raw_data[,c(2,6:8,3:5)]
head(trim_data)

# DMSO <- trim_data[,c(1:4)]
# DMSO <- knnImputation(DMSO)  # 采用 KNN 方法对数据缺失值进行补齐
# head(DMSO)
# 
# celp <- trim_data[,c(5:8)]
# celp <- knnImputation(celp)  # 采用 KNN 方法对数据缺失值进行补齐
# head(celp)
# 
# trim_data <- cbind(DMSO,celp)
# md.pattern(trim_data)

list <- which(duplicated(trim_data$Gene))
trim_data = trim_data[-list,]
rownames(trim_data) <- trim_data$Gene

trim_data <- trim_data[,-1]
head(trim_data)

write.csv(trim_data,"Protein expression in Label free.csv")

#	3.差异分析

data = log2(trim_data+1)
condition = factor(c(rep("Control",3),rep("Treatment",3)),levels = c("Control","Treatment"))    

genelist = DGEList(counts = data, group = condition)  
design <- model.matrix(~condition)                            
colnames(design) <- levels(condition)  
rownames(design) <- colnames(data)
design

fit <- lmFit(data, design)
fit <- eBayes(fit, trend=TRUE)
DEP_result <- topTable(fit,coef=2,n=Inf)     
head(DEP_result)

cut_off_pvalue = 0.05
cut_off_logFC = 0.263

DEP_result$change = ifelse(DEP_result$adj.P.Val < cut_off_pvalue & abs(DEP_result$logFC) >= cut_off_logFC, 
                           ifelse(DEP_result$logFC> cut_off_logFC ,'Up(1903)','Down(667)'),'Stable(2247)')
# DEP_result$label = ifelse(rownames(DEP_result) %in% c("Pkm","Hmgb1","Ldha"), rownames(DEP_result),'')
table(DEP_result$change)

# Down(667) Stable(2247)     Up(1903) 
# 667         2247         1903   


DEP_up <- row.names(DEP_result[DEP_result$change == "Up(1903)",]) # n = 1903
head(DEP_up)

DEP_down <- row.names(DEP_result[DEP_result$change == "Down(667)",]) # n = 667
head(DEP_down)


write.csv(DEP_result,"Different expressed protein between Treatment vs control in AA-Kideny.csv")


# options(ggrepel.max.overlaps = Inf)
# 
# pdf("Clep vs DMSO volcano plot.pdf",width = 6,height = 5)

ggplot(data=DEP_result, aes(x=logFC, y =-log10(adj.P.Val),colour=change)) +
  geom_point(alpha=0.8, size=1)+
  scale_color_manual(values=c("#546de5","#000000","#ff4757"))+
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.8,alpha=0.8)+
  geom_vline(xintercept = c(0.263,-0.263),lty=4,lwd=0.8,alpha=0.8)+
  labs(x="log2 (fold change) [Treatment vs Control]",y="-log10 (FDR)")+
  ggtitle("Label free proteomics \n Signigicant protein (1903 Up, 667 down)")+
  theme_bw() + 
  xlim(-5,5) + ylim(0, 7) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) 
# + geom_label_repel(data = DEP_result, aes(x=logFC, y =-log10(adj.P.Val),label = label),
#                    size = 4,box.padding = unit(2, "lines"),
#                    point.padding = unit(0.8, "lines"),
#                    segment.color = "black",
#                    show.legend = FALSE,max.overlaps = Inf)
# dev.off()

# 4.绘制热图

# 4.1 all


pheatmap::pheatmap(data,treeheight_row=0, treeheight_col=0,show_rownames = F,cluster_cols = F,
                   cellwidth = 15, cellheight = 1,main = "Protein expression ratio of Celp/DMSO_mean 
                   vs Celp+C/DMSO_mean(all protein)",
                   angle_col = 45,scale = "row", border_color = "black")


# 4.2 up

pheatmap::pheatmap(data[row.names(DEP_result[DEP_result$change == "Up(262)",]),],
                   treeheight_row=0, treeheight_col=0,show_rownames = F,cluster_cols = F,
                   cellwidth = 15, cellheight = 1,main = "Protein expression ratio of Celp/DMSO_mean 
                   vs Celp+C/DMSO_mean(up regulate protein)",
                   angle_col = 45,scale = "row", border_color = "black")



# 5.GO富集分析
# （4）分别对上下调基因进行GO富集分析

# ①进行GO分析
Go_BP_up <- enrichGO(gene = DEP_up, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                     ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)

Go_BP_down <- enrichGO(gene = DEP_down, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                       ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)


# ②进行可视化
p1 <- dotplot(Go_BP_up, showCategory = 10,title = "The GO enrichment(BP) analysis of Treatment vs Control up DEP in Label free")
p2 <- dotplot(Go_BP_down, showCategory = 10,title = "The GO enrichment(BP) analysis of Treatment vs Control down DEP in Label free")

pdf("GO enrichment analysis of Treatment vs Control DEG in Label free.pdf",width = 12,height = 10)

p1/p2

dev.off()

# ③导出富集结果
write.csv(Go_BP_up@result,"The GO enrichment(BP) analysis of Treatment vs Control up DEG in scc_integrated.csv" )
write.csv(Go_BP_down@result,"The GO enrichment(BP) analysis of Treatment vs Control down DEG in scc_integrated.csv" )


