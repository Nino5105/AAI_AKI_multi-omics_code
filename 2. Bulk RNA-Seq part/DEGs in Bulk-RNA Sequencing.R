library(limma)
library(edgeR)
library(ggplot2)
library(patchwork)
library(clusterProfiler)
library(org.Mm.eg.db)

raw_data <- read.table("AA_Kidney_Bulk_RNA.txt",check.names = F,header = T)
head(raw_data)
trim_data <- raw_data[,c(1,6:12)]
colnames(trim_data) <- c("Geneid","Length","KA1","KA2","KA3","KC1","KC2","KC3")
head(trim_data)

gene.df <- bitr(trim_data$Geneid, fromType ="ENSEMBL",
                toType = "SYMBOL",
                OrgDb = org.Mm.eg.db) # 对基因id进行转换（ENSEMBL → SYMBOL）

merge_data <- data.frame(gene.df,trim_data[match(gene.df$ENSEMBL,trim_data$Geneid),])

table(!duplicated(merge_data$SYMBOL))
dulpicated_list <- which(duplicated(merge_data$SYMBOL))
merge_data_rd <- merge_data[-dulpicated_list,]
row.names(merge_data_rd) <- merge_data_rd$SYMBOL

counts <- merge_data_rd[,c(8:10,5:7)]
head(counts)

write.csv(counts,"Bulk RNA counts.csv")

condition = factor(c(rep("Control",3),rep("Treatment",3)),levels = c("Control","Treatment"))
condition

genelist = DGEList(counts = counts, group = condition)
design <- model.matrix(~condition)
colnames(design) <- levels(condition)
rownames(design) <- colnames(counts)
design

keep <- rowSums(cpm(genelist)>1)>=2 
genelist.filted <- genelist[keep,,keep.lib.sizes=FALSE] 
genelist.filted

#         group lib.size norm.factors
# KC1   Control 27615473            1
# KC2   Control 32488202            1
# KC3   Control 35034073            1
# KA1 Treatment 37398755            1
# KA2 Treatment 43961408            1
# KA3 Treatment 32708642            1

genelist.norm <- calcNormFactors(genelist.filted)  
logCPM <- cpm(genelist.norm, log=TRUE, prior.count=3)

fit <- lmFit(logCPM, design)
fit <- eBayes(fit, trend=TRUE)
DEP_result <- topTable(fit,coef=2,n=Inf)
DEP_result
DEP_result$type = ifelse(DEP_result$adj.P.Val < 0.05 & abs(DEP_result$logFC) >= 1,
                         ifelse(DEP_result$logFC > 1 ,'Up(2871)','Down(1794)'),'Stable(9936)')
table(DEP_result$type)
# Down(1794) Stable(9936)     Up(2871) 
#   1794         9936         2871 

DEG_up <- row.names(DEP_result[DEP_result$type == "Up(2871)",]) # n = 2871
head(DEG_up)

DEG_down <- row.names(DEP_result[DEP_result$type == "Down(1794)",]) # n = 1794
head(DEG_down)


write.csv(DEP_result,"DEG in Bulk-RNA Seq.csv" )


ggplot(data=DEP_result, aes(x=logFC, y =-log10(adj.P.Val),colour=type)) +
  geom_point(alpha=0.6, size=1)+
  scale_color_manual(values=c("#546de5","#000000","#ff4757"))+
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.8,alpha=0.8)+
  geom_vline(xintercept = c(1,-1),lty=4,lwd=0.8,alpha=0.8)+
  labs(x="log2 (fold change)",y="-log10 (FDR)")+
  ggtitle("Bulk RNA-seq \n Signifianct genes")+ 
  theme_bw() + xlim(-10,10) + ylim(0, 8) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())  # 6x5

Go_BP_up <- enrichGO(gene = DEG_up, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                     ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)

Go_BP_down <- enrichGO(gene = DEG_down, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                       ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)


p1 <- dotplot(Go_BP_up, showCategory = 10,title = "The GO enrichment(BP) analysis of Treatment vs Control up DEG in Bulk-RNA Seq")
p2 <- dotplot(Go_BP_down, showCategory = 10,title = "The GO enrichment(BP) analysis of Treatment vs Control down DEG in Bulk-RNA Seq")

pdf("GO enrichment analysis of Treatment vs Control DEG in Bulk-RNA Seq.pdf",width = 10,height = 10)

p1/p2

dev.off()

write.csv(Go_BP_up@result,"The GO enrichment(BP) analysis of AA vs Con up DEG in scc_integrated.csv" )
write.csv(Go_BP_down@result,"The GO enrichment(BP) analysis of TAA vs Con down DEG in scc_integrated.csv" )
