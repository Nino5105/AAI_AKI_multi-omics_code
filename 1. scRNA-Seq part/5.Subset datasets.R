# (0) 加载R包
library(dplyr)
library(Seurat)
library(patchwork)
library(ggsci)
library(ggplot2) 
options(future.globals.maxSize = 20000 * 1024^2)

#（1）读取数据
scc_integrated <-readRDS("1.Data/4.SCT_combine.sample.anno.rds")
table(scc_integrated@active.ident)

table(ALL@active.ident)
# PT        DLH        ALH        DCT      CD-IC      CD-PC      Novel       Endo       Podo       Peri 
# 17384       1193       2285       1365        274        491        626       1566       1167       1151 
# Fibro      Macro     Neutro    B lymph  T lymph/NK 
# 1056       8867        767        742      13277 

#（2）提取各类细胞
PT <- subset(scc_integrated,idents = "PT")
saveRDS(PT,file="1.Data/5.1 PT.rds")

Other_Epi <- subset(scc_integrated,idents = c("DLH","ALH","DCT","CD-IC","CD-PC","Podo"))
saveRDS(Other_Epi,file="1.Data/5.2 Other_Epi.rds")

Stroma <- subset(scc_integrated,idents = c("Endo","Peri","Fibro"))
saveRDS(Stroma,file="1.Data/5.3 Stroma.rds")

T_NK_lymph <- subset(scc_integrated,idents =c("T lymph/NK"))
table(T_NK_lymph@active.ident)
saveRDS(T_NK_lymph,file="1.Data/5.4 SCT_T_NK_lymph.rds")

Other_immune <- subset(scc_integrated,idents = c("Macro","B lymph","Neutro"))
saveRDS(Other_immune,file="1.Data/5.5 Other_immune.rds")




