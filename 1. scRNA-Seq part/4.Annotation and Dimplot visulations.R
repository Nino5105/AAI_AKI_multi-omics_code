# 4.Annotation and Dimplot visulations

# 4.1 library packages
library(Seurat)
library(ggplot2)
library(ggsci)

# 4.2 load scRNA datasets
scc_integrated <-readRDS("1.Data/4.SCT_combine.sample.anno.rds")
scc_integrated@active.ident <- scc_integrated$id
DimPlot(scc_integrated,label = T,label.size = 5)

new.cluster.ids <- c("T lymph/NK","Macro","T lymph/NK","PT","PT","PT","Macro","PT","T lymph/NK","PT","PT",
                     "ALH","Endo","DCT","PT","T lymph/NK","Peri","Macro","PT","Neutro","B lymph",
                     "DLH","Macro","ALH","Novel","T lymph/NK","Macro","Podo","Podo","CD-PC","Fibro", 
                     "DLH","Fibro","CD-IC","T lymph/NK","Podo","Endo","Fibro")

names(new.cluster.ids) <- levels(scc_integrated)
scc_integrated <- RenameIdents(scc_integrated, new.cluster.ids)
levels(scc_integrated) <- c("PT","DLH","ALH","DCT","CD-IC","CD-PC",
                            "Novel","Endo","Podo","Peri","Fibro",
                            "Macro","Neutro","B lymph","T lymph/NK")

table(scc_integrated@active.ident) 

# PT        DLH        ALH        DCT      CD-IC      CD-PC      Novel       Endo       Podo       Peri      Fibro 
# 17384       1193       2285       1365        274        491        626       1566       1167       1151       1056 
# Macro     Neutro    B lymph T lymph/NK 
# 8867        767        742      13277 

DimPlot(scc_integrated,label = T,label.size = 5)

saveRDS(scc_integrated,file="1.Data/4.SCT_combine.sample.anno.rds")

DimPlot(scc_integrated_new,label = T,label.size = 5,cols = c("PT" = '#E95C59',
                                          "DLH" = '#53A85F',
                                          "ALH" = '#F1BB72',
                                          "DCT" = '#F3B1A0',
                                          "CD-IC" = '#D6E7A3',
                                          "CD-PC" = '#57C3F3',
                                          "Novel" = '#E63863',
                                          "Endo" = '#E4C755',
                                          "Podo" = '#E59CC4',
                                          "Peri" = '#AB3282',
                                          "Fibro" = '#23452F',
                                          "Myeloid" = '#BD956A',
                                          "Neutro" = '#8C549C',
                                          "B lymph" = '#58A4C3',
                                          "T lymph/NK" = "#00BFC4"),group.by = "cell_type") + 
  theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title =  paste0("Cell annotation"))

# 4.3 Dimplot_Anno_each_type

# DimPlot(scc_integrated_new,label = T,label.size = 5,split.by = "type",cols = c("PT" = '#E95C59',
#                                                          "DLH" = '#53A85F',
#                                                          "ALH" = '#F1BB72',
#                                                          "DCT" = '#F3B1A0',
#                                                          "CD-IC" = '#D6E7A3',
#                                                          "CD-PC" = '#57C3F3',
#                                                          "Novel" = '#E63863',
#                                                          "Endo" = '#E4C755',
#                                                          "Podo" = '#E59CC4',
#                                                          "Peri" = '#AB3282',
#                                                          "Fibro" = '#23452F',
#                                                          "Myeloid" = '#BD956A',
#                                                          "Neutro" = '#8C549C',
#                                                          "B lymph" = '#58A4C3',
#                                                          "T lymph/NK" = "#00BFC4"),group.by = "cell_type") + 
#   theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title =  paste0("Cell annotation"))

# 4.4 Dimplot_Anno_each_sample

DimPlot(scc_integrated,label = T,label.size = 5,ncol = 3,split.by = "sample",cols = c("PT" = '#E95C59',
                                                                           "DLH" = '#53A85F',
                                                                           "ALH" = '#F1BB72',
                                                                           "DCT" = '#F3B1A0',
                                                                           "CD-IC" = '#D6E7A3',
                                                                           "CD-PC" = '#57C3F3',
                                                                           "Novel" = '#E63863',
                                                                           "Endo" = '#E4C755',
                                                                           "Podo" = '#E59CC4',
                                                                           "Peri" = '#AB3282',
                                                                           "Fibro" = '#23452F',
                                                                           "Macro" = '#BD956A',
                                                                           "Neutro" = '#8C549C',
                                                                           "B lymph" = '#58A4C3',
                                                                           "T lymph/NK" = "#00BFC4")) + theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title =  paste0("Cell annotation"))

# 4.5 Count the proportions of each group and each sample in different cell types

cellnum <- as.data.frame(table(Idents(scc_integrated),scc_integrated$orig.ident))
cellnum

table(Idents(scc_integrated),scc_integrated$orig.ident)
#             KA1  KA2  KA3  KC1  KC2  KC3
# PT          913 1058  875 5558 3537 5443
# DLH         255  343  266  118  106  105
# ALH         326  338  425  227  574  395  #-8.9%
# DCT         151  220  202  117  353  322  #-13.5%
# CD-IC        41   38   40   10   95   50  #-2.3%
# CD-PC        46   67   42   75  138  123  #-53.8%
# Novel       210  159  117   50   23   67
# Endo        298  299  457  116  193  203
# Podo        158  185  242  179  182  221 
# Peri        430  292  202   95   22  110
# Fibro       392  321  173   65   34   71
# Macro      2431 2048 2417  700  353  918
# Neutro      277  153  116   95   60   66
# B lymph      83   93   58  157   97  254
# T lymph/NK 4057 3712 3929  719  505  355

prop.table(table(Idents(scc_integrated),scc_integrated$orig.ident),margin = 2)

# KA1         KA2         KA3         KC1         KC2         KC3
# PT         0.090683353 0.113446279 0.091517624 0.671174979 0.563934949 0.625416523
# DLH        0.025327771 0.036778898 0.027821358 0.014249487 0.016900510 0.012064805
# ALH        0.032379817 0.036242762 0.044451417 0.027412148 0.091517857 0.045386648
# DCT        0.014998014 0.023589964 0.021127497 0.014128728 0.056281888 0.036998736
# CD-IC      0.004072308 0.004074630 0.004183663 0.001207584 0.015146684 0.005745145
# CD-PC      0.004568931 0.007184216 0.004392846 0.009056877 0.022002551 0.014133058
# Novel      0.020858164 0.017049110 0.012237214 0.006037918 0.003667092 0.007698495
# Endo       0.029598729 0.032060905 0.047798347 0.014007970 0.030771684 0.023325290
# Podo       0.015693286 0.019837015 0.025311160 0.021615747 0.029017857 0.025393542
# Peri       0.042709575 0.031310315 0.021127497 0.011472044 0.003507653 0.012639320
# Fibro      0.038935240 0.034419901 0.018094342 0.007849294 0.005420918 0.008158106
# Macro      0.241458085 0.219601115 0.252797824 0.084530854 0.056281888 0.105480869
# Neutro     0.027512912 0.016405747 0.012132622 0.011472044 0.009566327 0.007583592
# B lymph    0.008243941 0.009972121 0.006066311 0.018959063 0.015465561 0.029185338
# T lymph/NK 0.402959873 0.398027021 0.410940278 0.086825263 0.080516582 0.040790532


# scc_integrated <- readRDS("SCT_combine.sample.anno.rds")
DimPlot(scc_integrated)
table(scc_integrated@active.ident)
scc_integrated$cell_type <- scc_integrated@active.ident
table(scc_integrated$cell_type)

levels <- c("PT","DLH","ALH","DCT","CD-IC","CD-PC","Novel","Endo","Podo","Peri","Fibro","Macro","Neutro","B lymph","T lymph/NK")

dfsam1 <- as.data.frame(table(scc_integrated$type,scc_integrated$cell_type))
dfsam1$Var1 <- factor(dfsam1$Var1,levels = rev(c("Control","Treatment")))
dfsam1$Var2 <- factor(dfsam1$Var2,levels = rev(levels))

p1 <- ggplot(data = dfsam1,aes(x = Var2,y=Freq,fill = Var1)) +
  geom_bar(stat="identity",position = "fill",width = 0.7)+
  scale_y_continuous(labels = scales::percent)+
  scale_fill_manual(values = rev(c("#00BFC4","#F8766D")), name = "" )+
  theme_classic() + 
  labs(y = 'Fraction of different type',x="") + 
  coord_flip()  + 
  NoLegend()
p1

dfsam2 <- as.data.frame(table(scc_integrated$sample,scc_integrated$cell_type))
dfsam2$Var1 <- factor(dfsam2$Var1,levels = rev(c("KC1","KC2","KC3","KA1","KA2","KA3")))
dfsam2$Var2 <- factor(dfsam2$Var2,levels = rev(levels))

p2 <- ggplot(data = dfsam2,aes(x = Var2,y=Freq,fill = Var1)) +
  geom_bar(stat="identity",position = "fill",width = 0.7)+
  scale_y_continuous(labels = scales::percent)+
  theme_classic() + 
  labs(y = 'Fraction of different sample',x="") +
  coord_flip() + 
  scale_fill_nejm() + 
  NoLegend()

p1 /p2 #3x8


saveRDS(scc_integrated,file="1.Data/4.SCT_combine.sample.anno.rds")


# 4.6 renew the metadata information

All_meta_data <- read.csv("1.Data/ALL_meta_data.csv",row.names = 1,check.names = F)
scc_integrated_new <- scc_integrated
scc_integrated_new@meta.data <- All_meta_data
Idents(scc_integrated_new) <- scc_integrated_new$cell_type
DimPlot(scc_integrated_new,group.by = "cell_type_sub",label = T)
table(scc_integrated_new$cell_type)
# ALH    B lymph      CD-IC      CD-PC        DCT        DLH       Endo      Fibro 
# 2416        742        272        481       1239       1201       1414        722 
# Myeloid     Neutro      Novel    Peri        Podo      PT     T lymph/NK 
# 8867        767        626       1637       1166      17384      13277 

saveRDS(scc_integrated_new,file="1.Data/4.SCT_combine.sample.anno_sub.rds")
