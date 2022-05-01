#!/usr/bin/env Rscript
# -*- encoding: utf-8 -*-
#' """
#' @Manuscript  : Integrated Proteogenomic Characterization of Cholangiocarcinoma
#' @Time        : 2022/03/11 10:57:18
#' @Author      : RanPeng 
#' @Version     : R version 4.1.0
#' @Contact     : 2502388440@hotmail.com
#' @License     : (C)Copyright 2021-2022, DingLab-CHINA-SHANGHAI
#' @Description : Cohort description and Gene mutation of CCA
#' @Section     : Results - Multi-platform profiling of tumour samples
#' """

#=========================================================================
# Panel=sFig3 - Volcano plot of Overall survival for HRs (hazard ratio)
#=========================================================================

setwd("/Users/ranpeng/OneDrive - mail.ecust.edu.cn/项目文件/CCA-project/Figure/Fig3")
library(ggplot2)
library(ggrepel)
library(dplyr)

genes <- read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig3/protein_os/results/protein_filter_os_3_2021-3-22.csv", header = TRUE)
colnames(genes)
# genes$Significant <- ifelse(genes$CCA_pval < 0.05 & genes$CCA_HR > 1 , "P < 0.05", "Not Sig")
ggplot(genes, aes(x = log2(CCA_HR), y = -log10(CCA_pval))) +
  geom_point(aes(color = Significant)) +
  scale_size_continuous(range = c(1,4))+
  scale_color_manual(values = c("#EA686B", "#78A9CE","gray")) + # ("gray", "#7AA9CE")
  theme_bw(base_size = 12) + theme(legend.position = "right") +
  geom_text_repel(
    data = subset(genes, Significant == "good prognosis"|Significant == "bad prognosis"),
    aes(label = Gene),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle("Candidate biomarkers volcano")+ 
  ylim(2,8) +
  theme(axis.text.x=element_text(angle=0,hjust = 0.5,colour="black",family="Times",size=16), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Arial",size=16,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Arial",size = 20,face="plain"), #设置y轴标题的字体属性
        axis.title.x=element_text(family="Arial",size = 20,face="plain"), #设置y轴标题的字体属性
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.6), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="plain", family="Arial", colour="black",  #设置图例的子标题的字体属性
                                 size=16),
        legend.title=element_text(face="plain", family="Arial", colour="black", #设置图例的总标题的字体属性
                                  size=18),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())

## data filter
genes$Significant = ifelse(genes$CCA_pval < 0.05 & genes$CCA_HR > 1.2 & 
                             genes$iCC_pval < 0.05 & genes$iCC_HR > 1.1 & 
                             genes$eCC_pval < 0.05 & genes$eCC_HR > 1.1,"bad prognosis", 
                           ifelse(genes$CCA_pval < 0.05 & genes$CCA_HR < 0.8 & 
                                    genes$iCC_pval < 0.05 & genes$iCC_HR < 0.9 & 
                                    genes$eCC_pval < 0.05 & genes$eCC_HR < 0.9,"good prognosis","not sig"))


filter_os = genes%>% filter( Significant == "bad prognosis" | Significant == "good prognosis") 
write.csv(genes,'/Users/ranpeng/Desktop/CCA/Data/Fig3/protein_os/results/protein_filter_os_3_2021-3-22.csv')


data = read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig3/protein_os/results/217_tumor/protein_integra_os_change_pval.csv", header = TRUE)
data$eCC_pval = ifelse(data$Gene == data$Unnamed..0_y.1, data$log_rank_p, data$eCC_pval)
write.csv(data,"/Users/ranpeng/Desktop/CCA/Data/Fig3/protein_os/results/217_tumor/protein_integra_os_change_pval_1.csv")

#=========================================================================
#               Panel=sFig3 - dotplot biomarker
#=========================================================================

library(ggplot2)
library(tidyverse)
data = read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig3/protein_os/RNA-biomarker/biomaker.csv")
data$median = apply(data[2:length(data)],1,median)
data = data %>%filter(median > 0 ) %>% arrange(median)
list = data[,1]
data$Gene <- factor(data$Gene, levels = list, ordered=TRUE)

## ---------------  Melt data into column format ------------------
df <- gather(data, "variable", "value",2:length(data))

# Convert the variable dose from a numeric to a factor variable
ToothGrowth = df
ToothGrowth$Gene <- as.factor(ToothGrowth$Gene)
head(ToothGrowth)
# Basic dot plot
## 
p<-ggplot(ToothGrowth, aes(x=Gene, y=value)) + 
  geom_boxplot()+
  # ggtitle("Plot of length \n by dose") +
  xlab("Candidate biomarkers") + ylab("Log2-FC(Tumor/Normal)")+
  geom_hline(yintercept =-0.1,color = 'red', linetype="dashed") +
  ylim(-8,20)

p+ coord_flip() +
  theme(axis.text.x=element_text(angle=0,hjust = 0.5,colour="black",family="Times",size=16), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Arial",size=16,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Arial",size = 20,face="plain"), #设置y轴标题的字体属性
        axis.title.x=element_text(family="Arial",size = 20,face="plain"), #设置y轴标题的字体属性
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.6), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="plain", family="Arial", colour="black",  #设置图例的子标题的字体属性
                                 size=16),
        legend.title=element_text(face="plain", family="Arial", colour="black", #设置图例的总标题的字体属性
                                  size=18),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())

#=========================================================================
# Panel=sFig3 - PCA analysis for RNA, Proteome and phosphoproteomics 
#=========================================================================

library("factoextra")
library("FactoMineR")
library("ggthemes")
matrix = read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig3/RNA/RNA_TN_del50%.csv",
                  row.names = 1,header = 1)
matrix = t(matrix)
# group information
df = read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig3/RNA/RNA_ana.csv")
group = df$type
matrix = log2(matrix+1)
# matrix = t(scale(t(matrix)))
# matrix = as.matrix(matrix)

res.pca <- PCA(matrix,#[,-c(1)],
               graph = FALSE,scale.unit = FALSE)
# matrix$group <- c(rep('Normal',124),rep('Tumor',124))
# rep('NT',5),rep('TC',5),rep('T',5))#分组

fviz_pca_ind(res.pca,
             label = "none", # hide individual labels "none""all"
             repel=TRUE,
             col.ind = group,
             legend.title = "Group",
             #palette = "npg",
             #habillage = annotation$BAP1, # color by groups
             axes = c(1, 2),
             palette = c("#177cb0","#dc3023"),
             addEllipses = TRUE, # Concentration ellipses
             ellipse.level = 0.95) +
  theme_base(base_size = 12) +
  scale_shape_manual(values = c(16,17,19,19,19,19)) #自定义点的形状，分别为15， 19， 17。

#=========================================================================
# Panel=sFig3 - Volcano plot of differential expression proteins 
#=========================================================================

# setwd("/Users/ranpeng/OneDrive - mail.ecust.edu.cn/项目文件/CCA-project/Figure/Fig2")
library(ggplot2)
library(ggrepel)
genes <- read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig3/RNA/results/Tumor/diff_gene_results_pvalue.csv", header = TRUE)
colnames(genes)
genes$Significant <- ifelse(genes$adj.P.Val < 0.01 & genes$logFC > 1 , "tumor_up", 
                            ifelse(genes$adj.P.Val < 0.01 & genes$logFC < -1 , "Normal_up","Not Sig"))
ggplot(genes, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = Significant)) +
  scale_size_continuous(range = c(1,4))+
  scale_color_manual(values = c("#7AA9CE","gray","#EA686B")) + # ("gray", "#7AA9CE") c( "#EA686B","#7AA9CE","gray"))
  theme_bw(base_size = 12) + theme(legend.position = "right") +
  geom_text_repel(
    data = subset(genes, adj.P.Val<0.0001),
    aes(label = gene),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle("volcano of RNA")+
  theme(axis.text.x=element_text(angle=0,hjust = 0.6,colour="black",family="Arial",size=16), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.title.x=element_text(angle=0,hjust = 0.5,colour="black",family="Arial",size=16),
        axis.text.y=element_text(family="Arial",size=16,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Arial",size = 20,face="plain"), #设置y轴标题的字体属性
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.8), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="plain", family="Arial", colour="black",  #设置图例的子标题的字体属性
                                 size=12),
        legend.title=element_text(face="plain", family="Arial", colour="black", #设置图例的总标题的字体属性
                                  size=12),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+
  geom_hline(yintercept = 1.301,color = 'gray', linetype="dashed")+
  geom_vline(xintercept = -0.2,color = 'gray', linetype="dashed") +
  geom_vline(xintercept = 0.2,color = 'gray', linetype="dashed") 
# panel.border = element_rect(colour = "black", fill=NA, size=1)

#=========================================================================
# Panel=Fig3 - Transcription factor activity inference 
#=========================================================================

## install packagea
# if (!requireNamespace("BiocManager", quietly=TRUE))
#   + install.packages("BiocManager")
#  BiocManager::install("mixtools")
#  BiocManager::install("bcellViper")
#  BiocManager::install("viper")

library(dorothea)
library(bcellViper)
library(dplyr)
library(viper)

## -------------------------  load data  --------------------------------
## transform file 
# dataDirectory <- system.file("extdata", package="Biobase")
# dataDirectory <-"/Users/ranpeng/Desktop/CCA/CCA_script/tf_activities/dset.demo.csv"
exprsData <- read.csv("/Users/ranpeng/Desktop/CCA/Metadata/RNA/RNA_FPKM_codegene_deduplicate_del75%.csv",header=T,row.names = 1)
# exprsData = t(exprsData)
exprs <- as.matrix(exprsData)
dim(exprs)
# Create expression matrix
minimalSet <- ExpressionSet(assayData=exprs)
# Phenotypic data
pData <- read.csv("/Users/ranpeng/Desktop/CCA/Metadata/RNA/anatation.csv",header=T,row.names = 1)
dim(pData)
rownames(pData)
summary(pData)

index = intersect(colnames(exprs),rownames(pData))
pData = pData[index,]
all(rownames(pData)== colnames(exprs))
phenoData <- new("AnnotatedDataFrame", data=pData)

# Assembling an ExpressionSet
exampleSet <- ExpressionSet( assayData=exprs,
                             phenoData=phenoData)

# accessing expression data from bcellViper
data(bcellViper, package = "bcellViper")
load("/Users/ranpeng/Desktop/CCA/CCA_script/tf_activities/liver_regulon.Rdata")
regulon = data
# adjfile <- system.file("aracne", "bcellaracne.adj", package = "bcellViper")
# regul <- aracne2regulon(adjfile, dset, verbose = FALSE)
# print(regul)
dset = exampleSet
# 6.1 Generating the gene expression signatures
signature <- rowTtest(dset, "type", "Tumor", "Normal")
signature <- (qnorm(signature$p.value/2, lower.tail = FALSE) *
                + sign(signature$statistic))[, 1]
# 6.2 NULL model by sample permutations
nullmodel <- ttestNull(dset, "type", "Tumor", "Normal", per = 1000,
                       repos = TRUE, verbose = FALSE)

# 6.3 msVIPER
mrs <- msviper(signature, regulon, nullmodel, verbose = FALSE)
df=summary(mrs)
plot(mrs, cex = .7)
write.csv(df,"/Users/ranpeng/Desktop/CCA/Data/Fig3/TF_active/results/TF_inference75%_liver_mrs_ld.csv")

# 6.3.1 Leading-edge analysis
mrs_ld <- ledge(mrs)
df = summary(mrs_ld)

# 7.1 Bootstrap msVIPER
signature <- bootstrapTtest(dset, "type", "Tumor", "Normal", verbose = FALSE) ## "description", c("CB", "CC"), "N",
mrs <- msviper(signature, regulon, nullmodel, verbose = FALSE)
mrs <- bootstrapmsviper(mrs, "mode")
plot(mrs, cex = .7)

# 7.2 Shadow analysis
mrshadow <- shadow(mrs, regulators = 25, verbose = FALSE)
summary(mrshadow)

# 7.3 Synergy analysis
mrs <- msviperCombinatorial(mrs, regulators = 25, verbose = FALSE)
mrs <- msviperSynergy(mrs, verbose = FALSE)
summary(mrs)
plot(mrs, 25, cex = .7)

# 8 Virtual Inference of Protein-activity by Enriched Regulon analysis (VIPER)
vpres <- viper(dset, regulon, verbose = FALSE)
dim(vpres)
tmp <- rowTtest(vpres, "description", c("CB", "CC"), "N")
data.frame(Gene = rownames(tmp$p.value), t = round(tmp$statistic, 2),
           "p-value" = signif(tmp$p.value, 3))[order(tmp$p.value)[1:10], ]

# 8.1 Running VIPER with a null model
vpsig <- viperSignature(dset, "description", "N", verbose = FALSE)
vpres <- viper(vpsig, regulon, verbose = FALSE)
summary(vpres)

pos <- pData(vpres)[["description"]] %in% c("M", "CB", "CC")
d1 <- exprs(vpres)[, pos]
colnames(d1) <- pData(vpres)[["description"]][pos]
dd <- dist(t(d1), method = "euclidean")
heatmap(as.matrix(dd), Rowv = as.dendrogram(hclust(dd, method = "average")), symm = T)

dd <- viperSimilarity(d1)
heatmap(as.matrix(as.dist(dd)), Rowv = as.dendrogram(hclust(as.dist(dd),method = "average")), symm = T)

#=========================================================================
# Panel=Fig3 - Acessing (human) dorothea Transcription factor and Target data
#=========================================================================
library(dorothea)
# library(bcellViper)
library(dplyr)
library(viper)

## -----------------------"load data"------------------------------
# acessing (human) dorothea regulons
# for mouse regulons: data(dorothea_mm, package = "dorothea")
data(dorothea_hs, package = "dorothea")

## ------------------"run viper", message=FALSE---------------------
regulons = dorothea_hs %>%
  filter(confidence %in% c("A","B","C"))

## 运行实例数据
demo=read.csv("/Users/ranpeng/Desktop/CCA/Metadata/RNA/RNA_FPKM_codegene_deduplicate_del50%.csv",header=T,row.names = 1,stringsAsFactors=FALSE)
# demo = apply(demo, 2, as.numeric())
description<-demo$description
esptm<-as.matrix(demo)
tf_activities <- run_viper(  esptm, regulons, 
                             options =  list( method = "scale", 
                                              minsize = 5, 
                                              eset.filter = FALSE, cores = 1,
                                              verbose = FALSE))

write.csv(tf_activities,"/Users/ranpeng/Desktop/CCA/Data/Fig3/TF_active/TF_anatation_50%.csv")
# summary(tf_activities)
mrs <- msviper(esptm, regulon, nullmodel=NULL, minsize = 4, verbose = FALSE)
summary(mrs)
plot(mrs, cex = .7)

## TF active correlation with target gene



