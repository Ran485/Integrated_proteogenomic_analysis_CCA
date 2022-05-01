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
#       Panel=Fig2 - Arm events analysis
#=========================================================================
require("ggrepel")
require("ggplot2")
rm(list=ls(all=TRUE)) 

genes <- read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig2/Figure-2a/Forcal_arm_os(2020-12-27).csv", header = TRUE)
colnames(genes)
# genes$significant <- ifelse((genes$iCC_sig > 1.301| genes$eCC_sig > 1.301 ), "P_value < 0.05", "Not Sig")
genes$significant <-ifelse((genes$log_p > 1.301 & genes$log2HR > 0 ), 'red',
                           ifelse((genes$log_p > 1.301 & genes$log2HR < 0 ), 'royalblue', 'grey'))

p = ggplot(genes, aes(x = log2HR, y = log_p, color = log_p)) +
  geom_point(aes(color = log_p,size = Frequency)) +
  # scale_colour_gradient(low = "blue", high = "yellow") +
  # scale_color_manual(values = c("grey",'red',"royalblue", 'blue')) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  geom_text_repel( 
    data = subset(genes, log_p > 1.301 & abs(log2HR) > 0),
    aes(label = gene),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  ) + geom_hline(yintercept = 1.301,color = 'red', linetype="dashed")+
  geom_vline(xintercept = 0,color = 'red', linetype="dashed") +
  scale_fill_gradient(low = "blue", high = "yellow")

p + xlim(-4,4) + ylim(0,2.5)
# scale_x_reverse() + scale_y_reverse() 

#=========================================================================
#       Panel=Fig2B - Survival analysis
#=========================================================================
## Load the required packages
library(survival)
library(survminer)
library(tidyverse)
library(DT)
setwd('/Users/ranpeng/Desktop/2021-2-3')
rm(list=ls(all=TRUE))
### Load data
CCA = read.csv('/Users/ranpeng/Desktop/2021-3-7/iCC_eCC_compare/iCC_mutation_os_delna.csv')
# CCA = as.data.frame(t(CCA))
# eCC = CCA[CCA$Pathology == "eCCA",]
dim(CCA)
exprSet = CCA

# exprSet[exprSet== 0] = NA
# summary(exprSet)

## 对需要做生存分析的样本分组，把连续变量变成分类变量，这里选择的基因是BAP1,这里使用中位数
# group = ifelse(exprSet$BAP1 > median(exprSet$BAP1),'High','Low')
group = ifelse(exprSet$BAP1 > 0 ,'Mut','WT')
# ifelse(exprSet$average == 0 ,'WT','Amp'))
# exprSet<- CCA[!(CCA[,4]=='5'|CCA[,4]=='4'|CCA[,4]=='6'),]
sfit <- survfit(Surv(time, status) ~ group , data=exprSet)
ggsurvplot(sfit, conf.int=F, pval=TRUE,risk.table = TRUE, risk.table.col = "strata",
           palette = c("#FF9E29", "#86AA00", "#F94F21", "#916CA0", "#599BAD", "#DBD289", "#FB9A99", "#A4D873", "#99CAE0"))

#=========================================================================
#      Panel=Fig2 - mutiOmicsViz analysis
#=========================================================================

library("multiOmicsViz")
rm(list=ls(all=TRUE))
setwd("/Users/ranpeng/Desktop/CCA/Data/Fig2/multiomics/results")

sourceOmics <- read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig2/multiomics/multiomics_CNV(2021-4-8).csv",header=TRUE,row.names=1,
                        stringsAsFactors=FALSE,check.names=FALSE)
targetOmics1 <- read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig2/multiomics/multiomics_RNA.csv",header=TRUE,row.names=1,
                         stringsAsFactors=FALSE,check.names=FALSE)
targetOmics1 <- read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig2/multiomics/multiomics_protein.csv",header=TRUE,row.names=1,
                         stringsAsFactors=FALSE,check.names=FALSE)
targetOmics3 <- read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig2/multiomics/multiomics_phospho_1.csv",header=TRUE,row.names=1,
                         stringsAsFactors=FALSE,check.names=FALSE)
targetOmics1 = log2(targetOmics1+1)
## Calculate the muti-omics correlation
# targetOmics1 = log10(targetOmics1)
targetOmics1 = scale(targetOmics1)
# sig <- calculateCorForTwoMatrices(matrix1 = sourceOmics,matrix2 = targetOmics1,fdr=0.01)
# write.csv(sig,file = 'sig.csv',row.names = T)

targetOmicsList <- list()
targetOmicsList[[1]] <- targetOmics1
targetOmicsList[[2]] <- targetOmics2
targetOmicsList[[3]] <- targetOmics2
outputfile <- paste(tempdir(),"/heatmap",sep="")
## define function
custom.function = function(x) { multiOmicsViz(sourceOmics,"CNA","All",targetOmicsList,
                                              c("RNA","protein","phospho"),"All", x, 
                                              paste("CCA_CNA correlation with protein", x, sep = "_"), 
                                              nThreads = 3, legend=TRUE)
}

ilitern_num = c(0.05)
# --------------------- parallel compute ------------------------
# library(parallel)
# core.num <- detectCores(logical = F) #detect core numbers
# cl <- makeCluster(core.num - 1) # set core numbers in parellel is 16
# system.time(parSapply(cl, ilitern_num, custom.function(x)) )
# stopCluster(cl) ## stop cluster

custom.function(ilitern_num)
multiOmicsViz(sourceOmics,"CNA","All",targetOmicsList,
              "RNA","All", 0.05, "CCA_CNA correlation with RNA_0.05_scale", nThreads = NULL, legend=TRUE)

## Calculate signature correlation
matrix1 <- read.csv("/Users/ranpeng/Desktop/2020-11-28/mutiOmicsViz/multiomics_CNA_1.csv",header=TRUE,row.names=1,
                    stringsAsFactors=FALSE,check.names=FALSE)
matrix2 <- read.csv("/Users/ranpeng/Desktop/2020-11-28/mutiOmicsViz/multiomics_phospho_1.csv",header=TRUE,row.names=1,
                    stringsAsFactors=FALSE,check.names=FALSE)
sig <- calculateCorForTwoMatrices(matrix1=matrix1,matrix2=matrix2,fdr=0.01)


#=========================================================================
#      Panel=sFig - CNV-protein-cis volcano
#=========================================================================

setwd("/Users/ranpeng/OneDrive - mail.ecust.edu.cn/项目文件/CCA-project/Figure/Fig2")
library(ggplot2)
library(ggrepel)
genes <- read.csv("/Users/ranpeng/OneDrive - mail.ecust.edu.cn/项目文件/CCA-project/Data/Fig2/Figure-2a/CNV/6p21.2_del_cor.csv", header = TRUE)
colnames(genes)
genes$Significant <- ifelse(genes$cor_P < 0.05 & genes$cor_r > 0.2 , "P < 0.05", "Not Sig")
ggplot(genes, aes(x = cor_r, y = -log10(cor_P))) +
  geom_point(aes(color = Significant,size = Mutation_num)) +
  scale_size_continuous(range = c(1,4))+
  scale_color_manual(values = c("gray", "#EA686B")) + # ("gray", "#7AA9CE")
  theme_bw(base_size = 12) + theme(legend.position = "right") +
  geom_text_repel(
    data = subset(genes, cor_P < 0.05),
    aes(label = X6p21.2_del),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle("X6p21.2_del")

#=========================================================================
# Panel=Fig2d - arm events overall survival analysis and fisher exact test
#=========================================================================
data  = read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig2/Figure-2d/Arm_significant_os.csv")
colnames(data)

data$X14p_1 = ifelse(data$X14p> 0.3, "X14p_1_Amp", "X14p_1_WT")
# ifelse(data$X1q21.2< -0.3, "del", "WT"))
data$X16p_1 = ifelse(data$X16p > 0.3, "X16p_1_Amp","X16p_1_WT")
# ifelse(data$X6q14.3< -0.3, "del", "WT"))
data$X4p_1 = ifelse(data$X4p >0.3, "X4p_1_Amp","X4p_1_WT")
# ----------------- fisher exact test -------------------
test = table(data$X14p_1,data$X16p_1)
fisher.test(test)
# write.csv(data,"/Users/ranpeng/Desktop/CCA/Data/Fig2/Figure-2d/Arm_significant_os.csv")
test = table(data$X1q21.2_1,data$X6q14.3_1)
fisher.test(test)

# group = ifelse(data$X1q21.2_1 > 0.3, "Del", "WT")
# ifelse(data$X1q21.2< -0.3, "del", "WT"))
# data = na.omit(data)
# ---------- arm events overall survival analysis ------------
df= data
library(survival)
library(survminer)
sfit <- survfit(Surv(time, status) ~ group , data = df)
ggsurvplot(sfit, conf.int=F, pval=TRUE,risk.table = TRUE, risk.table.col = "strata",
           palette = c("#FF9E29", "#86AA00", "#F94F21", "#916CA0", "#599BAD", "#DBD289", "#FB9A99", "#A4D873", "#99CAE0"))

#=========================================================================
# Panel=Fig2 - focal arm events cis-effect on corresponding proteins and RNA
#=========================================================================

library(Hmisc)
library(dplyr)
## anatation
data1 = read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig2/CNV_Protein_RNA_cis/1.csv",header = T,row.names = 1)

## protein or phospho
data2 = read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig2/multiomics/multiomics_protein_delna30%.csv",header = T,row.names = 1)
data2 = t(data2)
index = intersect(rownames(data1),rownames(data2))
data2 = data2[index,]
# z_score<- function(data){
#   data = data
#   data = t(scale(t(data)))
#   # data = cbind(row_names,data);colnames(data)[1:2] = c("Gene_ID","phospho_ID")
#   data = delete.na(data)
#   ## 删除NA值
#   return(data)
# }
# data2 = scale(data2)
# sd(data2[,1])
## -------- 数据处理-磷酸化数据需要进行log转化 -------
data2[is.na(data2)] = 0
# data2<- log2(data2 + 1)

result<-c()
for(i in colnames(data2)){
  newdata = data.frame(cbind(data1[,1],data2[,i]))
  ## ————————  Correlation test  ——————————————
  res <- rcorr(as.matrix(newdata),type=c("spearman")) #pearson spearman
  cor_r = res$r[2]; p_val = res$P[2]; gene = i
  result.linshi<-cbind(gene,cor_r,p_val)
  result<-rbind(result,result.linshi)
}
result = data.frame(result)
# write.csv(result,"/Users/ranpeng/Desktop/CCA/Data/Fig2/Arm/6q14.3/SH3BGRL2_protein_cis/SH3BGRL2_protein_KEGG_cor.csv")
result$cor_r = as.numeric(result$cor_r);result$p_val = as.numeric(result$p_val)
result_filter <- filter(result, p_val < 0.05 & abs(cor_r) > 0.2)
(result_pos_cor <- filter(result, p_val < 0.05 & cor_r > 0.2))
(result_neg_cor <- filter(result, p_val < 0.05 & cor_r < -0.2))
# write.csv(result,"/Users/ranpeng/Desktop/CCA/Data/Fig6/arm_cor/6p22.2_cor.csv")

#=========================================================================
# Panel=Fig2 - boxplot for RNA and Protein expression
#=========================================================================

## protein
df = read.csv("/Users/ranpeng/Desktop/CCA/Data/Protein/protein_Tumor_217patients_Delna30%",row.names = 1)
df = data.frame(t(df))
ana = read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig3/anatation_boxplot.csv")
## RNA
df = read.csv("/Users/ranpeng/Desktop/CCA/Metadata/RNA/RNA_FPKM_codegene_deduplicate.csv",row.names = 1)
df = data.frame(t(df))
ana = read.csv("/Users/ranpeng/Desktop/CCA/Metadata/RNA/anatation.csv")
## multiomics
df = read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig2/RNA/1p36_33_arm_merge_protein.csv",row.names = 1)
df = data.frame(t(df))
# protein
ana = read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig7/RNA/ana_protein_binary.csv")
# RNA
ana = read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig7/RNA/arm_ana.csv")

## --------------- customed define function ----------------
Z_Score <- function(x){
  res <- (x - mean(x)) / sd(x)
  return(res)}

search_boxplot <- function(gene){
  # search = paste(gene,"$|^",gene,".y",sep = '')
  search = paste("^",gene,"$", sep="") ## 完全匹配
  index=grep(search, colnames(df))  
  # index=grep(gene, colnames(df)) 
  print(colnames(df)[index])
  result_data=df[,index];result_data = data.frame(result_data);result_data$id = rownames(df)
  
  ## merge group data
  ana$id = gsub("#",".",ana$id)
  result_data = merge(ana,result_data)
  result_data$result_data = as.numeric(result_data$result_data)
  ## filter low expression row
  # result_data = subset(result_data,result_data < 4 & result_data > -4)
  ## data transform 不需要则注释掉
  result_data$result_data = log10(result_data$result_data)
  ## z-score transform
  result_data$result_data = Z_Score(result_data$result_data)
  # result_data = subset(result_data,result_data < 0.5 & result_data > -0.5)
  ## boxplot
  library(ggpubr)
  ggboxplot(result_data, x = "type", y = "result_data",combine = F,
            # add = c("mean_se"),
            color = "type" ,
            fill = "type", alpha = 0.12,
            order =  c("Normal","Tumor"), #,"Cluster_3","Cluster_4","Cluster_5","Cluster_6"),
            palette = c("#86AA00","#F94F21","#916CA0","#599BAD","#DBD289"), ##palette = c("#F94F21","#EFC000","#00AEBA"),
            width = 0.5 ,x.text.angle = 45)+
    # rotate_x_text(angle=15,hjust=.5,vjust=1) + #FF9E29",
    stat_compare_means(label.y = c()) +                            # Global p-value
    # stat_compare_means(aes(group = Type),  label.y = c()) +
    labs(title= gene, x="", y = "Relative Expression")  +
    geom_jitter(alpha= 0.7,aes(colour = type),width = 0.15,height = 0) +
    theme_test()+
    theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))  # 控制坐标轴文本的角度
  
}
colnames(df)
# Specify the corresponding individual genes for visualization
search_boxplot("NADK")

#=========================================================================
# Panel=Fig2 - pathway enrichment result volcano plot
#=========================================================================

# setwd("/Users/ranpeng/OneDrive - mail.ecust.edu.cn/项目文件/CCA-project/Figure/Fig2")
library(ggplot2)
library(ggrepel)
genes <- read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig4/BAP1_mutation/data/BAP1_pro_cor_pathway_1.csv", header = TRUE)
colnames(genes)
# genes$Significant <- ifelse(genes$p < 0.05 & genes$HR > 1 , "Bad prognosis",
#                             ifelse(genes$p < 0.05 & genes$HR < 1 , "good prognosis","Not Sig"))
genes$Significant = ifelse(genes$Significant == "select" , "select","Not Sig")
ggplot(genes, aes(x = cor_r, y = -log10(p_val))) +
  geom_point(aes(color = Significant)) +
  scale_size_continuous(range = c(1,4))+
  scale_color_manual(values = c("#515151", "red")) + # ("gray", "#7AA9CE") c( "#EA686B","#7AA9CE","gray"))
  theme_bw(base_size = 12) + theme(legend.position = "right") +
  geom_text_repel(
    data = subset(genes, Significant == "select"),
    aes(label = gene),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle("significant pathway cor with AHNAK2 Mut")+
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
  geom_vline(xintercept = -0.15,color = 'gray', linetype="dashed") +
  geom_vline(xintercept = 0.15,color = 'gray', linetype="dashed") 
# panel.border = element_rect(colour = "black", fill=NA, size=1)


#=========================================================================
# Panel=Fig2 - CNV cis effect on protein / RNA
#=========================================================================

RNA = read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig2/CNV_Protein_RNA_cis/CNV_on_RNA_cor_cis.csv",row.names = 2)
protein = read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig2/CNV_Protein_RNA_cis/CNV_on_protein_cor_cis.csv",row.names = 2)
phosprotein = read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig2/CNV_Protein_RNA_cis/CNV_on_phosphoprotein_cor_cis.csv",row.names = 2)
## 取基因的交集
inner_gene = intersect(rownames(RNA),rownames(protein))
RNA = RNA[inner_gene,]
protein = protein[inner_gene,]
data = cbind(RNA,protein)

x <- read.delim(pipe("pbpaste"),header = T) #,row.names = 1) ## 读取剪切板的数据
gene = x$gene
df_sel = data[gene,]
write.csv(df_sel,"/Users/ranpeng/Desktop/CCA/Data/Fig2/CNV_Protein_RNA_cis/CNV_on_RNA_protein_inner_CAG_cor_cis.csv")

## ----------- draw figure1 arm barplot ---------------------
df = read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig2/Arm/arm_figure1.csv",row.names = 1)
df_count = df
df_count$Amp = rowSums(df > 0.1) ## 统计每行中大于0.1 的个数
df_count$Del = rowSums(df < -0.1)
df_count$Amp_Frequency = df_count$Amp/139
df_count$Del_Frequency = df_count$Del/139
data = df_count[,c("Amp_Frequency","Del_Frequency","Amp","Del")]
write.csv(data,"/Users/ranpeng/Desktop/CCA/Data/Fig2/Arm/Chromosome Arm/Arm_Number.csv")

#=========================================================================
# Panel=Fig2 - pathway moleculor volcano plot
#=========================================================================

# setwd("/Users/ranpeng/OneDrive - mail.ecust.edu.cn/项目文件/CCA-project/Figure/Fig2")
library(ggplot2)
library(ggrepel)
genes <- read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig7/Xcell_signature_cor/Neutrophils/V1/Neutrophils_cor_molecular.csv", header = TRUE)
colnames(genes)
# genes$Significant <- ifelse(genes$p < 0.05 & genes$HR > 1 , "Bad prognosis",
#                             ifelse(genes$p < 0.05 & genes$HR < 1 , "good prognosis","Not Sig"))
# genes$Significant = ifelse(genes$Significant == "select" , "select","Not Sig")
genes$size = ifelse(genes$Significant == "select" , 4, 1)
genes$pathway = as.character(genes$pathway)
options(ggrepel.max.overlaps = 30) ## 调整可以重叠的个数
# options(digits = 3)
ggplot(genes, aes(x = cor_r, y = -log10(p_val),color=pathway,size = size)) +
  geom_point() + ## aes(color = Significant)
  scale_size_continuous(range = c(1,2))+
  scale_color_brewer(palette="Paired") +
  # scale_color_manual(values = c("#515151", "red","#7AA9CE")) + # ("gray", "#7AA9CE") c( "#EA686B","#7AA9CE","gray"))
  theme_bw(base_size = 12) + theme(legend.position = "right") +
  geom_text_repel(
    data = subset(genes, Significant == "select"),
    aes(label = gene),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines",
    )
  )+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle("Significant pathway cor with PARP1 Protein")+
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
  geom_vline(xintercept = -0.15,color = 'gray', linetype="dashed") +
  geom_vline(xintercept = 0.15,color = 'gray', linetype="dashed") +
  scale_x_continuous(limits = c(-0.3, 0.6))+
  scale_y_continuous(breaks=seq(0, 15, 3))
# panel.border = element_rect(colour = "black", fill=NA, size=1)

## count Amp and Del number
df = read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig2/CNV_Protein_RNA_cis/CAG_24.csv",row.names = 1)
df$Amp_count = rowSums(df > 0.1)
df$Del_count = rowSums(df < -0.1)
df_count = df[,c("Amp_count","Del_count")]
write.csv(df_count,"/Users/ranpeng/Desktop/CCA/Data/Fig2/CNV_Protein_RNA_cis/CAG_24_count.csv")

#=========================================================================
# Panel=Fig2 - Prognosis Results of SH3BGRL2 in RNA level Forest
#=========================================================================
#install.packages("forestplot")
library(forestplot)
rs_forest <- read.csv('/Volumes/Samsung_T5/downloads/CCA/Data/Fig2/Arm/6q14.3/SH3BGRL3_RNA_os/SH3BGRL2_RNA_OS.csv',header = T)
# tiff('Figure 1.tiff',height = 6000,width = 7000,res= 600)
forestplot(labeltext = as.matrix(rs_forest[,1:6]),
           # 设置用于文本展示的列，此处我们用数据的前四列作为文本，在图中展示
           mean = rs_forest$HR,   #设置均值
           lower = rs_forest$HRL, #设置均值的lowlimits限
           upper = rs_forest$HRH, #设置均值的uplimits限
           is.summary = c(F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,T,F,F,F,F,F,F),
           #该参数接受一个逻辑向量，用于定义数据中的每一行是否是汇总值，若是，则在对应位置设置为TRUE，若否，则设置为FALSE；设置为TRUE的行则以粗体出现
           zero = 1, #设置参照值，此处我们展示的是HR值，故参照值是1，而不是0
           boxsize = 0.3, #设置点估计的方形大小
           lineheight = unit(17,'mm'),#设置图形中的行距
           colgap = unit(5,'mm'),#设置图形中的列间距
           lwd.zero = 3,#设置参考线的粗细
           lwd.ci = 2,#设置区间估计线的粗细
           col=fpColors(box='red', summary= "#8B008B",lines = 'black',zero = '#7AC5CD'),
           #使用fpColors()函数定义图形元素的颜色，从左至右分别对应点估计方形，汇总值，区间估计线，参考线
           xlab="The estimates",#设置x轴标签
           lwd.xaxis=2,#设置X轴线的粗细
           lty.ci = "solid",
           # xlog=TRUE,
           graph.pos = 4)#设置森林图的位置，此处设置为4，则出现在第四列
