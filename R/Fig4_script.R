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
#       Panel=Fig4 - Heatmap of proteome subtypes
#=========================================================================

library(pheatmap)
library(RColorBrewer)
# rm(list = ls())
matrix=read.csv("/Volumes/Samsung_T5/downloads/CCA/Data/Fig4/NMF_heatmap.csv",header=T,skip = 0,row.names = 1)
#chending_color<-colorRampPalette(c("#3090D1","#FFFFFF","red"))(50)
chending_color<-colorRampPalette(c("#0690F7","#FFFFFF","red"))(50)
#chending_color <- colorRampPalette(brewer.pal(8, "PiYG"))(50)
# Generate annotations for rows and columns
annotation_col = read.csv ('/Users/ranpeng/Desktop/data/Fig4//Fig4A_clinical_feature_anatation_20220326.csv',header = T,row.names = 1)
# annotation_row = read.csv ('/Users/ranpeng/Desktop/胆管癌/data/2020-8-21/annotation_row_all.csv',header = T,row.names = 1)
ann_colors <- list(
  Age = c("#FFFFCB", "#C1E598", "#8BF28B", "#3AC160", "#006837"),
  Gender = c(Female = "#FCBA83", Male = "#E24F3B"),
  Pathology = c(iCCA = "#FFC51A", eCCA = "#a767ff"),
  iCCA_small_moleculor = c(smallBD = "#89BCDF", largeBD = "#2080C3", NA_1 = "#ededed", NA_2 = "#ededed"),
  CCP_Class = c(cluster_1 = "#1F78B4", cluster_2 = "#33A02C", cluster_3 = "#92CFEA", NA_1 = "#ededed"),
  Proteome_Subtype = c(S_I = "#FF800D", S_II = "#0E65AD", S_III = "#FF0027"),
  Grade = c(GX = "#ffe0e9", G1 = "#ffc2d4", G1_2 = "#ff9ebb", G2 = "#ff7aa2", G2_3 = "#e05780", G3 = "#b9375e"),
  TNM_Stage = c(I = "#0272BC", II = "#ABD373", III = "#FCED68", IV = "#F8941D"),
  iccStage = c(IA = "#caf0f8", IB = "#90e0ef", II = "#48cae4", IIIA = "#0096c7", IIIB = "#0077b6", IV = "#023e8a", NA_1 = "#ededed"),
  eccStage = c(I = "#caf0f8", II = "#90e0ef", IIIA = "#48cae4", IIIB = "#0096c7", IIIC = "#0077b6", IVB = "#023e8a", NA_1 = "#ededed"),
  Lymphovascular_invasion = c(No = "#fca311", Yes = "#14213d", Unknown = "#EDEDED"),
  Perineural_invasion = c(No = "#fbd87f", Yes = "#c62c43", Unknown = "#EDEDED"),
  Presence_of_fluke_infection = c(No = "#EDEDED", Yes = "#000000"),
  Type_2_diabetes = c(No = "#EDEDED", Yes = "#14213d"),
  Choledocholithiasis = c(No = "#020202", Yes = "#72D672", Unknown = "#EDEDED"),
  Cholecystolithiasis = c(No = "#020202", Yes = "#72D672", Unknown = "#EDEDED"),
  Cholelithiasis = c(No = "#EDEDED", Yes = "#f5cb5c"),
  Hypertension = c(No = "#EDEDED", Yes = "#086375"),
  CA199 = c(Normal = "#ffa5ab", Elevated = "#bc4877", NA_1 = "#ededed"),
  Tbil = c(Normal = "#1282a2", Elevated = "#034078", NA_1 = "#ededed"),
  ALT = c(Below = "#eff48e", Normal = "#d2e603", Elevated = "#3e978b", NA_1 = "#ededed"),
  GGT = c(Normal = "#68adff", Elevated = "#2b69c9", NA_1 = "#ededed"),
  HBV = c(Ever = "#376770", Never = "#9ac9d2", NA_1 = "#ededed"),
  PT = c(Below = "#e1f48f", Normal = "#f9d073", Elevated = "#e5643c", NA_1 = "#ededed")
)
matrix = as.matrix(matrix)
pheatmap(matrix,cluster_rows=0,cluster_cols=0,
         clustering_distance_cols = "correlation",fill = T,
         clustering_distance_rows = "correlation",border_color ="white", na_col = "white",
         col=chending_color,show_rownames=T,show_colnames=F,display_numbers=F,
         width = 4.85,height = 5,fontsize_number=12,number_color="black",number_format="%.1f",
         annotation_col = annotation_col , annotation_colors = ann_colors) #, annotation_row = annotation_row)


#=========================================================================
#  Pathology distribution across the proteome subtypes
#=========================================================================
source("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/R/stacked_barplot.R")

## ------------------------ load Pathology data -----------------------------
data = read.csv('/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/Fig1/clinical_merge20220406.csv')
colnames(data)
## ----------------- 指定 category_var_colname 进行分析 --------------------
stacked_barplot(data, target = 'Proteome_Subtype',category_var_colname = 'Pathology')


#=========================================================================
#  Panel=Fig4 - Analysis of differentially expressed proteins
#=========================================================================

library(limma)
library(ComplexHeatmap)
library(Hmisc)
library(circlize)
library(BBmisc)

options(stringsAsFactors = F)
warnings('off')
rm(list = ls())
## --------- multi_group_DEP function ----------------------------
Z_Score <- function(x){
  res <- (x - mean(x)) / sd(x)
  return(res)}

Get_subtype_info <- function(meta){
  f <- factor(meta$type);table(f)
  subtype_info = data.frame(table(f))
  rownames(subtype_info) = subtype_info[,1]
  subtype_list = rownames(subtype_info) ## 获取分组信息
  return(subtype_list)
}

## -------------- load file ----------------------
if(TRUE){
  input_path <- '/Users/ranpeng/Desktop/CCA/Data/Fig4/data_wash/Fig_4_merge_CNV_20210325.csv'
  anno_path <- '/Users/ranpeng/Desktop/CCA/Data/Fig4/anatation.csv'
  input_data <- read.csv(input_path, header=T, row.names=1, stringsAsFactors = F,)
  meta <- read.csv(anno_path)
  meta$id = gsub("#",".",meta$id)
  exprSet <- as.matrix(input_data)
  
  ## 如果存在空值进行下面的处理
  # tm1<-testmatrix[,-which(apply(testmatrix,2,function(x)all(is.na(x))))]  ## 删除所有是空值的列
  # tm2<-tm1[-which(apply(testmatrix,1,function(x)all(is.na(x)))),] ## 删除所有是空值的行
  exprSet<-exprSet[,-which(apply(exprSet,2,function(x) all(is.na(x))))]
  meta <- meta[which(meta$id %in% colnames(exprSet)),]
  ## 数据处理,磷酸化数据需要进行log转化
  exprSet[is.na(exprSet)] = 0
  # exprSet<- log2(exprSet + 1)
  ## 数据存储路径
  setwd('/Users/ranpeng/Desktop/CCA/Data/Fig4/gene_mutation_results')
  outpath = getwd()
}
### 相关卡值
AdjPvalueCutoff= 0.05; topNum = 40; height=16; width = 14; heatmap_breakNum = 2 ;
## the heatmap_annotation color
Tar_group ="#2080C3";The_remaining ="#89BCDF"    ## Tar_group ="#2080C3";The_remaining ="#89BCDF" 常用颜色
## the heatmap color
heatmap_col = c('#2080C3', '#f7f7f7', 'red') ##c('#2166ac', '#f7f7f7', 'red')) #E7E6E1 亮灰色
### 获取分组信息
subtype_list = Get_subtype_info(meta)

## ----------- using limma to conduct DE analysis ------------
for (i in subtype_list)
{
  print("-----------------------------------------------------")
  cat(i,"differential protein analysis is starting", sep = " ")
  print(i)
  # print("-----------------------------------------------------")
  Group = meta
  Group$type = ifelse(Group$type == i, i,"The_remaining" )  ## 表示亚型—1，和剩余其他亚型的比较
  group <- factor(Group$type)
  design <- model.matrix(~0 + group)
  # gsub("group","", colnames(design))
  colnames(design) <- c(i,"The_remaining") ## 更改名字
  rownames(design) <- colnames(exprSet)
  
  fit <- lmFit(exprSet, design)
  Contrasts = paste(i,"The_remaining",sep = "-")
  cont.matrix=makeContrasts(Contrasts,levels = design)
  fit2=contrasts.fit(fit,cont.matrix)
  fit2=eBayes(fit2)
  diff_gene <- topTable(fit2, adjust='BH', number=Inf, p.value=AdjPvalueCutoff )
  diff_gene$cluster <- ifelse(diff_gene$t > 0 , i, "The_remaining")
  # return(diff_gene)
  ##------ to plot we take top 20 pathways by logFC -------------------------------
  diff_gene <- sortByCol(diff_gene, c("logFC")) # ,asc = FALSE 升序排列
  topGene <- c(rev(rownames(diff_gene))[1:topNum],rownames(diff_gene)[1:topNum])
  mat <- exprSet[topGene, ]
  
  top_Gene = data.frame(topGene)
  topGeneset = data.frame(mat)
  
  #----------- crete dir  --------------------------
  outpath <-paste(getwd(),i ,sep = '/')
  if (file.exists(outpath)){
    print('Path already exists')
  } else {
    dir.create(file.path(outpath))
  }
  
  #----------- Print results objects to workspace  --------------------------
  
  write.csv(diff_gene, paste(getwd(), i,paste('diff_gene_results_pvalue',".csv",sep = ''),sep = '/'))
  write.csv(top_Gene,paste(getwd(), i, paste('top_Gene',".csv",sep = ''),sep = '/'))
  write.csv(topGeneset,paste(getwd(), i, paste('topGeneset_results',".csv",sep = ''),sep = '/'))
  
  ##-----------------------format pathway names ----------------------------------
  # rownames(mat) <- gsub("_", " ", rownames(mat))
  rownames(mat) <- gsub("KEGG ", "", rownames(mat))
  
  ##-----------------------format pathway names 大小写转换 -----------------------
  pathway_name = rownames(mat)
  # pathway_name = tolower(pathway_name) ## 将x中的字符全部转换为小写字母
  # pathway_name = capitalize(pathway_name) ## 将y中的字符全部转换为大写字母
  rownames(mat) = pathway_name
  if (TRUE) {
    
    Group$type = ifelse(Group$type != "The_remaining", "Tar_group", "The_remaining")
    topAnn = HeatmapAnnotation(df = Group[ ,"type", drop=F], 
                               col = list(type = c("Tar_group" = Tar_group,"The_remaining"= The_remaining)),
                               # col = c("#2080C3", "#89BCDF"),
                               annotation_height = unit(1, "cm"))
    
    heat.col <-  colorRamp2(c(-heatmap_breakNum, 0, heatmap_breakNum), heatmap_col)   #c('#2166ac', '#f7f7f7', 'red'))  ## 原始图例c('#2166ac', '#f7f7f7', '#b2182b'))
    ## ------  对数据求Z-score  -----------
    mat<- apply(mat, 1, Z_Score)
    mat = t(mat)
    ht <- Heatmap(mat, name="Z-score", col = heat.col, top_annotation = topAnn, 
                  cluster_rows = F, cluster_columns = F, show_column_names =F , show_row_names = TRUE,
                  row_names_side = "right")
    ## pdf oupput pathway
    pdf(paste(getwd(), i ,paste('top_20_pathway', '.pdf' ,sep = ''), sep = '/'), height=height, width = width)
    maxChar <- rownames(mat)[nchar(rownames(mat))==max(nchar(rownames(mat)))]
    
    padding <- unit.c(unit(2, "mm"), 
                      grobWidth(textGrob(maxChar))-unit(50, "mm"),
                      unit(c(2, 2), "mm"))
    draw(ht, padding = padding, merge_legends = TRUE)
    dev.off()
    
  }
  cat(i ,"differential protein analysis is Done", sep = " ","\n")
}

#=========================================================================
#       Panel=Fig4 - Heatmap of BAP1 mutation analysis
#=========================================================================

library(pheatmap)
library(RColorBrewer)
matrix=read.csv("/Users/ranpeng/OneDrive - mail.ecust.edu.cn/项目文件/CCA-project/Data/Fig4/Fig4_mutation_seelct.csv",header=T,
                skip =0,row.names = 1)
# matrix[is.na(matrix)] = 1
# matrix <- log2(matrix+1)
# matrix[matrix == 1] <- NA

#color for subtype
annotation_col = read.csv ('/Users/ranpeng/Desktop/CCA/Data/Fig6/iCC_eCC_anatation.csv',header = T,row.names = 1)
rownames(annotation_col) = colnames(matrix)
# ann_colors = list(type = c(iCC_N= "#48cae4", eCC_N = "#f48c06", iCC_T = "#0077b6" , eCC_T = "#dc2f02")) ## #0772BC 蓝色
# ann_colors = list(Subtype = c(S_I = "#F94F21",S_II = "#EFC000",S_III = "#00AEBA"))
chending_color = c(colorRampPalette(c("#1E90FF", "white"))(3),  ## #184D87  常用#1E90FF ##14417C ## 0772BC ## #395F90黑蓝   ## 淡蓝#1E90FF
                   colorRampPalette(c("white", "red"))(3) )

breaks = unique(c(seq(-2, 0, length = 4), 0, seq(0, 2, length = 4)))
## 选取基因进行可视化
# df = read.csv("/Users/ranpeng/Desktop/2021-3-7/KSEA_output/select.csv")
x <- read.delim(pipe("pbpaste"),row.names = 1) ## 读取剪切板的数据
gene = x$topGene
# gene = c("ROCK1","MAPK3")
data = matrix[gene,]
# data[is.na(data)] = 0
write.csv(data,"/Users/ranpeng/Desktop/2021-3-7/iCC_eCC_compare/phospho/select_pathway_2.csv")
# matrix = log2(matrix+1)
pheatmap(x, scale = "row",cluster_rows=0,cluster_cols=0,
         clustering_distance_cols = "correlation",fill = T, breaks=breaks,
         clustering_distance_rows = "correlation",border_color ="gray", na_col = "#EDEEEF",
         col=chending_color,show_rownames=T,show_colnames=F,display_numbers=F,
         width = 4.85,height = 5, fontsize_number=12, number_color="black",number_format="%.1f",
         annotation_col = annotation_col)#, annotation_colors = ann_colors) #, annotation_row = annotation_row)

#=========================================================================
#       Panel=Fig4 - Group boxplot for continuous numeric variables
#=========================================================================

# library
library(tidyverse)
library(viridis)
# ------------  load dataset ---------
data <- read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig4/iCC_subtype_mapping/immune_score.csv")
# Create dataset
# data <- data.frame(
#   individual=paste( "Mister ", seq(1,60), sep=""),
#   group=c( rep('A', 10), rep('B', 30), rep('C', 14), rep('D', 6)) ,
#   value1=sample( seq(10,100), 60, replace=T),
#   value2=sample( seq(10,100), 60, replace=T),
#   value3=sample( seq(10,100), 60, replace=T)
# )
# Transform data in a tidy format (long format) 数据溶化
data <- data %>% gather(key = "observation", value="value", -c(1,2)) 
Data = data
### ---------- 小提琴图结合箱线图用于表达数据 ------------
## 定义画图的格式
Group = observation; value = value; Attribute = Proteome_Subtype

P2 <- ggplot(Data, aes(x=observation, y=value,color=Proteome_Subtype)) + 
  geom_violin(trim=FALSE) + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  geom_boxplot(width=0.2,position=position_dodge(0.9),outlier.colour = "red")+ #绘制箱线图
  # scale_fill_manual(values = c("#0E9F87", "#3C5588"))+ #设置填充的颜色
  scale_color_manual(values=c("#0E9F87", "#3C5588")) +
  theme_bw()+ #背景变为白色
  stat_summary(fun.y=mean, geom="point", shape=23, size=2,position = position_dodge(width = 0.9)) + # violin plot with mean points
  theme(axis.text.x=element_text(angle=15,hjust = 1,colour="black",family="Times",size=20), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Arial",size=16,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Arial",size = 20,face="plain"), #设置y轴标题的字体属性
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.6), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="plain", family="Arial", colour="black",  #设置图例的子标题的字体属性
                                 size=16),
        legend.title=element_text(face="plain", family="Arial", colour="black", #设置图例的总标题的字体属性
                                  size=18),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  ylab("Value")+xlab("") + #设置x轴和y轴的标题
  # stat_compare_means(aes(group = Attribute), label = "p.format")  # 添加统计显著型水平
  stat_compare_means(aes(group = Proteome_Subtype),method = "wilcox",label = "p.format", label.y = 4.5)+        # Add global annova p-value
  stat_compare_means(aes(group = Proteome_Subtype),label = "p.signif", method = "t.test")  # Pairwise comparison against all

P2 

## ------------------ ggboxplot visualization ------------------
rm(list = ls())
#先加载包
library(ggplot2)
library(magrittr)
library(ggpubr)
library(cowplot)
library(tidyr)
data<-read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig4/imsig_score.csv",header=T)
data <- read.delim(pipe("pbpaste"),header = T)#,row.names = 1) ## 读取剪切板的数据
head(data)
compare_means(ARG1~ Type,  data = data, method = "anova")
# Bar plot of mean +/-seggbarplot
df = data[,3:length(data)]
df = as.data.frame(t(scale(t(df))))
df_scale = cbind(data[,1:2],df)

sig4 = subset(data,data$Proliferation < 1)

colnames(data)
# my_comparisons <- list( c("S_I", "S_II"), c("S_I", "S_III"), c("S_I", "S_III") )
b<-ggboxplot(data, x = "BAP1", y = "log2_TMB",combine = F,
             add = c("mean_se"),
             color = "BAP1" ,
             # fill = "Proteome_Subtype", alpha = 0.12,
             # order =  c("S_I","S_II","S_III"),
             palette = c("#FF9E29", "#0E65AD","#F94F21"),   #c("#FF9E29","#86AA00","#F94F21","#916CA0","#599BAD","#DBD289"), ##palette = c("#F94F21","#EFC000","#00AEBA"),
             width = 0.5 ,x.text.angle = 45)+
  rotate_x_text(angle=45,hjust=.5,vjust=1) +
  stat_compare_means(label.y = c()) +                # Global p-value
  # tat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  # stat_compare_means(label.y = 50) +    # Add global p-value
  # stat_compare_means(aes(group = Type),  label.y = c()) +
  labs(title="BAP1 Mutation", x="", y = "Mutation number(Log2)")  +
  geom_jitter(alpha= 0.5,aes(colour = BAP1),width = 0.15,height = 0) +
  theme_test()+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))  ## 控制坐标轴文本的角度
  theme(axis.text.x=element_text(angle=45,hjust = 1,colour="black",family="Arial",size=16), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Arial",size=16,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Arial",size = 20,face="plain"), #设置y轴标题的字体属性
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.8), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="plain", family="Arial", colour="black",  #设置图例的子标题的字体属性
                                 size=12),
        legend.title=element_text(face="plain", family="Arial", colour="black", #设置图例的总标题的字体属性
                                  size=12),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank(),
        # panel.border = element_rect(colour = "black", fill=NA, size=1)
        )  #不显示网格线
        # coord_flip() # 坐标轴转置
print(b)

## -------------- for-loop plot muiltiple boxplot ----------------
# a = log2(ToothGrowth[,7:length(ToothGrowth)])
a = t(scale(t((data[,7:length(data)]))))
df_scale = cbind(data[,1:6],a)
df_scale[df_scale == NA] = 0.00001
data = df_scale
# 需要先运行上面的步骤
b = list()
for(i in c(7:ncol(data)))
{
  gene_data<-data[,c(5,i)]
  # gene_data %>% drop_na() ## delletion NA
  b[[i]]<-ggboxplot(gene_data, x=colnames(gene_data)[1],
                    y = colnames(gene_data)[2],combine = F, add = c("mean_se"),
                    color = "Proteome_sub", fill = "Proteome_sub", alpha = 0.12,
                    order = c("S_I","S_II","S_III"),
                    palette = c("#0E9F87", "#3C5588","#F94F21"), ### c("#FF9E29","#86AA00","#F94F21","#916CA0","#599BAD","#DBD289"),
                    width = 0.5,x.text.angle = 45)+
    stat_compare_means(label.y = c()) +            # Global p-value
    # stat_compare_means(aes(group = Type)) +
    labs(title= colnames(gene_data)[2], x="", y = "Relative Expression")  +
    geom_jitter(alpha= 0.5,aes(colour = Proteome_sub),width = 0.15,height = 0) +
    theme_test() +
    theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))  ## 控制坐标轴文本的角度
  
}
## 生成多个图关键的一步
plot_grid(plotlist = b)

#=========================================================================
#       Panel=Fig4 - volcano plot of mutation overall survival
#=========================================================================

setwd("/Users/ranpeng/OneDrive - mail.ecust.edu.cn/项目文件/CCA-project/Figure/Fig")
library(ggplot2)
library(ggrepel)
genes <- read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig4/iCC_subtype_mapping/mutation_OS.csv", header = TRUE)
colnames(genes)
genes$Significant <- ifelse(genes$p < 0.05 & genes$HR > 1 , "Bad prognosis", 
                            ifelse(genes$p < 0.05 & genes$HR < 1 , "good prognosis","Not Sig"))
ggplot(genes, aes(x = log2(HR), y = -log10(p))) +
  geom_point(aes(color = Significant)) +
  scale_size_continuous(range = c(1,4))+
  scale_color_manual(values = c( "#EA686B","#7AA9CE","gray")) + # ("gray", "#7AA9CE")
  theme_bw(base_size = 12) + theme(legend.position = "right") +
  geom_text_repel(
    data = subset(genes, p < 0.05),
    aes(label = X),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle("Mutation correlation with Protein")

#=========================================================================
#     Panel=Fig4 - Fisher exact test of Significant mutation genes
#=========================================================================

library(dplyr)
df = read.csv('/Users/ranpeng/Desktop/CCA/Data/Fig4/iCC_subtype_mapping/significant_mutation_gene.csv',row.names = 1,
              skip = 1, header = T)
head(df)
df_LargeBD = df[,1:93]
df_SmallBD = df[,94:length(df)]

## 函数定义
# f<-function(x) {sum(x==1)}
# apply(A,2,f)

if(TRUE){
  ### ————————————— 基因突变计数  ——————————————————————————————-
  LargeBD_Mut_count = apply(df_LargeBD, 1, function(x) { length(which(x== 1))})
  LargeBD_WT_count= apply(df_LargeBD, 1, function(x){ length(which(x== 0))})
  SmallBD_Mut_count = apply(df_SmallBD, 1, function(x) { length(which(x== 1))})
  SmallBD_WT_count = apply(df_SmallBD, 1, function(x) { length(which(x== 0))})
  
  ### ————————————— 基因突变计数dataframe  ——————————————————————————————-
  LargeBD_fishier_count <- data.frame(LargeBD_Mut_count, LargeBD_WT_count,SmallBD_Mut_count,SmallBD_WT_count)
  
  # a = read.csv('/Users/ranpeng/Desktop/CCA-data/2021-1-18/HBV/HBV_Arm_fisher_test.csv')
  a = LargeBD_fishier_count
  gene_length = nrow(a)
  PV_result_ab = c()
  # PV_result_ac = c()
  # PV_result_bc = c()
  
  for (i in 1:gene_length) {
    ab<-matrix(c(a[i,1],a[i,2],a[i,3],a[i,4]),nrow = 2)
    y = fisher.test(ab)
    PV_result_ab <- c(PV_result_ab, y$p.value)
    
    # ac<-matrix(c(a[i,2],a[i,3],a[i,6],a[i,7]),nrow = 2)
    # y = fisher.test(ac)
    # PV_result_ac <- c(PV_result_ac, y$p.value)
    # 
    # bc<-matrix(c(a[i,4],a[i,5],a[i,6],a[i,7]),nrow = 2)
    # y = fisher.test(ac)
    # PV_result_bc <- c(PV_result_bc, y$p.value)
  }
  
  result_merge = a
  result_merge["PV_ab"] <-  PV_result_ab
  # result_merge["PV_ac"] <-  PV_result_ac
  # result_merge["PV_bc"] <-  PV_result_bc
  fisher_significant = subset(result_merge, result_merge[,5] < 0.05 )
  
}

df = read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig4/iCC_subtype_mapping/significant_mutation_gene.csv")
df$BAP1_y = ifelse(df$BAP1_y == 1, "Mut", "WT")
df$X6p22.2_Amp = ifelse(df$X6p22.2_Amp > 0.3, "Amp", 
                        ifelse(df$X6p22.2_Amp < -0.3,"Del","WT"))

df$X6p21.2_del = ifelse(df$X6p21.2_del > 0.3, "Amp", 
                        ifelse(df$X6p21.2_del < -0.3,"Del","WT"))
table(df$BAP1_y)

#=========================================================================
# Panel=Fig4 - Stack barplot for didplying distribution of categorical variables
#=========================================================================

library(vcd)
library(ggplot2)

data = read.csv('/Users/ranpeng/Desktop/CCA/Data/Fig4/BAP1_mutation/BAP1_mutation_barplot.csv')
colnames(data)
# data = df
## Calculate the counts
ggplot(data=data, mapping=aes(x=Proteome_Subtype,fill=BAP1))+
  geom_bar(stat="count",width=0.5,position='fill')+
  scale_fill_manual(values=c("#FF9E29","#86AA00","#F94F21","#916CA0","#599BAD","#DBD289"))+  ## '#999999','#E69F00'
  geom_text(stat='count',aes(label=..count..), color="white", size=3.5,position=position_fill(0.5))+
  theme_minimal() #+coord_flip()

## Calculate percentage
ggplot(data=data, mapping=aes(x=Proteome_Subtype,fill=BAP1))+
  geom_bar(stat="count",width=0.7,position='fill')+
  scale_fill_manual(values=c("#FF9E29","#86AA00","#F94F21","#916CA0","#599BAD","#DBD289"))+
  geom_text(stat='count',aes(label=scales::percent(..count../sum(..count..)))
            , color="white", size=3.5,position=position_fill(0.5))+
  # theme_minimal() #+coord_flip()
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle("BAP1 mutations across proteome subtypes")

## ------------- 统计突变的频率 --------------------------
df = read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig4/gene_mutation_results/CCA_CNAs_del.csv")
res = table(df$alteration,df$subtype)
fisher.test(res)

df_1 = data.frame(table(df$alteration,df$subtype))
ggplot(data=df_1, aes(x=Var2, y=Freq, fill=Var1)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=Freq), vjust=1.6, color="white",
            position = position_dodge(0.9), size=3.5)+
  scale_fill_brewer(palette="Paired")+
  theme_minimal() +
  xlab("") + ylab("Mutation Number")+
  theme(panel.border = element_blank(),# panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle("Mutation")

## ------------ mutation fisher -------------------------
df = read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig4/BAP1_mutation/fisher.csv")
res = table(df$Proteome_Subtype,df$BAP1)
fisher.test(res)

df = read.csv("/Users/ranpeng/Desktop/result-test1.csv")
df$cot = apply(df[,2:length(df)], 1, function(x)sum(is.na(x)))
filter = subset(df,df$cot<length(df)*0.7)
dim(filter)
