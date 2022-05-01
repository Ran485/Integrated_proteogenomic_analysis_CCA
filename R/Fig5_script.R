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
#  Panel=Fig5 - Analysis of differentially expressed proteins
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
  input_path <- '/Users/ranpeng/Desktop/CCA/Data/Fig5/HBV/data_wash/Fig_5_HBV_merge_phosphosite_20210317.csv'
  anno_path <- '/Users/ranpeng/Desktop/CCA/Data/Fig5/HBV/HBV_anatation.csv'
  input_data <- read.csv(input_path, header=T, row.names=1, stringsAsFactors = F,)
  meta <- read.csv(anno_path)
  meta$id = gsub("#",".",meta$id)
  exprSet <- as.matrix(input_data)
  ## 数据处理 ，磷酸化数据需要进行log转化
  exprSet[is.na(exprSet)] = 0
  exprSet<- log2(exprSet + 1)
  ## 如果存在空值进行下面的处理
  # tm1<-testmatrix[,-which(apply(testmatrix,2,function(x)all(is.na(x))))]  ## 删除所有是空值的列
  # tm2<-tm1[-which(apply(testmatrix,1,function(x)all(is.na(x)))),] ## 删除所有是空值的行
  # exprSet<-exprSet[,-which(apply(exprSet,2,function(x) all(is.na(x))))]
  # meta <- meta[which(meta$id %in% colnames(exprSet)),]
  ## 数据存储路径
  setwd('/Users/ranpeng/Desktop/CCA/Data/Fig5/HBV/results')
  outpath = getwd()
}
### 相关卡值
AdjPvalueCutoff= 1; topNum = 40; height=16; width = 14; heatmap_breakNum = 2 ;
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
# Panel=Fig5 - Analysis of HBV-associated mutation genes overall survivlal 
#=========================================================================

setwd("/Users/ranpeng/Desktop/CCA/Figure/Fig5")
library(ggplot2)
library(ggrepel)
genes <- read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig5/HBV_mutation_os.csv", header = TRUE)
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
  ggtitle("Significant Mutation")

#=========================================================================
# Panel=Fig5 - Analysis of focal CNV and protein correlation
#=========================================================================

library(Hmisc)
library(dplyr)
## anatation
data1 = read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig5/HBV/HBV_Amp_Del_TMB.csv", header = T, row.names = 1)
data1 = t(data1)
## protein or phospho
data2 = read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig1/RNA_protein_cor/protein_del30%.csv",header = T,row.names = 1)
data2 = t(data2)

## -------- 数据处理 ，磷酸化数据需要进行log转化 -------
data2[is.na(data2)] = 0
# data2<- log2(data2 + 1)

result<-c()
for(i in colnames(data2)){
  newdata = data.frame(cbind(data1[,i],data2[,i]))
  ## ————————  相关性检验  ——————————————
  res <- rcorr(as.matrix(newdata),type=c("spearman")) #pearson
  cor_r = res$r[2]; p_val = res$P[2]; gene = i
  result.linshi<-cbind(gene,cor_r,p_val)
  result<-rbind(result,result.linshi)
}
result = data.frame(result)
write.csv(result,"/Users/ranpeng/Desktop/CCA/Data/Fig5/HBV/RNA_protein_cor/results/HBV_Protein_cor_spearman.csv")
result$cor_r = as.numeric(result$cor_r);result$p_val = as.numeric(result$p_val)
result_filter <- filter(result, p_val < 0.05 & abs(cor_r) > 0.3)
(result_pos_cor <- filter(result, p_val < 0.05 & cor_r > 0.3))
(result_neg_cor <- filter(result, p_val < 0.05 & cor_r < -0.3))
write.csv(result_neg_cor,"//Users/ranpeng/Desktop/CCA/Data/Fig5/HBV/RNA_protein_cor/results/HBV_Protein_cor_spearman_molecular_r0.3_neg.csv")

#=========================================================================
#      Panel=Fig5 - ggboxplot for HBV, TMB and focal arm events
#=========================================================================

rm(list = ls())
#先加载包
library(ggplot2)
library(magrittr)
library(ggpubr)
library(cowplot)
library(tidyr)

data<-read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig5/HBV/HBV_select_arm_os.csv",header=T)
head(data)
# Bar plot of mean +/-seggbarplot
df = data[,3:length(data)]
df = as.data.frame(t(scale(t(df))))
df_scale = cbind(data[,1:2],df)

sig4 = subset(data,data$ImmuneScore_Xcell < 0.075)
data$Amp_log = log2(data$AMP)
colnames(data)
b<-ggboxplot(data, x = "HBV_status", y = "X6p22.2.._Amplification",combine = F,
             add = c("mean_se"),
             color = "HBV_status" ,
             # fill = "Proteome_Subtype", alpha = 0.12,
             order =  c("HBV","Non_HBV"),
             palette = c("#0E9F87", "#3C5588"),   #c("#FF9E29","#86AA00","#F94F21","#916CA0","#599BAD","#DBD289"), ##palette = c("#F94F21","#EFC000","#00AEBA"),
             width = 0.5 ,x.text.angle = 45)+
  rotate_x_text(angle=45,hjust=.5,vjust=1) +
  stat_compare_means(label.y = c()) +                      # Global p-value
  # stat_compare_means(aes(group = Type),  label.y = c()) +
  labs(title="Del_log", x="", y = "Total Deletion")  +
  geom_jitter(alpha= 0.5,aes(colour = HBV_status),width = 0.15,height = 0) +
  theme_test()+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))  ## 控制坐标轴文本的角度
# coord_flip()
print(b)

#=========================================================================
#   Panel=Fig5 - Survival analysis of HBV-associated focal arm events
#=========================================================================
## ------------ Panel=Fig5i - arm events os ---------------
data  = read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig5/HBV/HBV_select_arm_os.csv")
colnames(data)

data$X6p22.2.._Amplification_1 = ifelse(data$X6p22.2.._Amplification> 0.3, "Amp", 
                                        ifelse(data$X6p22.2.._Amplification< -0.3, "Del", "WT"))
data$X5q33.1.._Deletion_1 = ifelse(data$X5q33.1.._Deletion > 0.3, "Amp",#"WT")
                                   ifelse(data$X5q33.1.._Deletion< -0.3, "Del", "WT"))
data$X4p_1 = ifelse(data$X4p >0.3, "X4p_1_Amp","X4p_1_WT")

test = table(data$HBV_status,data$X6p_Amp)
# test[1,1] = 1
fisher.test(test)
write.csv(data,"/Users/ranpeng/Desktop/CCA/Data/Fig5/HBV/HBV_select_arm_os.csv")
test = table(data$X1q21.2_1,data$X6q14.3_1)

group = ifelse(data$X6p22.2.._Amplification >0.3, "Amp",
               ifelse(data$X6p22.2.._Amplification< -0.3, "Del", "WT"))
# data = na.omit(data)
df= data
df$time = as.numeric(df$time)
library(survival)
library(survminer)
sfit <- survfit(Surv(time, status) ~ group , data = df)

ggsurvplot(sfit, conf.int=F, pval=TRUE,risk.table = TRUE, risk.table.col = "strata",
           palette = c("#FF9E29", "#86AA00", "#F94F21", "#916CA0", "#599BAD", "#DBD289", "#FB9A99", "#A4D873", "#99CAE0"))

