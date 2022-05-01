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
#  Panel=Fig6 - Analysis of differentially expressed proteins
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
  input_path <- '/Users/ranpeng/Desktop/CCA/Data/Fig6/data_wash/Fig6_merge_iCCA_eCCA_protein_20210317.csv'
  anno_path <- '/Users/ranpeng/Desktop/CCA/Data/Fig6/iCCA_eCCA_anatation.csv'
  input_data <- read.csv(input_path, header=T, row.names=1, stringsAsFactors = F,)
  meta <- read.csv(anno_path)
  meta$id = gsub("#",".",meta$id)
  exprSet <- as.matrix(input_data)
  ## 数据处理 ，磷酸化数据需要进行log转化
  # exprSet[is.na(exprSet)] = 0
  # exprSet<- log2(exprSet + 1)
  ## 如果存在空值进行下面的处理
  # tm1<-testmatrix[,-which(apply(testmatrix,2,function(x)all(is.na(x))))]  ## 删除所有是空值的列
  # tm2<-tm1[-which(apply(testmatrix,1,function(x)all(is.na(x)))),] ## 删除所有是空值的行
  # exprSet<-exprSet[,-which(apply(exprSet,2,function(x) all(is.na(x))))]
  # meta <- meta[which(meta$id %in% colnames(exprSet)),]
  ## 数据处理 ，磷酸化数据需要进行log转化
  exprSet[is.na(exprSet)] = 0
  # exprSet<- log2(exprSet + 1)
  ## 数据存储路径
  setwd('/Users/ranpeng/Desktop/CCA/Data/Fig6/results')
  outpath = getwd()
}
### 相关卡值
AdjPvalueCutoff= 1; topNum = 10; height=16; width = 14; heatmap_breakNum = 2 ;
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
    ht
    draw(ht, padding = padding, merge_legends = TRUE)
    dev.off()
    
  }
  cat(i ,"differential protein analysis is Done", sep = " ","\n")
}


#=========================================================================
#  Panel=Fig6 - Volcano of focal arm events and protein cis-effect
#=========================================================================

setwd("/Users/ranpeng/Desktop/CCA/Figure/Fig6")
library(ggplot2)
library(ggrepel)
genes <- read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig6/Fig6_focal_arm_cor.csv", header = TRUE)
colnames(genes)
genes$Significant <- ifelse(genes$p-value < 0.05 & genes$Correlation > 0.2 , "Positive correaltion", 
                            ifelse(genes$p-value < 0.05 & genes$Correlation < 0.2 , "Negtive correaltion","Not Sig"))
ggplot(genes, aes(x = Correlation, y = -log10(p-value))) +
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
  ggtitle("Significant cis-effect correlation")

#=========================================================================
# Panel=Fig6 - Lasso and survival analysis for pathology of iCCA and eCCA
#=========================================================================

## ----------- lasso and survival--------------------
#filter potential useful sig genes using univariate cox regression.
library('DESeq2')
library('survival')
library('survminer')
library('dplyr')
library('glmnet')
library('ggplot2')
library('GGally')
library('rms')
library('survivalROC')
library('plotROC')

# 1.对差异基因进行批量单因素COX回归。
uni_cox_in_bulk <- function(gene_list, survival_info_df){
  library('survival')
  gene_list <- gsub(gene_list, pattern = '-', replacement = '_')
  uni_cox <- function(single_gene){
    formula <- as.formula(paste0('Surv(overall_survival, censoring_status)~', single_gene))
    surv_uni_cox <- summary(coxph(formula, data = df))
    ph_hypothesis_p <- cox.zph(coxph(formula, data = df))$table[1,3]
    if (surv_uni_cox$coefficients[,5]<0.05 & ph_hypothesis_p>0.05){  #get the pvalue
      single_cox_report <- data.frame('uni_cox_sig_genes'=single_gene,
                                      'beta'=surv_uni_cox$coefficients[,1],
                                      'Hazard_Ratio'=exp(surv_uni_cox$coefficients[,1]),
                                      'z_pvalue'=surv_uni_cox$coefficients[,5],
                                      'Wald_pvalue'=as.numeric(surv_uni_cox$waldtest[3]),
                                      'Likelihood_pvalue'=as.numeric(surv_uni_cox$logtest[3]))
      single_cox_report
    }
  }
  uni_cox_list <- lapply(gene_list, uni_cox)
  do.call(rbind, uni_cox_list)
}
df = read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig6/biomarker/iCC_os.csv",header = T,row.names = 1)
# df = t(df)
##  去除空值
df<- df[-which(is.na(df$time)),]
colnames(df)[1:3] = c("overall_survival","time_day","censoring_status")
uni_cox_df <- uni_cox_in_bulk(gene_list = colnames(df)[4:length(df)], survival_info_df = df)

# 2.Lasso回归筛选具有代表性的变量。

#about glmnet: x should be in format of matrix, and time&status in y should be in double format.
x <- as.matrix(df[,colnames(df)[4:length(df)]])
y <- df[,c('overall_survival', 'censoring_status')]
names(y) <- c('time', 'status')
y$time <- as.double(y$time)
y$status <- as.double(y$status)
y <- as.matrix(survival::Surv(y$time, y$status))
lasso_fit <- cv.glmnet(x, y, family='cox', type.measure = 'deviance')
coefficient <- coef(lasso_fit, s=lasso_fit$lambda.min)
Active.Index <- which(as.numeric(coefficient) != 0)
active.coefficients <- as.numeric(coefficient)[Active.Index]
sig_gene_multi_cox <- rownames(coefficient)[Active.Index]

# 3.利用筛选出来的少数candidates建立多因素COX回归模型。
#perform the multi-variates cox regression using qualified genes.
formula_for_multivariate <- as.formula(paste0('Surv(overall_survival, censoring_status)~', paste(sig_gene_multi_cox, sep = '', collapse = '+')))
multi_variate_cox <- coxph(formula_for_multivariate, data = df)
#check if variances are supported by PH hypothesis.
ph_hypo_multi <- cox.zph(multi_variate_cox)
#The last row of the table records the test results on the GLOBAL model. Delete it.
ph_hypo_table <- ph_hypo_multi$table[-nrow(ph_hypo_multi$table),]
#Remove variances not supported by ph hypothesis and perform the 2nd regression.
formula_for_multivariate <- as.formula(paste0('Surv(overall_survival, censoring_status)~', paste(rownames(ph_hypo_table)[ph_hypo_table[,3]>0.05], sep = '', collapse = '+')))
multi_variate_cox_2 <- coxph(formula_for_multivariate, data = df)
# ALAD 移除
# 关于检测变量间的共线性，需要综合考虑：变量间的相关系数<0.5和vif平方根<2均被提出是可行的检验方法。
# 我们进行相关系数，vif的计算，并进行相关性矩阵的可视化。

# 4.check the co-linearity between samples.
correlation <- cor(df[,rownames(ph_hypo_table)[ph_hypo_table[,3]>0.05]], method = 'pearson')
library('GGally')
ggpairs(df[,rownames(ph_hypo_table)[ph_hypo_table[,3]>0.05]], 
        axisLabels = 'show')+
  theme_bw()+
  theme(panel.background = element_rect(colour = 'black', size=1, fill = 'white'), 
        panel.grid = element_blank())
library('rms')
vif <- rms::vif(multi_variate_cox_2)
#Some people said if the square root of VIF >2, they might be co-linear.
sqrt(vif) < 2

# -------- 接下来用森林图可视化这个模型 -------------
ggforest(model = multi_variate_cox_2, data = df, main = 'Hazard ratios of candidate genes', fontsize = 1)

# 5.多方面评价cox模型的预测能力。
C_index <- multi_variate_cox_2$concordance['concordance']
if(C_index >= 0.9){
  print('High accuracy')
}else{ 
  if(C_index < 0.9 & C_index >= 0.7){
    print('Medium accuracy')
  }else{
    print('Low accuracy')
  }
}


#=========================================================================
#     Panel=Fig6 - Fisher exact test for clinical features
#=========================================================================

df = read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig6/iCC_eCC_Arm/eCC_ana.csv")
# df = df[-which(df$Fluke_infection == "Unknown"),]
res = table(df$X2q11.2,df$X1q21.1)
fisher.test(res)
## wilcoxl.test
wilcox.test(MicroenvironmentScore_Xcell~Pathology,data=df) ##三种方法总体疗效有差别
## kruskal.test
kruskal.test(months~group,data=krusdata) ##三种方法总体疗效有差别

a = read.delim(pipe("pbpaste"),header = T,row.names = 1) ## 读取剪切板的数据
gene_length = nrow(a)
PV_result_ab = c()
# PV_result_ac = c()
# PV_result_bc = c()

for (i in 1:gene_length) {
  ab<-matrix(c(a[i,1],a[i,2],a[i,3],a[i,4]),nrow = 2)
  y = fisher.test(ab)
  PV_result_ab <- c(PV_result_ab, y$p.value)
}

result_merge = a
result_merge["PV_ab"] <-  PV_result_ab
# result_merge["PV_ac"] <-  PV_result_ac
# result_merge["PV_bc"] <-  PV_result_bc
fisher_significant = subset(result_merge, result_merge[,5] < 0.05 )
write.csv(fisher_significant,"/Users/ranpeng/Desktop/CCA/Data/Fig6/results/eCC_fisher_test_significance.csv")

#=========================================================================
#     Panel=Fig6 - Merge ssGSEA pathway results
#=========================================================================

if(TRUE){
  library("plyr")  #加载获取rbind.fill函数
  library("dplyr")  
  ##-----------------------load gmt file -----------------------------------------
  KEGG <- read.table('/Users/ranpeng/Desktop/CCA/Data/Fig6/GSVA_pathway/KEGG/ssGSEA_results.csv',sep = ",",header = T)
  Hall_mark <- read.table('/Users/ranpeng/Desktop/CCA/Data/Fig6/GSVA_pathway/Hall_mark/ssGSEA_results.csv',sep = ",",header = T)
  Reactome <- read.table('/Users/ranpeng/Desktop/CCA/Data/Fig6/GSVA_pathway/Reactome/ssGSEA_results.csv',sep = ",",header = T)
  
  #不等长数据合并
  list1<-list()
  list1[[1]]=KEGG
  list1[[2]]=Hall_mark
  list1[[3]]=Reactome
  # list1[[4]]=GO_BP
  merge_pathway = do.call(rbind.fill,list1)
  rownames(merge_pathway) = merge_pathway[,1]
  merge_pathway = merge_pathway[,-1]
}

pathway_name = read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig6/GSVA_pathway/select_pathway.csv")
data = merge_pathway [pathway_name$pathway,]









