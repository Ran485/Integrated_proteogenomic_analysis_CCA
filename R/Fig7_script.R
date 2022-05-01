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
#      Panel=Fig7 - Heatmap xCell signatures and Pathways
#=========================================================================

library(pheatmap)
library(RColorBrewer)
# rm(list = ls())
matrix=read.csv("/Volumes/Samsung_T5/downloads/CCA/Data/Fig7/Xcell_signatures.csv",header=T,skip = 0,row.names = 1)
chending_color<-colorRampPalette(c("#0690F7","#FFFFFF","red"))(50)
# Generate annotations for rows and columns
annotation_col = read.csv ('/Volumes/Samsung_T5/downloads/CCA/Data/Fig7/Immune_subtypes_anatation.csv',header = T,row.names = 1)
# Clinical anatation colors
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
#      Panel=Fig7 - Survival analysis for immune subtypes
#=========================================================================
library(survival)
library(survminer)

data = read.csv("/Volumes/Samsung_T5/downloads/CCA/Data/Fig7/immune_cluster_index.csv")
colnames(data)
exprSet = data
sfit <- survfit(Surv(time, status) ~ CCP_Class , data = exprSet)
ggsurvplot(sfit, conf.int=F, pval=TRUE,risk.table = TRUE, risk.table.col = "strata",
           palette = c("#FF9E29", "#86AA00", "#F94F21", "#916CA0", "#599BAD", "#DBD289", "#FB9A99", "#A4D873", "#99CAE0"))

#=========================================================================
# Panel=Fig7 - ggboxplot for numerical continuous variables of immune signatures
#=========================================================================
rm(list = ls())
#先加载包
library(ggplot2)
library(magrittr)
library(ggpubr)
library(cowplot)
library(tidyr)

Immune_signatures<-read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig7/RNA/Xcell_immune_boxplot.csv",header=T)
Immune_signatures = Immune_signatures[-which(is.na(Immune_signatures$pro.B.cells)),]
head(Immune_signatures)
fcompare_means(ARG1~ Type,  data = Immune_signatures, method = "anova")
# a<-log(ToothGrowth[2:21],2)
# write.table(a,"log2(select).csv",sep=",")
# Bar plot of mean +/-seggbarplot
df = Immune_signatures[,3:length(Immune_signatures)]
df = as.data.frame(t(scale(t(df))))
df_scale = cbind(Immune_signatures[,1:2],df)

sig4 = subset(Immune_signatures,Immune_signatures$ImmuneScore_Xcell < 0.075)
Immune_signatures = data.frame(data)
colnames(Immune_signatures)
Immune_signatures$CD4 = as.numeric(Immune_signatures$CD4)
data = subset(Immune_signatures, CCP_Class == "Cluster_5" | CCP_Class == "Cluster_6" )#|CCP_Class == "Cluster_3" |CCP_Class == "Cluster_4" )

b<-ggboxplot(data, x = "CCP_Class", y = "ImmuneScore",combine = F,
             # add = c("mean_se"),
             color = "CCP_Class" ,
             fill = "CCP_Class", alpha = 0.12,
             order =  c("Cluster_1","Cluster_2","Cluster_3","Cluster_4","Cluster_5","Cluster_6"),
             palette = c("#FF9E29","#86AA00","#F94F21","#916CA0","#599BAD","#DBD289"), ##palette = c("#F94F21","#EFC000","#00AEBA"),
             width = 0.5 ,x.text.angle = 45)+
  rotate_x_text(angle=45,hjust=.5,vjust=1) +
  stat_compare_means(label.y = c()) +                                         # Global p-value
  # stat_compare_means(aes(group = Type),  label.y = c()) +
  labs(title="ImmuneScore", x="", y = "RNA_FPKM")  +
  geom_jitter(alpha= 0.7,aes(colour = CCP_Class),width = 0.15,height = 0) +
  theme_test()+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))  ## 控制坐标轴文本的角度
# coord_flip()
print(b)

## ------------ for-loop ----------------

Immune_signatures[Immune_signatures == 0.00001] = NA
a = log2(Immune_signatures[,4:length(Immune_signatures)]+1)
df_scale = cbind(Immune_signatures[,1:3],a)
df_scale[df_scale == NA] = 0.00001

# ToothGrowth = df_scale
# 需要先运行上面的步骤
for(i in c(4:ncol(df_scale)))
{
  gene_data<-df_scale[,c(3,i)]
  # gene_data %>% drop_na() ## delletion NA
  b[[i]]<-ggboxplot(gene_data, x=colnames(gene_data)[1],
                    y = colnames(gene_data)[2],combine = F, add = c("mean_se"),
                    color = "CCP_Class", fill = "CCP_Class", alpha = 0.12,
                    order = c("Cluster_1","Cluster_2","Cluster_3","Cluster_4","Cluster_5","Cluster_6"),
                    palette = c("#FF9E29","#86AA00","#F94F21","#916CA0","#599BAD","#DBD289"),
                    width = 0.5,x.text.angle = 45)+
    stat_compare_means(label.y = c()) +            # Global p-value
    # stat_compare_means(aes(group = Type)) +
    labs(title= colnames(gene_data)[2], x="Type", y = "Xcell Score")  +
    geom_jitter(alpha= 0.5,aes(colour = CCP_Class),width = 0.15,height = 0) +
    theme_test() +
    theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))  ## 控制坐标轴文本的角度
  
}
## 生成多个图关键的一步
plot_grid(plotlist = b)
ggarrange(a,a,a + rremove("x.text"),ncol = 2, nrow = 2 )

#=========================================================================
# Panel=Fig7 - Density plot for immune signatures
#=========================================================================
## scatter + density
data = read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig7/Neutrophil_cor.csv",row.names = 1)
library(cowplot) 
# data$type = ifelse(data$location == 1, "iCC", "eCC")
colnames(data)
# Main plot
pmain <- ggplot(data, aes(x = Neutrophils_RAW, y = CA199, color = CCP_Class),alpha = 0.6)+
  geom_point(alpha = 0.8)+
  ggpubr::color_palette("jco")
# Marginal densities along x axis
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = data, aes(x = Neutrophils_RAW, fill = CCP_Class),
               alpha = 0.7, size = 0.2)+
  ggpubr::fill_palette("jco")
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = data, aes(x = CA199, fill = CCP_Class),
               alpha = 0.7, size = 0.2)+
  coord_flip()+
  ggpubr::fill_palette("jco")
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
p <- ggdraw(p2)
p + theme(panel.grid.major =element_blank(),
          # panel.grid.minor = element_blank(),
          # panel.background = element_blank(),
          axis.line = element_line(colour = "black"))

p + classic()


#=========================================================================
# Panel=Fig7 - Scatter plot for Clinical and xCell signature overall survival
#=========================================================================
library(ggplot2)
library(ggrepel)

genes <- read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig7/Xcell_signature_cor/Neutrophils/focal_arm_CNV/Neutrophils_cor_CNV_thresh.csv", header = TRUE)
colnames(genes)
genes$Significant <- ifelse(genes$p_val < 0.05 & genes$cor_r > 0.2 , "Positive cor",
                            ifelse(genes$p_val < 0.05 & genes$cor_r < -0.2 , "Negtive cor","Not Sig"))
genes$size = ifelse(genes$Significant == "Positive cor" | genes$Significant == "Negtive cor" , abs(genes$cor_r), 0.2)

ggplot(genes, aes(x = cor_r , y = -log10(p_val), size = size)) +
  geom_point(aes(color = Significant)) +
  scale_size_continuous(range = c(2,4))+
  scale_color_manual(values = c("#7AA9CE","gray","#EA686B")) + # ("gray", "#7AA9CE") c( "#EA686B","#7AA9CE","gray")) c("#515151", "red"))
  theme_bw(base_size = 12) + theme(legend.position = "right") +
  geom_text_repel(
    data = subset(genes, Significant == "Positive cor" | Significant == "Negtive cor"),
    aes(label = gene),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle("Xcell signature cor with CA199")+
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
# scale_x_continuous(limits = c(-0.25, 0.3))+
# scale_y_continuous(breaks=seq(0, 4, 1))
# panel.border = element_rect(colour = "black", fill=NA, size=1)

#=========================================================================
#     Panel=Fig7 - Correlation analysis of clinical features
#=========================================================================

rm(list = ls())
#先加载包
library(ggplot2)
library(magrittr)
library(ggpubr)
library(cowplot)
library(tidyr)

data = read.csv('/Users/ranpeng/Desktop/CCA/Data/Fig7/Xcell_signature_cor/clinical_sort_os.csv')
data$Group <- as.factor(data$Group)
colnames(data)
df = data
df$PLT_log = log2(df$PLT)

p1 <-ggscatter(df, x = "PLT_log", 
               y = "GMP",
               color="#B2DF89",
               add = "reg.line") + 
  # geom_point()+
  # geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  labs(#title= paste(colnames(gene_data)[2],colnames(gene_data)[3],sep = '_'),
    x="PLT_log", y = "GMP") +
  # scale_color_brewer(palette="paired") +
  theme_minimal()+
  stat_cor(aes(color = Group ,alpha = 0.3),method = "pearson",alpha = 1) + ### 添加相关性 pearson spearman
  theme(#panel.border = element_blank(),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(), 
    #axis.line = element_line(colour = "black")
  )+
  ggtitle("TMB correlation with signature")

p1 + theme(panel.grid.major =element_blank(),
           # panel.grid.minor = element_blank(),
           panel.background = element_blank(),
           axis.line = element_line(colour = "black")
)

p1 + theme_define2

theme_define1 <- theme(panel.border = element_blank(),
                       # panel.grid.major = element_blank(),
                       panel.background = element_blank(), 
                       axis.line = element_line(colour = "black"))

theme_define2 <- theme(axis.text.x=element_text(angle=0,hjust = 0.6,colour="black",family="Arial",size=16), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
                       axis.title.x=element_text(angle=0,hjust = 0.5,colour="black",family="Arial",size=16),
                       axis.text.y=element_text(family="Arial",size=16,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
                       axis.title.y=element_text(family="Arial",size = 20,face="plain"), #设置y轴标题的字体属性
                       panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.8), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
                       legend.text=element_text(face="plain", family="Arial", colour="black",  #设置图例的子标题的字体属性
                                                size=12),
                       legend.title=element_text(face="plain", family="Arial", colour="black", #设置图例的总标题的字体属性
                                                 size=12),
                       panel.grid.major = element_blank(),   #不显示网格线
                       panel.grid.minor = element_blank())


## -------------- clincal info correlation with Xcell signature ----------
#计算两个matrix correlation
library(reshape)
library(Hmisc)
library(dplyr)
setwd('/Users/ranpeng/Desktop/Xcell_score/Xcell_score_clinical_correlation')
data1=read.csv("/Users/ranpeng/Desktop/Xcell_score/Xcell_score_clinical_correlation/Xcell_score_immune_cluster.csv",header=T,row.names=1)
data2=read.csv("/Users/ranpeng/Desktop/Xcell_score/Xcell_score_clinical_correlation/clionical_sort.csv",header=T,row.names=1)
dim(data1)
dim(data2)

ov_sample <- intersect(colnames(data1), colnames(data2))
matrix1 <- data1[, ov_sample]
matrix2 <- data2[, ov_sample]
#corrArray <- cor(t(matrix1), t(matrix2), method = "spearman")
corrp <-rcorr(t(matrix1), t(matrix2), type = c("spearman"))

p1<-corrp$P
r1<-corrp$r
dim(p1)
## 融化数据
b1<-p1[30:59,1:29]  ## 取右下角的矩阵data1 为CNV ;data2 为protein。
b2<-r1[30:59,1:29]  ## [row(data1)+1 : sum(row(data1)+row(data2)) , 1:row(data1)]
p_val=melt(b1)
cor =melt(b2)

## 给数据列名重新命名
colnames(p_val)[1]="clinical" ;colnames(p_val)[2]="Xcell_score";colnames(p_val)[3]="P_value"
colnames(cor)[1]="clinical" ;colnames(cor)[2]="Xcell_score";colnames(cor)[3]="cor_r"

merge = merge(p_val, cor, all = F)
# load('/Users/ranpeng/Desktop/CCA-data/2020-12-3/CNV_protein_correlation/raw_data/cor_merge_protein_2020_12_03.Rdata')
#merge_1 = inner_join(p_val,cor,by="protein")
# save(merge,file = "cor_merge_protein_2020_12_03.Rdata")
## 数据筛选 P<0.05 R>0.1
#c1<-as.matrix(merge,stringsAsFactors=FALSE)

newdata<-subset(merge, merge[,3] < 0.05 )#& abs(merge[,4]) > 0.3)
newdata<-subset(merge, merge[,3] < 0.05 & abs(merge[,4]) > 0.25)
write.csv(merge,'/Users/ranpeng/Desktop/CCA/Data/Fig7/Xcell_signature_cor/Xcell_score_immune_cluster_cor_result.csv')

## Xcell signature correlation
my_data = read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig7/RNA/RNA_immune_ICI_target_cor.csv",row.names = 1)
my_data = my_data[,-1]
# my_data = my_data[-which(is.na(my_data$ImmuneScore)), ]
res <- cor(my_data, method = "spearman") #"spearman" pearson
round(res, 2)
col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = res, col = col, symm = TRUE)

chending_color = c(colorRampPalette(c("blue", "white"))(0),  ## #184D87  常用#1E90FF ##14417C ## 0772BC ## #395F90 黑蓝   ## 淡蓝#1E90FF
                   colorRampPalette(c("white", "red"))(15) ) 

breaks = unique(c(seq(0, 0, length = 31), 0, seq(0.1, 0.8, length = 16)))
pheatmap(res, #scale = "row",
         cluster_rows=0,cluster_cols=0,
         clustering_distance_cols = "correlation",fill = T, breaks=breaks,
         clustering_distance_rows = "correlation",border_color ="Black", na_col = "#EDEEEF",
         col=chending_color,show_rownames=T,show_colnames=T,display_numbers= T,
         width = 4.85,height = 5, fontsize_number=18, number_color="black",number_format="%.2f")

#=========================================================================
#  Panel=Fig7 - Focal 1p36.33 Amplication cis-effect on RNA and Protein
#=========================================================================
#先加载包

library(ggplot2)
library(magrittr)
library(ggpubr)
library(cowplot)
library(tidyr)
ToothGrowth<-read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig7/RNA/arm_merge.csv",header=T)
colnames(ToothGrowth)

x_lab = "type"; y_lab = "CA199"
b<-ggboxplot(ToothGrowth, x = x_lab, y = y_lab ,combine = F,
             # add = c("mean_se"),
             color = x_lab ,
             fill = x_lab, alpha = 0.12,
             # order =  c("Cluster_1","Cluster_2","Cluster_3","Cluster_4","Cluster_5","Cluster_6"),
             palette = c("#FF9E29","#86AA00","#F94F21","#916CA0","#599BAD","#DBD289"), ##palette = c("#F94F21","#EFC000","#00AEBA"),
             width = 0.5 ,x.text.angle = 45)+
  rotate_x_text(angle=45,hjust=.5,vjust=1) +
  stat_compare_means(label.y = c()) +                                         # Global p-value
  # stat_compare_means(aes(group = Type),  label.y = c()) +
  labs(title = y_lab, x="", y = "Realtive Expression")  +
  geom_jitter(alpha= 0.7,aes(colour = type),width = 0.15,height = 0) +
  theme_test()+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))  ## 控制坐标轴文本的角度
# coord_flip()
print(b)


### ---------- 小提琴图结合箱线图用于表达数据 ------------
# ## 数据类型为长格式
# Attribute  Group     value
# 1         Amp CDK11A  9.559318
# 2         Amp CDK11A  5.154437
# 3         Amp CDK11A  5.539539
Data = read.csv("/Volumes/Samsung_T5/downloads/CCA/Data/Fig7/Neutrophils/Protein_boxplot20210605.csv")
## 过滤离群值
# Data = subset(Data,value <= 5)
## 运用盖帽法处理离群值
Data$value = ifelse(Data$value > 4, 4, Data$value)

P2<- ggplot(Data, aes(x=Group, y=value,color=Attribute)) + 
  # geom_violin(trim=FALSE) + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  geom_boxplot(width=0.6,position=position_dodge(0.9),outlier.colour = "red")+ #绘制箱线图
  # geom_jitter(alpha= 0.7,aes(x=Group, y=value,colour = Attribute),width = 0.15,height = 0) +
  # scale_fill_manual(values = c("#0E9F87", "#3C5588"))+ #设置填充的颜色
  scale_color_manual(values=c("#0E9F87", "#3C5588")) +
  theme_bw()+ #背景变为白色
  stat_summary(fun.y=mean, geom="point", shape=23, size=2,position = position_dodge(width = 0.9)) + # violin plot with mean points
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle("Cis efffect on Protein")+
  theme(axis.text.x=element_text(angle=30,hjust = 0.6,colour="black",family="Arial",size=16), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.title.x=element_text(angle=0,hjust = 0.5,colour="black",family="Arial",size=16),
        axis.text.y=element_text(family="Arial",size=16,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Arial",size = 20,face="plain"), #设置y轴标题的字体属性
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.8), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="plain", family="Arial", colour="black",  #设置图例的子标题的字体属性
                                 size=12),
        legend.title=element_text(face="plain", family="Arial", colour="black", #设置图例的总标题的字体属性
                                  size=12),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  ylab("Value")+xlab("") + #设置x轴和y轴的标题
  # stat_compare_means(aes(group = Attribute), label = "p.format")  # 添加统计显著型水平
  stat_compare_means(aes(group = Attribute),method = "wilcox",label = "p.format", label.y = 4.5) +       # Add global annova p-value
  stat_compare_means(aes(group = Attribute),label = "p.signif", method = "wilcox")  # Pairwise comparison against all

P2 

