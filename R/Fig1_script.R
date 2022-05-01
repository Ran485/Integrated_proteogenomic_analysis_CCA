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

#=======================================================================================
# load the R packages
library(pheatmap)
library(RColorBrewer)

rm(list=ls(all=TRUE)) ## 清空数据

# change the contents of variable baseDir to the root analysis folder 
baseDir <- "~/Desktop/Integrated_multiomics_analysis_CCA"
# load directory structure downloaded from github site
source(paste0(baseDir,"/R/load_directory_structure.R"))
#=========================================================================
# Panel=Fig1A - Cohort design and Multi-platform profiling of tumour samples
#=========================================================================

matrix=read.csv("../../Data/Fig1/Figure-1a/heatmap_all.csv",header=T,row.names = 1)
#chending_color<-colorRampPalette(c("#0689EA","#FFFFFF","#962436"))(20)
chending_color<-colorRampPalette(c("#0690F7","#FFFFFF","red"))(50)
#chending_color<-colorRampPalette(rev(c("#EE0000","#white","#436eee")))(20)
#chending_color <- colorRampPalette(brewer.pal(8, "PiYG"))(50)
# Generate annotations for rows and columns
annotation_col = read.csv ('/Users/ranpeng/Desktop/CCA/Data/Fig1/Figure-1a/Figure-1a-ana.csv',header = T,row.names = 1)
# annotation_row = read.csv ('/Users/ranpeng/Desktop/胆管癌/data/2020-8-21/annotation_row_all.csv',header = T,row.names = 1)

ann_colors = list(
  Patient = c(select= "#078DC6"),
  Proteome_N = c(select = "#7CD8AE", Miss = "#E4E4DF"),
  Proteome_T = c(select = "#0EA55D", Miss = "#E4E4DF"),
  Phosphoproteome_N = c(select = "#F78566", Miss = "#E4E4DF"),
  Phosphoproteome_T = c(select = "#E55625", Miss = "#E4E4DF"),
  WES_N = c(select = "#CDB3DD", Miss = "#E4E4DF"),
  WES_T = c(select = "#B076D6", Miss = "#E4E4DF")
)

pheatmap(matrix,cluster_rows=0,cluster_cols=0,
         clustering_distance_cols = "correlation",fill = T,
         clustering_distance_rows = "correlation",border_color ="white", na_col = "white",
         col=chending_color,show_rownames=T,show_colnames=F,display_numbers=F,
         width = 4.85,height = 5,fontsize_number=12,number_color="black",number_format="%.1f",
         annotation_col = annotation_col, annotation_colors = ann_colors)


#=========================================================================
# Panel=stable1 - Cox Proportional-Hazards Model for clinical Features
#=========================================================================
library("survival")
library("survminer")
library ("readxl")
library("tidyverse")
library("dplyr")

# data("lung")
# head(lung)

# load clinical metadata (Supplementary Table 1)
metadata  <- data.frame(read_excel(paste0(dataDir,"Fig1/clinical_information_cox_model.xlsx"),sheet = 1))
# Converting CA-199, CEA, Tumor size, ALT, Tbil features to categorical variable
metadata$Preoperative.CA19.9.U.ml. = ifelse(metadata$Preoperative.CA19.9.U.ml. >= 37, 1, 
                                            ifelse(metadata$Preoperative.CA19.9.U.ml. < 37, 0, NA))
metadata$Preoperative.CEA.ng.ml. = ifelse(metadata$Preoperative.CEA.ng.ml. >= 5, 1, 
                                          ifelse(metadata$Preoperative.CEA.ng.ml. < 5, 0, NA))
metadata$Tumor.size..cm. = ifelse(metadata$Tumor.size..cm. >= 5, 1, 
                                  ifelse(metadata$Tumor.size..cm. < 5, 0, NA))
metadata$ALT..aminoleucine.transferase..U.L. = ifelse(metadata$ALT..aminoleucine.transferase..U.L. >= 40, 1, 
                                                      ifelse(metadata$ALT..aminoleucine.transferase..U.L. < 40, 0, NA))
metadata$TB..total.bilirubin..µmol.L.  = ifelse(metadata$TB..total.bilirubin..µmol.L.  >= 20.4, 1, 
                                                ifelse(metadata$TB..total.bilirubin..µmol.L.  < 20.4, 0, NA))

# ----------------  Univariate Cox regression ----------------------
exclude_col <- c("Patient.No.", "Absolute.tumor.cellularity","Histology.tumor.cellularity","Survival.1.dead.0.alive.",
                 "Tumor..T..proteomic.ID","Overall.survival.month.","Overall.survival.day.")
covariates <- colnames(metadata)
covariates <- setdiff(covariates, exclude_col)
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(Overall.survival.day., Survival.1.dead.0.alive.)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = metadata)})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
res = as.data.frame(res)
res
# write.csv(res, paste0(outTableDir, "Fig1_Univariate Cox regression_converted.csv"))
# ----------------  Multivariate Cox regression analysis ----------------------
# A Cox regression of time to death on the time-constant covariates is specified as follow
res$p.value =as.numeric(res$p.value)
res_p0.05 = filter(res, p.value < 0.05)
res.cox <- coxph(Surv(Overall.survival.day., Survival.1.dead.0.alive.) ~ 
                   Preoperative.CA19.9.U.ml. +
                   Preoperative.CEA.ng.ml. +
                   ALB..albumin...g.L. +
                   TB..total.bilirubin..µmol.L. +
                   ALT..aminoleucine.transferase..U.L. +
                   PT..prothrombin.time..s. +
                   # Histology + 
                   iCCA.classification +
                   HBV.infection.status +
                   HBcAb +
                   Perineural.invasion..PNI. +
                   TNM.Stage +
                   Grade +
                   Proteome.subtype +
                   Immune.cluster,
                   # Proteome.subtype11 +
                   # CCP_Class11,
                   data =  metadata)
summary(res.cox)

#=========================================================================
# Panel=stable1 - Baseline characteristics and Demographics of the CCA cohort.
#=========================================================================
## Loading configuration files
source(paste0(baseDir,"/R/clinical_triliner_table_config.R"))
library("table1")
metadata_triliner = data.frame(read_excel(paste0(dataDir,"Fig1/Clinical_information_triliner_table.xlsx"),sheet = 1))
colnames(metadata_triliner)
str(metadata_triliner)
# Setting factor
metadata_triliner$Survival.1.dead.0.alive. <- factor(metadata_triliner$Survival.1.dead.0.alive.,levels = c(0,1),labels = c('alive','dead'))
metadata_triliner$Phosphoproteome.cluster <- factor(metadata_triliner$Phosphoproteome.cluster,levels = c(1,2,3),labels = c('Cluster_1','Cluster_2','Cluster_3'))
a = metadata_triliner

# Set up units
units(a$Overall.survival.month.) <- 'months'
units(a$Tumor.size..cm.) <- 'cm'
# units(a$Medical_History_Year) <- 'years'
# units(a$Residential_cycle) <- 'weeks'
# units(a$Height) <- 'cm'
# units(a$Weight) <- 'Kg'


# Set up lables
label(a$Overall.survival.month.) <- "survival"
label(a$Tumor.size..cm.) <- "Tumor size"
# label(a$wt.loss) <- "体重减低"
# label(a$sex) <- "性别"
# label(a$normdata) <- "正态数据"
# label(a$exp) <- "指数分布数据"
# table1(~status+age+wt.loss+normdata|sex,data = a)
# table1(~status+age+wt.loss+normdata|sex,data = a,render = rd,render.continuous=my.render.cont)


# Original function plot the chart
table1(~ Gender+
         Age+
         Histology+
         Overall.survival.month.+ 
         iCCA.classification+
         Tumor.size..cm.+
         TNM.Stage+
         TNM.Stage+
         Grade+
         HBV.infection.status+ HBsAg + HBcAb+
         Preoperative.AFP.ng.mL.+
         Preoperative.CA19.9.U.ml.+
         Preoperative.CEA.ng.ml.+
         ALB..albumin...g.L.+
         TB..total.bilirubin..µmol.L.+
         ALT..aminoleucine.transferase..U.L.+
         γ.GT..γ.glutamyltransferase..U.L.+
         TBA..total.bile.acid..µmol.L.+
         X..PLT..platelet.counts...10.9.L.+
         PT..prothrombin.time..s.+
         Cholelithiasis+
         Jaundice+
         Type.II.diabetes +
         Presence.of.fluke.infection+
         Perineural.invasion..PNI.+
         Lymphovascular.invasion..LVI.+
         Proteome.subtype+
         Phosphoproteome.cluster+
         CCP_Class+
         Absolute.tumor.cellularity+
         Histology.tumor.cellularity,data = a, 
         topclass="Rtable1-zebra",render = rd)

# Plota the diagram between iCCA and eCCA
# Sample stratification by `Histology` clinical feature
# table1(~status+age+wt.loss+normdata+exp|sex,data = a,render = rd,
#        droplevels = T,overall = F,
#        render.continuous=my.render.cont,
#        #关于新增列的名字可以自己在list()列表里面改
#        extra.col=list(`Normality`=normality,`homogeneity of var`=hom.var,`P.method`=pmethod,`P.value`=p))

table1(~  Gender+
          Age+
          # Histology+
          Overall.survival.month.+ 
          iCCA.classification+
          Tumor.size..cm.+
          TNM.Stage+
          TNM.Stage+
          Grade+
          HBV.infection.status+ HBsAg + HBcAb+
          Preoperative.AFP.ng.mL.+
          Preoperative.CA19.9.U.ml.+
          Preoperative.CEA.ng.ml.+
          ALB..albumin...g.L.+
          TB..total.bilirubin..µmol.L.+
          ALT..aminoleucine.transferase..U.L.+
          γ.GT..γ.glutamyltransferase..U.L.+
          TBA..total.bile.acid..µmol.L.+
          X..PLT..platelet.counts...10.9.L.+
          PT..prothrombin.time..s.+
          Cholelithiasis+
          Jaundice+
          Type.II.diabetes +
          Presence.of.fluke.infection+
          Perineural.invasion..PNI.+
          Lymphovascular.invasion..LVI.+
          Proteome.subtype+
          Phosphoproteome.cluster+
          CCP_Class+
          Absolute.tumor.cellularity+
          Histology.tumor.cellularity|Histology,
          data = a,
          topclass="Rtable1-zebra",render = rd,
          droplevels = T,overall = F,
          render.continuous=my.render.cont,
          #关于新增列的名字可以自己在list()列表里面改
          extra.col=list(`P.method`=pmethod,`P.value`=p))

## ================ Panel=Fig1B - gene mutation oncoprint ========================

## --------------- Synonymous and Non_Synonymous mutation stacked plot ------------------
library(plyr)

df = read.csv("/Users/ranpeng/Desktop/Integrated_multiomics_analysis_CCA/data/Fig1/fig1b/mutation_barplot_input.csv", row.names = 1)
colnames(df)
df$count = ifelse(df$count >= 800, 800, df$count)
# Create the barplot
ggplot(data=df, aes(x=num, y=count, fill=mutation_type)) +
  geom_bar(stat="identity")+
  # geom_text(aes(y=label_ypos, label=len), vjust=1.6, 
  #           color="white", size=3.5)+
  scale_fill_brewer(palette="Paired")+
  # xlim(0,25) +
  ylim(0,800) +
  labs(title="Significantly mutated genes", x="Mutation frequency (%)", y = "-Log10 (FDR)")+
  # theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  #       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(angle=0,hjust = 0.6,colour="black",family="Arial",size=16), 
        axis.title.x=element_text(angle=0,hjust = 0.5,colour="black",family="Arial",size=20),
        axis.text.y=element_text(family="Arial",size=16,face="plain"),
        axis.title.y=element_text(family="Arial",size = 20,face="plain"), 
        # panel.border = element_blank(),
        axis.line = element_line(colour = "black",size=0.8),
        legend.text=element_text(face="plain", family="Arial", colour="black", size=12),
        legend.title=element_text(face="plain", family="Arial", colour="black", size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

## --------------- Clinical information heatmap ------------------
library(pheatmap)
library(RColorBrewer)
# rm(list = ls())
matrix=read.csv("/Users/ranpeng/Desktop/Integrated_multiomics_analysis_CCA/data/Fig1/fig1b/CCA_CNV_plot_20220326_del.csv",header=T,skip = 0,row.names = 1)
#chending_color<-colorRampPalette(c("#3090D1","#FFFFFF","red"))(50)
chending_color<-colorRampPalette(c("#0690F7","#FFFFFF","red"))(50)
#chending_color <- colorRampPalette(brewer.pal(8, "PiYG"))(50)
# Generate annotations for rows and columns
# annotation_col = read.csv ('/Users/ranpeng/Desktop/CCA/Data/Fig1/Figure-1b/figure1_anatation.csv',header = T,row.names = 1)
annotation_col = read.csv ('/Users/ranpeng/Desktop/Integrated_multiomics_analysis_CCA/data/Fig1/fig1b/Fig1b_clinical_feature_anatation_20220326.csv',header = T,row.names = 1)
# ann_colors = list(type = c(iCC_N= "#48cae4", eCC_N = "#f48c06", iCC_T = "#0077b6" , eCC_T = "#dc2f02"))
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


# Specify colors
ann_colors = list(
  Patient = c(select= "#078DC6"),
  Proteome_N = c(select = "#7CD8AE", Miss = "#E4E4DF"),
  Proteome_T = c(select = "#0EA55D", Miss = "#E4E4DF"),
  Phosphoproteome_N = c(select = "#F78566", Miss = "#E4E4DF"),
  Phosphoproteome_T = c(select = "#E55625", Miss = "#E4E4DF"),
  WES_N = c(select = "#CDB3DD", Miss = "#E4E4DF"),
  WES_T = c(select = "#B076D6", Miss = "#E4E4DF")
)

## --------------- oncoprint_mutation ----------------------

library(ComplexHeatmap)
library(RColorBrewer)

onco_mat = read.csv('/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/Fig1/fig1b/onco_matrix(2022_03_26).csv',header = T, row.names = 1)
# metadata_sorted = read.csv('/Users/ranpeng/Desktop/2020-12-15/metadata_sorted.csv',header = T)
sample_order<- colnames(onco_mat)

# ================ oncoplot colors ========================

# col = c("Missense_Mutation" = "#EE1D24",
#         "In_Frame_Del" = "#F8941D",
#         "Nonsense_Mutation" = "#FCED68",
#         "Frame_Shift_Del" = "#ABD373",
#         "Splice_Site" = "#39B54A",
#         "Frame_Shift_Ins" = "#00A99D",
#         "In_Frame_Ins" = "#0272BC" ,
#         "Multi_Hit" = "#92278F")

# col = c("Missense_Mutation" = "#0085C7",
#         "In_Frame_Del" = "#02A595",
#         "Nonsense_Mutation" = "#34AD40",
#         "Frame_Shift_Del" = "#A3D357",
#         "Splice_Site" = "#F4E55B",
#         "Frame_Shift_Ins" = "#F8941D",
#         "In_Frame_Ins" = "#F04124" ,
#         "Multi_Hit" = "#92278F")

col = c("Missense_Mutation" = "#EE1D24",
        "In_Frame_Del" = "#F8941D",
        "Nonsense_Mutation" = "#FCED68",
        "Frame_Shift_Del" = "#ABD373",
        "Splice_Site" = "#39B54A",
        "Frame_Shift_Ins" = "#00A99D",
        "In_Frame_Ins" = "#0272BC" ,
        "Multi_Hit" = "#92278F",
        "FGFR2_fusion" = "#64D8C2")

# just for demonstration
alter_fun = list(
  background = alter_graphic("rect", fill = "#EFF8FD"),	#E5E5E5
  Missense_Mutation = alter_graphic("rect", fill = col["Missense_Mutation"]),
  In_Frame_Del = alter_graphic("rect", fill = col["In_Frame_Del"]),
  Nonsense_Mutation = alter_graphic("rect", height = 1, fill = col["Nonsense_Mutation"]),
  Frame_Shift_Del = alter_graphic("rect", fill = col["Frame_Shift_Del"]),	
  Splice_Site = alter_graphic("rect", fill = col["Splice_Site"]),
  Frame_Shift_Ins = alter_graphic("rect", fill = col["Frame_Shift_Ins"]),
  In_Frame_Ins = alter_graphic("rect", height = 1, fill = col["In_Frame_Ins"]),
  Multi_Hit = alter_graphic("rect", height = 1, fill = col["Multi_Hit"]),
  FGFR2_fusion = alter_graphic("rect", height = 1, fill = col["FGFR2_fusion"])
  # NA_1 = alter_graphic("rect", height = 0.4, fill = col["NA_1"])
  
)
column_title = "OncoPrint of the CCA"
heatmap_legend_param = list(title = "Alternations", at = c("Missense_Mutation", "In_Frame_Del", "Nonsense_Mutation","Frame_Shift_Del", "Splice_Site", "Frame_Shift_Ins","In_Frame_Ins", "Multi_Hit", "FGFR2_fusion"), 
                            labels = c("Missense_Mutation", "In_Frame_Del", "Nonsense_Mutation","Frame_Shift_Del", "Splice_Site", "Frame_Shift_Ins","In_Frame_Ins", "Multi_Hit", "FGFR2_fusion"))
mix.hp<-oncoPrint(onco_mat,
                  alter_fun = alter_fun, col = col, 
                  remove_empty_columns = F, remove_empty_rows = F, #column_order  = sample_order,
                  pct_side = "right", row_names_side = "left", #row_order  = 1:nrow(onco_mat),
                  column_title = column_title, heatmap_legend_param = heatmap_legend_param, show_column_names = F)
draw(mix.hp, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")

## --------------- oncoprint significant CNA ----------------------

onco_CNA_mat = read.csv('/Users/ranpeng/Desktop/Integrated_multiomics_analysis_CCA/data/Fig1/fig1b/CCA_CNV_oncoplot_20220326.csv',header = T, row.names = 1)
# metadata_sorted = read.csv('/Users/ranpeng/Desktop/2020-12-15/metadata_sorted.csv',header = T)
sample_order<- colnames(onco_CNA_mat)
row_order<- rownames(onco_CNA_mat)
# ================ oncoplot colors ========================

col = c("Amp" = "#EA1B2A",
        "gain" = "#FC7581",
        "loss" = "#4E9CF9",
        "Deletion" = "#1B6DCE")

# just for demonstration
alter_fun = list(
  background = alter_graphic("rect", fill = "#EFF8FD"),	#E5E5E5
  Amp = alter_graphic("rect", fill = col["Amp"]),
  gain = alter_graphic("rect", fill = col["gain"]),
  loss = alter_graphic("rect", height = 1, fill = col["loss"]),
  Deletion = alter_graphic("rect", fill = col["Deletion"])
  )

column_title = "CNA OncoPrint of the CCA"
heatmap_legend_param = list(title = "Alternations", at = c("Amp", "gain", "Nonsense_Mutation","loss", "Deletion"), 
                            labels = c("Amp", "gain", "Nonsense_Mutation","loss", "Deletion"))
mix.hp<-oncoPrint(onco_CNA_mat,
                  alter_fun = alter_fun, col = col, 
                  remove_empty_columns = F, remove_empty_rows = F, column_order  = sample_order,
                  row_order = row_order,
                  pct_side = "right", row_names_side = "left", #row_order  = 1:nrow(onco_mat),
                  column_title = column_title, heatmap_legend_param = heatmap_legend_param, show_column_names = F)
draw(mix.hp, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")


## --------------- Panel=Fig1F - Mutation signature -----------------------

#BiocManager::install("sigminer", dependencies = TRUE)
## 加载包
library(sigminer)
library(NMF)
library(BSgenome.Hsapiens.UCSC.hg19)
## 数据输入
#laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools", mustWork = TRUE)
#laml <- read_maf(maf = laml.maf)
# setwd('/Users/ranpeng/Desktop/CCA/Data/Fig1/Figure-1f')
laml.maf = 'Somatic.snv_indel.maf'
laml <- read_maf(maf = laml.maf)
#laml <- read.maf(maf = laml.maf, clinicalData = laml.clin, gisticAllLesionsFile="all_lesions.conf_90.txt", gisticAmpGenesFile="amp_genes.conf_90.txt", gisticDelGenesFile="del_genes.conf_90.txt", gisticScoresFile="scores.gistic")
samplesummary<-getSampleSummary(laml)


## 生成突变分类矩阵
#使用 sig_tally() 对突变进行归类整理，针对 MAF 对象，支持设定 mode 为 ‘SBS’, ‘DBS’, ‘ID’ 以及 ‘ALL’

mats <- mt_tally <- sig_tally(
  laml,
  ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
  useSyn = TRUE,
  mode = "ALL"
)
str(mats, max.level = 1)

# 针对整个数据集的分类就可以画图，Signatures 其实就是它的拆分。
show_catalogue(mt_tally$SBS_96 %>% t(), mode = "SBS", style = "cosmic")
show_catalogue(mt_tally$SBS_6 %>% t(), mode = "SBS", style = "cosmic")

## ------------ 估计 Signature 数目 --------------------
# 这一步实际上是多次运行 NMF，查看一些指标的变化，用于后续确定提取多少个 Signatures
# 运行最好30-50次可以得到稳定结果
est <- sig_estimate(mt_tally$SBS_96, range = 2:8, nrun = 50, pConstant = 1e-9, verbose = TRUE)

show_sig_number_survey2(est$survey)
show_sig_number_survey(est$survey)
## --------------- 提取 signatures 和可视化 ----------------
sigs <- sig_extract(mt_tally$SBS_96, n_sig = 4, nrun = 200, pConstant = 1e-9)
# 生成的是一个带 Signature 类信息的列表：
str(sigs, max.level = 1)
## 很多信息存在里面，可以自己提取自己感兴趣的信息
p <- show_sig_profile(sigs, mode = "SBS", style = "cosmic")
p

## 计算下它与 COSMIC signatures 的相似性，评估病因，对于 SBS 有 2 个 COSMIC 数据库
## 版本 legacy（30 个，目前最常用的） 与 SBS v3
get_sig_similarity(sigs, sig_db = "legacy")
sim <- get_sig_similarity(sigs, sig_db = "SBS")
add_labels(p, x = 0.72, y = 0.15, y_end = 0.9, labels = sim, n_label = 5)

## 自动提取 signatures
## SignatureAnalyzer 提供的变异 NMF 提供了自动提取 Sigantures 的功能，无需要自己判断提取的 signature 数目，这个可以通过 sig_auto_extract() 实现
sigs2 <- sig_auto_extract(mt_tally$SBS_96, nrun = 200)

## 绘图展示
p <- show_sig_profile(sigs2, mode = "SBS", style = "cosmic")
p

## ————————————————————  Signature 活动图谱  ——————————————————————————
# sigminer 提供绝对和相对两种 Signature 活动度值。
get_sig_exposure(sigs2) %>% head()
get_sig_exposure(sigs2, type = "relative") %>% head()
show_sig_exposure(sigs, rm_space = TRUE, style = "cosmic")

## -------------- 根据已知的 Signatures 提取活动度 --------------------
examp_fit <- sig_fit(mt_tally$SBS_96[1:5, ] %>% t(), sig_index = "ALL", sig_db = "legacy")
head(examp_fit)

## 图形绘制
show_sig_fit(examp_fit, palette = NULL) + ggpubr::rotate_x_text()

## 设置散点图，绘制单样本结果
show_sig_fit(examp_fit,
             palette = NULL, plot_fun = "scatter",
             signatures = paste0("COSMIC_", c(1, 2, 4, 6, 19))
) + ggpubr::rotate_x_text()

## --------------- Panel=Fig1G - the indentify number of T/N  -----------------------
rm(list = ls()) ## 清除数据
library("ggpubr")
library("ggplot2")

# Load data
df = read.csv('/Users/ranpeng/Desktop/CCA/Data/Fig1/Figure-1g/NP_MBR.csv',header = T)
head(df)

sp <- ggscatter(df, x = "Exp", y = "number",
                color = "Type", 
                #palette = "jco",
                # palette = c("blue","red"),
                palette = c("#3A86FF", "#FF5400"), # 自定义颜色
                size = 1.5, 
                alpha = 0.6, ## 调整透明度
                #xlab = "Paired tumor and NAT samples (n = 200) with unpaired tumor samples (n = 21)", 
                #ylab = "Number of protein identifications",
                add = "loess", 
                conf.int = TRUE,
                add.params = list(size = 0.7, linetype = 'dashed',alpha = 0.3)
) +
  border() + ## 添加边框
  geom_line(aes(group = Exp), color = "gray", alpha = 0.3)+  ### 添加 sample 连线
  stat_cor(aes(color = Type,alpha = 0.3),method = "spearman",alpha = 1) ### 添加相关性
#geom_smooth(aes(x = Exp,y = number,color=Type), alpha = 0.3)

sp
# sp + xlim(min(dt$A, 0)*1.2, max(dt$A)*1.2)
sp  + scale_y_continuous(limits = c(7500,10500)) ## 调整坐标轴的标度范围
sp + scale_x_continuous(breaks = c(4000,8000))


## --------------- Panel=sFig1C - the QC of protein profiling  -----------------------
# Density_T/N
rm=(list=ls())
for (file in 1:1){
  data = read.csv(paste("0",firN[file,],".csv",sep=''),header = T, row.names = 1)
  data <- density(as.matrix(log2(data)))
  plot(data,xlim=c(5,40),ylim=c(0,0.20),col="#3A86FF")
}
for (file in 2:7){
  data = read.csv(paste("0",firN[file,],".csv",sep=''),header = T, row.names = 1)
  data <- density(as.matrix(log2(data)))
  lines(data,col="#3A86FF")
}
for (file in 9:199){
  data = read.csv(paste("0",firN[file,],".csv",sep=''),header = T, row.names = 1)
  data <- density(as.matrix(log2(data)))
  lines(data,col="#3A86FF")
}

for (file in 1:220){
  data = read.csv(paste("0",firT[file,],".csv",sep=''),header = T, row.names = 1)
  data[data==0]=1e-5
  data <- density(as.matrix(log2(data)))
  lines(data,col="#FF5400")
}

for (file in 9:199){
  data = read.csv(paste("0",firN[file,],".csv",sep=''),header = T, row.names = 1)
  data[data==0]=1e-5
  data <- density(as.matrix(log2(data)))
  lines(data,col="#FF5400")
}

## --------------- Panel=sFig1D/F - the QC of protein profiling  -----------------------
# Load package

#install.packages("mnornt")
library("PerformanceAnalytics")
library("psych")
# Load data
my_data = read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig1/sFigure-1d/S-figure-1d_phospho_293T.csv" ,header = T, row.names = 1)
head(my_data)
data = log2(my_data[,2:13])
#data = abs(data)
#write.table(data,"293_log2.csv",row.names=FALSE,col.names=TRUE,sep=",")
pairs.panels(data, 
             hist.col="gray", 
             show.points=TRUE, 
             stars=FALSE, 
             gap=0.05, 
             pch=".", 
             ellipses=FALSE, 
             scale= FALSE,
             jiggle=TRUE,
             factor=2,
             main="Scatter plots and Spearmans correlation coefficients for 293T", 
             col="#ADFF2F", 
             pty="m", 
             font=2)

## --------------- Panel=sFig1E - the QC of protein profiling  -----------------------

library(tidyverse)
library(dplyr)

df = read_csv('/Users/ranpeng/Desktop/CCA/Data/Fig1/sFigure-1e/S-figure-1e.csv')
df[df == 0.00001] <- NA ## 替换数据框中所有的0值为NA
df %>% summarise_all(funs(sum(!is.na(.)))) ## 对每列进行计数
df_1 <- log2(df[2:249])

#p <- list()
plot(x = NA , y= NA ,xlab = 'Rank of quantified proteins',
     ylab = 'Quantative protein intensity[Log10]',
     xlim = c(0,7000), ylim = c(-10,20))

for (i in c(1:10)) 
{
  
  df_2 <-cbind(df[1],df_1[i]) ## 按列合并
  del_na <-df_2 %>% na.omit(df_2)  ## 去掉含有空值的列
  del_na <- arrange(del_na,desc(del_na[2]))  ## 从大到小排序
  n<- count(del_na) ## 计数
  del_na <- mutate(del_na, numbers = seq(1: n[[1]])) 
  points(x = del_na[[3]], y= del_na[[2]] ,xlab = 'Rank of quantified proteins',
         ylab = 'Quantative protein intensity[Log10]',
         pch=2, cex=0.5, col="#69b3a2" )
  ## pch= 2 调整形状
}

setwd("/Users/ranpeng/Desktop/CCA/Data/Fig1/sFigure-1e/results")
pdf(file="myplot.pdf")

## ggplot2
df_1[5]
df_2 <-cbind(df[1],df_1[5])
del_na <-df_2 %>% na.omit(df_2)
del_na <- arrange(del_na,desc(del_na[2]))
n<- count(del_na)
del_na <-mutate(del_na, numbers = seq(1: n[[1]]))
ggplot(data = del_na, aes(x = numbers , y= del_na[[2]])) + geom_point(size=1)
plot(x = del_na[[3]] , y= del_na[[2]])

## --------------- Panel=sFig1G - maftools    -----------------------
library(maftools)
# set path to MAF file
setwd('/Volumes/Samsung_T5/downloads/CCA/Data/Fig1/maf')
laml.maf = 'Somatic.snv_indel.maf'
laml.clin = 'clin_genomic_mapping.tsv'
laml = read.maf(maf = laml.maf, clinicalData = laml.clin)
# laml <- read.maf(maf = laml.maf, clinicalData = laml.clin, gisticAllLesionsFile = "all_lesions.conf_95.txt", gisticAmpGenesFile = "amp_genes.conf_95.txt", gisticDelGenesFile = "del_genes.conf_95.txt", gisticScoresFile = "scores.gistic")
samplesummary<-getSampleSummary(laml)

## ---------- Plotting MAF summary --------------
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE,top = 20)
## 可以将MAF文件的gene ，sample的 summary 的信息，输出到laml前缀的summary文件
# write.mafSummary(maf = laml, basename = 'laml')
## Drawing oncoplots
#oncoplot for top ten mutated genes.
oncoplot(maf = laml, top = 20, draw_titv = T)

## One can use any colors, here in this example color palette from RColorBrewer package is used
vc_cols = RColorBrewer::brewer.pal(n = 12, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)
vc_cols =c("#FC4963","#5AADFF","#66CDAA","#2B69C9","#CD0000","#86C4BA","#8ACC83","#6ACC5C","#E84A5F","#FBD46D","#96BB7C","#E1F48F")
names(vc_cols) =c('Amp','Del','In_Frame_Del', 'Missense_Mutation','In_Frame_Ins', 'Multi_Hit', 'Frame_Shift_Ins', 'Frame_Shift_Del','Complex_Event','Nonsense_Mutation')
print(vc_cols)

# oncoplot(maf = laml, colors = vc_cols, top = 15500, draw_titv = F,keepGeneOrder = TRUE,writeMatrix = TRUE)
# pdf('oncoplot_clin_top120_noamp.pdf',width=10,height=10) #'KRAS','FGFR3','SYNE1', 'TP53','LRP1B',
oncoplot(maf = laml, 
         genes = c('ARID1A','TP53', 'KMT2C','DPCR1','BAP1', 'ARID2', 
                   'KMT2A','MGA','KMT2D','APC','IDH1', 'IDH2', 'SMAD4',
                   'NCOR1','ATM','KRAS','PML','MSH6', 'EPHA2','TET2','PBRM1','RBM47','NF1','FGFR2'),
         showTumorSampleBarcodes=F,colors = vc_cols,keepGeneOrder = TRUE,writeMatrix = TRUE,
         sortByMutation = T, removeNonMutated = F,
         sortByAnnotation=TRUE, borderCol = "white", bgCol = "#F4F4F4", draw_titv = T,
         clinicalFeatures = c('Pathology','Perineural_invasion','Proteome_Subtype'))
dev.off()

oncoplot(maf=laml, 
         genes = c("TP53","PTEN ","IDH1","IDH2","KRAS","BRAF","PIK3CA","ARAF ","SMAD4","FGFR1–3","EGFR","BAP1"), 
         showTumorSampleBarcodes=TRUE,colors = vc_cols,keepGeneOrder = TRUE,writeMatrix = TRUE,
         sortByAnnotation=TRUE,
         borderCol=NULL)
## Oncostrip
oncostrip(maf = laml, genes = c('EPHA2','CDKN2B', 'CDKN2A'))

##  Transitions and Transversions
laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = laml.titv)

## Lollipop plots for amino acid changes
#lollipop plot for DNMT3A, which is one of the most frequent mutated gene in Leukemia.
lollipopPlot(maf = laml, gene = 'FGFR2', AACol = 'AAchange', showMutationRate = TRUE)

## Labelling points
lollipopPlot(maf = laml, gene = 'DNMT3A', showDomainLabel = FALSE, labelPos = 882)

## Rainfall plots

brca <- system.file("extdata", "brca.maf.gz", package = "maftools")
brca = read.maf(maf = brca, verbose = FALSE)
rainfallPlot(maf = brca, detectChangePoints = TRUE, pointSize = 0.6)

## 将突变载量与TCGA队列进行比较
laml.mutload = tcgaCompare(maf = laml, cohortName = 'CCA')

## Plotting VAF
plotVaf(maf = laml, vafCol = 'i_TumorVAF_WU')

## Genecloud
geneCloud(input = laml, minMut = 3)

## ———————————— 处理拷贝数数据 ——————————————————

## 8.1.读取并汇总gistic输出文件
all.lesions <- system.file("extdata", "all_lesions.conf_99.txt", package = "maftools")
amp.genes <- system.file("extdata", "amp_genes.conf_99.txt", package = "maftools")
del.genes <- system.file("extdata", "del_genes.conf_99.txt", package = "maftools")
scores.gis <- system.file("extdata", "scores.gistic", package = "maftools")

laml.gistic = readGistic( gisticAllLesionsFile="all_lesions.conf_95.txt", gisticAmpGenesFile="amp_genes.conf_95.txt", gisticDelGenesFile="del_genes.conf_95.txt", gisticScoresFile="scores.gistic", isTCGA = TRUE)

#GISTIC object
laml.gistic

## 8.2.1基因组图
gisticChromPlot(gistic = laml.gistic,fdrCutOff = 0.25,ref.build = "hg19",cytobandOffset = 0.01, markBands = "all",mutGenes = c("ARID1A", "DNMT3A", "RUNX1"))

## 8.2.2 Bubble plot
gisticBubblePlot(gistic = laml.gistic,fdrCutOff = 0.01,txtSize = 6)

## 8.2.3 oncoplot
vc_cols_1 =c("#FF788B","#1C95F9")
names(vc_cols_1) = c('Amp', 'Del')
pdf('/Users/ranpeng/Desktop/胆管癌/data/2020-10-4/gisticOncoPlot.pdf',width=8,height=5.5)
gisticOncoPlot(gistic = laml.gistic, clinicalData = getClinicalData(x = laml),colors = vc_cols_1,
               borderCol = "white", bgCol = "#F4F4F4", ## bgCol = "#F4F4F4",
               clinicalFeatures = c('Pathology','Perineural_invasion','Proteome_Subtype'), sortByAnnotation = TRUE, top = 20)
dev.off()
## 8.2.4 Visualizing CBS segments
tcga.ab.009.seg <- system.file("extdata", "TCGA.AB.3009.hg19.seg.txt", package = "maftools")
plotCBSsegments(cbsFile = tcga.ab.009.seg)

## ——————————————  数据分析  ————————————————
#exclusive/co-occurance event analysis on top 10 mutated genes. 
pdf('/Users/ranpeng/Desktop/胆管癌/data/2020-10-4/somaticInteractions.pdf',width=9,height=8.5)
somaticInteractions(maf = laml, genes = c('OBSCN','VPS13D','ARID1A','LAMA2','TERT', 'MUC17', 'MUC16','KRAS','FGFR3','SYNE1', 'TP53','ERBB3','LRP1B', 'FLG','BAP1',
                                          'MUC22','SYNE2','IDH1', 'SMAD4','PBRM1','FGFR2','ARID1B','IDH2','PIK3CA','BRAF','EGFR','ALB'), pvalue = c(0.05, 0.1))
dev.off()
## 基于位置聚类的癌症驱动基因检测
laml.sig = oncodrive(maf = laml,  minMut = 5,AACol = 'AAchange', pvalMethod = 'zscore',bgEstimate = TRUE, ignoreGenes = NULL)
head(laml.sig)
write.csv(laml.sig,'/Users/ranpeng/Desktop/Integrated_multiomics_analysis_CCA/output_tables/oncodriveCLUST_Mutsig.csv')
plotOncodrive(res = laml.sig, fdrCutOff = 0.25, useFraction = TRUE)

## 9.3 Adding and summarizing pfam domains
laml.pfam = pfamDomains(maf = laml, AACol = 'AAchange', top = 10)

## 9.4 Pan-Cancer comparison
#MutsigCV results for TCGA-AML

laml.mutsig <- system.file("extdata", "LAML_sig_genes.txt.gz", package = "maftools")
pancanComparison(mutsigResults = laml.mutsig, qval = 0.1, cohortName = 'LAML', inputSampleSize = 200, label = 1)

## 生存分析
#Survival analysis based on grouping of DNMT3A mutation status
mafSurvival(maf = laml, genes = c('MUC5B',"TP53"), time = 'time', col = c("maroon", "royalblue"),Status = 'status', isTCGA = F)

## Predict genesets associated with survival

#Using top 20 mutated genes to identify a set of genes (of size 2) to predict poor prognostic groups
prog_geneset = survGroup(maf = laml, top = 20, geneSetSize = 2, time = "time", Status = "status", verbose = FALSE)
print(prog_geneset)
# write.csv(prog_geneset,'/Users/ranpeng/Desktop/胆管癌/data/2020-10-4/genesets_associated_with_survival.csv')
## 利用函数mafSurvGroup可以绘制上述结果的KM曲线
mafSurvGroup(maf = laml, geneSet = c("TTN", "TP53"), time = "time", Status = "status")

## -------------- 对比两个队列(MAF)  -----------------

#Primary APL MAF
# primary.apl = system.file("extdata", "APL_primary.maf.gz", package = "maftools")
iCC.maf = 'group1_all.maf'
iCC = read.maf(maf = iCC.maf)
#Relapse APL MAF
# relapse.apl = system.file("extdata", "APL_relapse.maf.gz", package = "maftools")
eCC.maf = 'group2_all.maf'
eCC = read.maf(maf = eCC.maf)

##  Forest plots 森林图
#Considering only genes which are mutated in at-least in 5 samples in one of the cohort to avoid bias due to genes mutated in single sample.
pt_vs_rt <- mafCompare(m1 = iCC, m2 = eCC, m1Name = 'iCCA', m2Name = 'eCCA', minMut = 3)
print(pt_vs_rt)
forestPlot(mafCompareRes = pt_vs_rt, pVal = 0.1, color = c('royalblue', 'maroon'), geneFontSize = 0.8)
write.csv(pt_vs_rt$results,'/Users/ranpeng/Desktop/2020-12-15/Forest_plots.csv')
## Co-onco plots
genes = c("BAP1", "SVEP1", "KMT2B",'KMT2D','IDH1', "PML", "IQGAP1", "NRAS","SEC24B", "CAND1", "XRCC5", "THAP9" ,"RNF43")
coOncoplot(m1 = iCC, m2 = eCC, m1Name = 'iCCA', m2Name = 'eCCA', genes = genes, removeNonMutated = F,borderCol = "white", bgCol = "#F4F4F4")

## Lollipop plot-2棒棒糖图
lollipopPlot2(m1 = primary.apl, m2 = relapse.apl, gene = "PML", AACol1 = "amino_acid_change", AACol2 = "amino_acid_change", m1_name = "Primary", m2_name = "Relapse")

## ---------------  临床富集分析  ---------------
fab.ce = clinicalEnrichment(maf = laml, clinicalFeature = 'FAB_classification')

#Results are returned as a list. Significant associations p-value < 0.05
fab.ce$groupwise_comparision[p_value < 0.05]

plotEnrichmentResults(enrich_res = fab.ce, pVal = 0.05)


## 药物与基因的相互作用
dgi = drugInteractions(maf = laml, fontSize = 0.75)

dnmt3a.dgi = drugInteractions(genes = "MUC16", drugs = TRUE)
#> Number of claimed drugs for given genes:
#>      Gene N
#> 1: DNMT3A 7
#Printing selected columns.
dnmt3a.dgi[,.(Gene, interaction_types, drug_name, drug_claim_name)]

## 9.9致癌信号通路
OncogenicPathways(maf = laml)
## 抑癌基因用红色表示，癌基因用蓝色字体表示
PlotOncogenicPathways(maf = laml, pathways = "Hippo",fullPathway = FALSE)


## 9.11 突变Signatures
#Requires BSgenome object
library(BSgenome.Hsapiens.UCSC.hg19, quietly = TRUE)
laml.tnm = trinucleotideMatrix(maf = laml, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")

## 9.11.1 APOBEC Enrichment estimation.
plotApobecDiff(tnm = laml.tnm, maf = laml, pVal = 0.2)

## 9.11.3 Signature分析
library('NMF')
laml.sign = estimateSignatures(mat = laml.tnm, nTry = 6)
plotCophenetic(res = laml.sign)

laml.sig = extractSignatures(mat = laml.tnm, n = 3)
## 将检测到的signatures与COSMIC数据库中的已知signatures进行比较。
#Compate against original 30 signatures 
laml.og30.cosm = compareSignatures(nmfRes = laml.sig, sig_db = "legacy")
#Compate against updated version3 60 signatures 
laml.v3.cosm = compareSignatures(nmfRes = laml.sig, sig_db = "SBS")

## 检测到的signatures与验证过的signatures的相似性比较。
library('pheatmap')
pheatmap::pheatmap(mat = laml.og30.cosm$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures")

## Finally plot signatures
maftools::plotSignatures(nmfRes = laml.sig, title_size = 0.8, sig_db = "SBS")

## 9.11.4 Signature enrichment analysis

#  library("barplot3d")
#  #Visualize first signature
#  sig1 = laml.sig$signatures[,1]
#  barplot3d::legoplot3d(contextdata = sig1, labels = FALSE, scalexy = 0.01, sixcolors = "sanger", alpha = 0.5)

laml.se = signatureEnrichment(maf = laml, sig_res = laml.sig)

## 上述结果可进行和临床结果相似的可视化操作。
plotEnrichmentResults(enrich_res = laml.se, pVal = 0.05)

#=========================================================================
# Panel=sFigure1A - significant mutatuion sactter plot.
#=========================================================================
library ("readxl")
require("ggrepel")
require("ggplot2")
library("dplyr")

## ---------------  MutSig2CV SMGs   ---------------
MutSigCV_Mutsig <- read_excel("/Users/ranpeng/Desktop/Integrated_multiomics_analysis_CCA/data/Fig1/Significant_ mutation_gene.xlsx", sheet = 2)
colnames(MutSigCV_Mutsig)
MutSigCV_Mutsig = filter(MutSigCV_Mutsig,  n_nonsilent < 30 )
# genes$significant <- ifelse((genes$iCC_sig > 1.301| genes$eCC_sig > 1.301 ), "P_value < 0.05", "Not Sig")
MutSigCV_Mutsig$significant <-ifelse((MutSigCV_Mutsig$q < 0.1 & MutSigCV_Mutsig$n_nonsilent > 5 ), 'SMGs','Non_SMGs')
MutSigCV_Mutsig$Frequency = MutSigCV_Mutsig$n_nonsilent /138*100
MutSigCV_Mutsig$Mutation_count = MutSigCV_Mutsig$n_nonsilent
p = ggplot(MutSigCV_Mutsig, aes(x = Frequency, y = -log10(q), color = Frequency)) +
    geom_point(aes(color = Frequency,size = -log10(q))) +
    # scale_colour_gradient(low = "blue", high = "yellow") +
    scale_fill_gradient(low = "blue", high = "yellow") +
    # scale_color_manual(values = c("grey",'red',"royalblue", 'blue')) +
    theme_bw(base_size = 12) + theme(legend.position = "right") +
    geom_text_repel( 
      data = subset(MutSigCV_Mutsig, MutSigCV_Mutsig$q < 0.1 & MutSigCV_Mutsig$n_nonsilent >= 5  ),
      aes(label = gene),
      size = 6,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines")) + 
    geom_hline(yintercept = 1,color = 'red', linetype="dashed")+
    geom_vline(xintercept = 3.6,color = 'red', linetype="dashed")+
    xlim(0,25) #+
    # ylim(0,2.5) +
    # scale_y_continuous(limits = c(0, 3.5))
    # scale_x_continuous(breaks=seq(0, 30, 5))

p +
    ## Set title, labels
    # ggtitle("Significantly mutated genes")+
    labs(title="Significantly mutated genes", x="Mutation frequency (%)", y = "-Log10 (FDR)")+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.text.x=element_text(angle=0,hjust = 0.6,colour="black",family="Arial",size=16), 
          axis.title.x=element_text(angle=0,hjust = 0.5,colour="black",family="Arial",size=20),
          axis.text.y=element_text(family="Arial",size=16,face="plain"),
          axis.title.y=element_text(family="Arial",size = 20,face="plain"), 
          panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.8), 
          legend.text=element_text(face="plain", family="Arial", colour="black", size=12),
          legend.title=element_text(face="plain", family="Arial", colour="black", size=12),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
    # geom_hline(yintercept = 1.301,color = 'gray', linetype="dashed")+
    # geom_vline(xintercept = -0.57,color = 'gray', linetype="dashed") +
    # geom_vline(xintercept = 0.58,color = 'gray', linetype="dashed") 

## ---------------  oncodriveCLUST SMGs   ---------------
oncodriveCLUST_Mutsig <- read_excel("/Users/ranpeng/Desktop/Integrated_multiomics_analysis_CCA/data/Fig1/Significant_ mutation_gene.xlsx", sheet = 1)
colnames(oncodriveCLUST_Mutsig)
oncodriveCLUST_Mutsig = as.data.frame(oncodriveCLUST_Mutsig)
# oncodriveCLUST_Mutsig = filter(oncodriveCLUST_Mutsig,  n_nonsilent < 30 )
oncodriveCLUST_Mutsig$significant <-ifelse((oncodriveCLUST_Mutsig$fdr < 0.1 & oncodriveCLUST_Mutsig$total > 5 ), 'SMGs','Non_SMGs')
oncodriveCLUST_Mutsig$Frequency = oncodriveCLUST_Mutsig$total /138*100
oncodriveCLUST_Mutsig$Mutation_count = oncodriveCLUST_Mutsig$total
oncodriveCLUST_Mutsig$fdr = as.numeric(oncodriveCLUST_Mutsig$fdr)
p = ggplot(oncodriveCLUST_Mutsig, aes(x = Frequency, y = -log10(fdr), color = Frequency)) +
  geom_point(aes(color = Frequency,size = -log10(fdr))) +
  # scale_colour_gradient(low = "white", high = "red") +
  scale_fill_gradientn(colours = terrain.colors(7)) +
  # scale_color_manual(values = c("grey",'red',"royalblue", 'blue')) +
  theme_bw(base_size = 12) + theme(legend.position = "right") +
  geom_text_repel( 
    data = subset(oncodriveCLUST_Mutsig, oncodriveCLUST_Mutsig$fdr < 0.1 & oncodriveCLUST_Mutsig$total >= 5  ),
    aes(label = Hugo_Symbol),
    size = 6,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")) + 
  geom_hline(yintercept = 0.8,color = 'red', linetype="dashed")+
  geom_vline(xintercept = 3.6,color = 'red', linetype="dashed")+
  xlim(0,20) #+
# ylim(0,2.5) +
# scale_y_continuous(limits = c(0, 3.5))
# scale_x_continuous(breaks=seq(0, 30, 5))

p +
  ## Set title, labels
  # ggtitle("Significantly mutated genes")+
  labs(title="Significantly mutated genes", x="Mutation frequency (%)", y = "-Log10 (FDR)")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(angle=0,hjust = 0.6,colour="black",family="Arial",size=16), 
        axis.title.x=element_text(angle=0,hjust = 0.5,colour="black",family="Arial",size=20),
        axis.text.y=element_text(family="Arial",size=16,face="plain"),
        axis.title.y=element_text(family="Arial",size = 20,face="plain"), 
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.8), 
        legend.text=element_text(face="plain", family="Arial", colour="black", size=12),
        legend.title=element_text(face="plain", family="Arial", colour="black", size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


## --------------- SMGs scatter of the two combined data   ---------------
MutSigCV_Mutsig <- read_excel("/Users/ranpeng/Desktop/Integrated_multiomics_analysis_CCA/data/Fig1/SMGs.xlsx", sheet = 1)
colnames(MutSigCV_Mutsig)
MutSigCV_Mutsig = filter(MutSigCV_Mutsig,  n_nonsilent < 30 )
# genes$significant <- ifelse((genes$iCC_sig > 1.301| genes$eCC_sig > 1.301 ), "P_value < 0.05", "Not Sig")
MutSigCV_Mutsig$significant <-ifelse((MutSigCV_Mutsig$q < 0.1 & MutSigCV_Mutsig$n_nonsilent > 5 ), 'SMGs','Non_SMGs')
MutSigCV_Mutsig$Frequency = MutSigCV_Mutsig$n_nonsilent /138*100
MutSigCV_Mutsig$Mutation_count = MutSigCV_Mutsig$n_nonsilent
p = ggplot(MutSigCV_Mutsig, aes(x = Frequency, y = -log10(q), color = Frequency)) +
  geom_point(aes(color = Frequency,size = -log10(q), shape = type)) +
  # scale_colour_gradient(low = "blue", high = "yellow") +
  scale_fill_gradient(low = "blue", high = "yellow") +
  # scale_color_manual(values = c("grey",'red',"royalblue", 'blue')) +
  theme_bw(base_size = 12) + theme(legend.position = "right") +
  geom_text_repel( 
    data = subset(MutSigCV_Mutsig, MutSigCV_Mutsig$q < 0.1 & MutSigCV_Mutsig$n_nonsilent >= 5  ),
    aes(label = gene),
    size = 7,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")) + 
  geom_hline(yintercept = 1,color = 'red', linetype="dashed")+
  geom_vline(xintercept = 3.6,color = 'red', linetype="dashed")+
  xlim(0,25) #+
# ylim(0,2.5) +
# scale_y_continuous(limits = c(0, 3.5))
# scale_x_continuous(breaks=seq(0, 30, 5))

p +
  ## Set title, labels
  # ggtitle("Significantly mutated genes")+
  labs(title="Significantly mutated genes", x="Mutation frequency (%)", y = "-Log10 (FDR)")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(angle=0,hjust = 0.6,colour="black",family="Arial",size=16), 
        axis.title.x=element_text(angle=0,hjust = 0.5,colour="black",family="Arial",size=20),
        axis.text.y=element_text(family="Arial",size=16,face="plain"),
        axis.title.y=element_text(family="Arial",size = 20,face="plain"), 
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.8), 
        legend.text=element_text(face="plain", family="Arial", colour="black", size=12),
        legend.title=element_text(face="plain", family="Arial", colour="black", size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

## --------------- Panel= Fig1F - boxplot for mutation signature  -----------------------
library(ggplot2)
library(magrittr)
library(ggpubr)
library(cowplot)
library(tidyr)
library(reshape2)

## ---------- Violin plots represent the expression data ------------
# ## 数据类型为长格式
# Attribute  Group     value
# 1         Amp CDK11A  9.559318
# 2         Amp CDK11A  5.154437
# 3         Amp CDK11A  5.539539
# Data = read.csv("/Volumes/Samsung_T5/downloads/CCA/Data/Fig7/Neutrophils/Protein_boxplot20210605.csv")
Data = read.csv("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/Fig1/fig1f/AA_signature_mutation.csv")
head(Data)

melted_data<-melt(Data,id=c("CCA_No"), 
                  variable.name = "Group", 
                  value.name = "Attribute")
# write.csv(melted_data,"/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/Fig1/fig1f/AA_signature_mutation_1.csv")
## 过滤离群值
# Data = subset(Data,value <= 5)
## 运用盖帽法处理离群值
Data$value = ifelse(Data$value > 100, 100, Data$value)

P2<- ggplot(Data, aes(x=Group, y=log2(value+1),color=Attribute)) + 
  # geom_violin(trim=FALSE) + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  geom_boxplot(width=0.7,position=position_dodge(0.9),outlier.colour = "red")+ #绘制箱线图
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

## =================== AA mutation signature associated with immunescore ===============

# load boxplot_for_loop customed functions
source("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/R/boxplot_for_loop.R")
## ------------------------ input data -----------------------------
data <-read.csv("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/Fig1/fig1f/AA_signature_mutation.csv",header=T)
# data<-read.xlsx("/Users/ranpeng/Desktop/Desktop/项目文件/对外服务/罗俊一/clinical_info/clinical_info_arranged/Clinical_continuous_variables.xlsx",sheet=1,colNames = T)
colnames(data)
data$Sig4 = ifelse(data$Sig4 > 120, 120, data$Sig4)
## ----------------- 指定 category_var_colname 进行分析 --------------------
customed_boxplot(data, target = "TP53",category_var_name = "ImmuneScore",log2_transform = F, capping_outliers = T)


## ======================= Panel= Fig1E - RNA protein correlation  ===============

library(Hmisc)
library(dplyr)
## anatation
data1 = read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig1/RNA_protein_cor/RNA_1.csv",header = T,row.names = 1)
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
write.csv(result,"/Users/ranpeng/Desktop/CCA/Data/Fig1/RNA_protein_cor/results/Tumor/RNA_Protein_cor_spearman.csv")
result$cor_r = as.numeric(result$cor_r);result$p_val = as.numeric(result$p_val)
result_filter <- filter(result, p_val < 0.05 & abs(cor_r) > 0.2)
(result_pos_cor <- filter(result, p_val < 0.05 & cor_r > 0.2))
(result_neg_cor <- filter(result, p_val < 0.05 & cor_r < -0.2))
write.csv(result_neg_cor,"/Users/ranpeng/Desktop/CCA/Data/Fig7/Xcell_signature_cor/Th1_cell_cor_molecular_r0.2_neg.csv")


## ======================= GSEA analysis ===============

library(DESeq2)
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(pheatmap)

## geneID 转换
data = read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig7/Xcell_signature_cor/Neutrophils/")
names(data)[1] = "SYMBOL"
rownames(data) = data[,1]
geneid <- rownames(data)
genes  <- clusterProfiler::bitr(geneid, fromType = 'SYMBOL', toType = c('ENSEMBL','ENTREZID'), OrgDb = 'org.Hs.eg.db')


genes <- genes[!duplicated(genes$SYMBOL),]
genes <- genes[!duplicated(genes$SYMBOL),] %>% dplyr::inner_join(data, 'SYMBOL')
# (1) 建立输入对象并分析表达谱在GO BP和KEGG中的富集程度。
# Create input objects and analyze the degree of enrichment of expression profiles in GO BP and KEGG.
input_GSEA <- genes$cor_r
names(input_GSEA) <- genes$ENTREZID
input_GSEA <- sort(input_GSEA, decreasing = T)
GSEGO_BP <- gseGO(input_GSEA, ont = 'BP', OrgDb = org.Hs.eg.db, nPerm = 1000, pvalueCutoff = 0.05)
GSEGO_BP <- setReadable(GSEGO_BP, org.Hs.eg.db, keyType = "ENTREZID")
GSEA_KEGG <- gseKEGG(input_GSEA, organism = 'hsa', keyType = 'ncbi-geneid', nPerm = 1000, pvalueCutoff = 0.05)
GSEA_KEGG <- setReadable(GSEA_KEGG, org.Hs.eg.db, keyType = "ENTREZID")
# (2) GSEA结果的可视化和解读。
# Visualization and interpretation of GSEA results.

#further filter out the significant gene sets and order them by NES scores.
# GO pathway enrichment
GSEA_BP_df <- as.data.frame(GSEGO_BP) %>% dplyr::filter(abs(NES)>1 & pvalue<0.05 & qvalues<0.25)
GSEA_BP_df <- GSEA_BP_df[order(GSEA_BP_df$NES, decreasing = T),]
p <- gseaplot(GSEGO_BP, GSEA_BP_df$ID[1], by='all', title = GSEA_BP_df$Description[1], color.vline = 'gray50', color.line='red', color='black') #by='runningScore/preranked/all'
p <- p+annotate(geom = 'text', x = 0.87, y =0.85, color='red', fontface = 'bold', size= 4,
                label=paste0('NES= ', round(GSEA_BP_df$NES[1], 1), '\n', 'p.adj= ', round(GSEA_BP_df$p.adjust[1], 2)))+
  theme(panel.grid = element_line(colour = 'white'))
p
# KEGG pathway enrichment
GSEA_KEGG_df <- as.data.frame(GSEA_KEGG) %>% dplyr::filter(abs(NES)>1 & pvalue<0.05 & qvalues< 0.25)
GSEA_KEGG_df <- GSEA_KEGG_df[order(GSEA_KEGG_df$NES, decreasing = T),]
i= 1
p <- enrichplot::gseaplot2(GSEA_KEGG, geneSetID = GSEA_KEGG_df$ID[i], pvalue_table = F, ES_geom = 'line')+
  annotate('text', x = 0.87, y =0.85, color='red', fontface = 'bold', size= 4,
           label=paste0('NES= ', round(GSEA_KEGG_df$NES[i], 1), '\n', 'p.adj= ', round(GSEA_KEGG_df$p.adjust[1], 2)))+
  labs(title = GSEA_KEGG_df$Description[i])+theme(plot.title = element_text(hjust = 0.5))
p

GSEA_KEGG$Description
# enrichplot::gseaplot2(GSEA_KEGG, c(1,5,6,8))
enrichplot::gseaplot2(GSEA_KEGG, geneSetID = c(9,17,12), color = c('red','#56AF61','#3C679E')) #
# write.csv(GSEA_KEGG_df,"/Users/ranpeng/Desktop/CCA/Data/Fig7/Xcell_signature_cor/Neutrophils/GSEA/Neutrophils_GSEA.csv")

# (3) 自建基因集合进行GSEA分析。
# Self-constructed gene collections for GSEA analysis.
db_for_GSEA = read.csv('/Users/ranpeng/Desktop/CCA/CCA_script/tf_activities/High_CI_regulon_dorothea.csv',header = F, fill = TRUE)
# ------ ----------- Create gmt standard format pathway database ------------------
# original_gmt <- readLines('/Users/ranpeng/Desktop/CCA/CCA_script/tf_activities/High_CI_regulon_dorothea.gmt')
# strsplit_name <- function(gmt.list_layer){
#   GSgenes <- as.character(unlist(strsplit(gmt.list_layer, split = '\t', fixed = T)))
#   data.frame('Genesets'=rep(GSgenes[1], (length(GSgenes)-2)), 'Genes'=GSgenes[c(-1,-2)], stringsAsFactors = F)
# }
# database_list <- lapply(original_gmt, strsplit_name)
# db_for_GSEA <- do.call(rbind, database_list)
# db_for_GSEA$Genesets <- as.character(db_for_GSEA$Genesets)
# db_for_GSEA$Genes <- as.character(db_for_GSEA$Genes)
colnames(db_for_GSEA) = c("Genesets","Genes")
input_GSEA <- genes$logFC
names(input_GSEA) <- genes$SYMBOL
input_GSEA <- sort(input_GSEA, decreasing = T)
GSEA_test <- clusterProfiler::GSEA(input_GSEA, nPerm = 1000, pvalueCutoff = 0.5, TERM2GENE = db_for_GSEA)
GSEA_test_df <- as.data.frame(GSEA_test) %>% dplyr::filter(abs(NES)>1 & pvalue<0.05 & qvalues<0.25)
GSEA_test_df <- GSEA_test_df[order(GSEA_test_df$NES, decreasing = T),]
enrichplot::gseaplot2(GSEA_test, GSEA_test_df$Description[c(1,4, (length(GSEA_test_df$Description)-2),length(GSEA_test_df$Description))], color = c('red','orange','blue','navy'))
write.csv(GSEA_test_df,"/Users/ranpeng/Desktop/CCA/Data/Fig3/TF_active/TF_GSEA_50%_new.csv")
enrichplot::gseaplot2(GSEA_test, GSEA_test_df$Description[c((length(GSEA_test_df$Description)-7),length(GSEA_test_df$Description))], color = c('blue','navy'))


i= 2
p <- enrichplot::gseaplot2(GSEA_test, geneSetID = GSEA_test_df$ID[i], pvalue_table = F, ES_geom = 'line')+
  annotate('text', x = 0.87, y =0.85, color='red', fontface = 'bold', size= 4,
           label=paste0('NES= ', round(GSEA_test_df$NES[i], 1), '\n', 'p.adj= ', round(GSEA_test_df$p.adjust[1], 2)))+
  labs(title = GSEA_test_df$Description[i])+theme(plot.title = element_text(hjust = 0.5))
p

## ======================= gene mutation fisher test  ===============

library('openxlsx')
mydat<-read.xlsx("/Users/ranpeng/Desktop/CCA/Data/Fig1/Figure-1d/CCA_other_cohot_gene_mutation_info.xlsx",sheet=9,colNames = T)

result<-c()
for (i in (2:ncol(mydat))){
  newdat<-mydat[,i]
  test<-matrix(newdat,nrow=5,ncol=2,byrow=TRUE)
  qq<-fisher.test(test,simulate.p.value=TRUE)   #fisher.test()   chisq.test
  Odds_ratio<-qq[["estimate"]][["odds ratio"]]
  p.result<-qq$p.value
  print(p.result)
  result.linshi<-cbind(i,Odds_ratio,p.result)
  result<-rbind(result,result.linshi)
}
gene<-names(mydat[,-1])  ##提取列名，对应的row.names(mydat)
result<-cbind(gene,result)
# result$p.adjust = result$p.result /14
write.csv(result,"/Users/ranpeng/Desktop/CCA/Data/Fig1/Figure-1d/gene_mutaion_fisher_result.csv")

result = read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig1/Figure-1d/gene_mutaion_fisher_result.csv")
result$p.adjust = p.adjust(result$p_val, method = "bonferroni", n = length(result$p_val))
