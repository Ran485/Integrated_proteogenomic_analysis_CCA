#=========================================================================
# response R2-C3 - FGFR2 fusion
#=========================================================================
source("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/R/stacked_barplot.R")

## ------------------------ iCCA and eCCA FGFR2 fusion -----------------------------
data = read.csv('/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/FGFR2_data/FGFR_LS_BD.csv')
colnames(data)
## ----------------- 指定 category_var_colname 进行分析 --------------------
# data = data[data$iCC_eCC != 'eCCA',]
# data = data[data$iCCA_Pathology != 'NA_1',]
stacked_barplot(data, target = 'FGFR2_Alter',category_var_colname = 'iCCA_Pathology')


## ------------------------ FGFR2 fusion cohort compared in iCCA -----------------------------
data = read.csv('/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/FGFR2_cohort_compared.csv')
colnames(data)
## ----------------- 指定 category_var_colname 进行分析 --------------------
stacked_barplot(data, target = 'cohort',category_var_colname = 'FGFR2_fusion')

## ------------------------Frequency of FGFR2 partner -----------------------------
library ("readxl")
FGFR2_Fusions = read_excel("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/Fig1/FGFR2_fusions.xlsx", sheet = 6)
FGFR2_Fusions_Deng = FGFR2_Fusions[FGFR2_Fusions$cohort=='Deng_cohort',]
FGFR2_Fusions_Dong = FGFR2_Fusions[FGFR2_Fusions$cohort=='Dong_cohort',]
FGFR2_Fusions_Deng = as.data.frame(sort(table(FGFR2_Fusions_Deng$RightGene), decreasing = T))
FGFR2_Fusions_Dong = as.data.frame(sort(table(FGFR2_Fusions_Dong$RightGene), decreasing = T))
merge_Data = rbind(FGFR2_Fusions_Deng,FGFR2_Fusions_Dong)
write.csv(merge_Data,"/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/FGFR2_cohort_count.csv")

df = read.csv("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/FGFR2_Chormosome_count.csv", row.names = 1)
colnames(df)
# df$count = ifelse(df$count >= 800, 800, df$count)
# Create the barplot
# factor = as.data.frame(sort(table(FGFR2_Fusions$Chormosome), decreasing = T))
# write.csv(factor,"/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/FGFR2_Chormosome_count.csv")

df$Chr = factor(df$Chr, levels = df$Chr)
# ggplot(data=df, aes(x=Chr, y=num, fill=cohort)) +
ggplot(data=df, aes(x=Chr, y=num)) +
  geom_bar(stat="identity")+
  # geom_text(aes(y=label_ypos, label=len), vjust=1.6, 
  #           color="white", size=3.5)+
  scale_fill_brewer(palette="Paired")+
  # xlim(0,25) +
  ylim(0,22) +
  labs(title=" Distribution of FGFR2 rearrangement partners in the genome", x="Chromosome position", y = "Number of patients")+
  # theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  #       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(angle=45,hjust = 0.6,colour="black",family="Arial",size=14), 
        axis.title.x=element_text(angle=0,hjust = 0.5,colour="black",family="Arial",size=20),
        axis.text.y=element_text(family="Arial",size=16,face="plain"),
        axis.title.y=element_text(family="Arial",size = 20,face="plain"), 
        # panel.border = element_blank(),
        axis.line = element_line(colour = "black",size=0.8),
        legend.text=element_text(face="plain", family="Arial", colour="black", size=12),
        legend.title=element_text(face="plain", family="Arial", colour="black", size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

## ------------------------ Fisher exact test for FGFR2 fusions  -----------------------------
source("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/R/fisher_exact_test.R")

# data = read.csv('/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/Fig6/LS_arm_events.csv')
data = read.csv('/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/FGFR2_mutation_matrix.csv')
colnames(data)
features = colnames(data)[3:ncol(data)]
significant_res = fisher_exact_test(data, target='FGFR2_fusion', features = features)


## --------------- oncoprint_mutation for FGFR2 alterations ----------------------

library(ComplexHeatmap)
library(RColorBrewer)

onco_mat = read.csv('/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/Fig1/fig1b/onco_matrix_BAP1(2022_03_26).csv',header = T, row.names = 1)
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


## ----------------------- FGFR2 boxplot for RNA and Protein
source("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/boxplot_for_loop.R")
library('readxl')
## ------------------------ input data -----------------------------
# data <-read.csv("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/Fig1/fig1f/AA_signature_mutation.csv",header=T)
data<-read.xlsx("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/FGFR2_data/FGFR2_expression.xlsx",sheet=1,colNames = T)
colnames(data)
data$FGFR2_protein = ifelse(data$FGFR2_protein > 1, 1, data$FGFR2_protein)
## ----------------- 指定 category_var_colname 进行分析 --------------------
customed_boxplot(data, target = "type",category_var_name = "FGFR2_protein",
                 log2_transform = T, capping_outliers = F, violin_plot = T)

# category_list = c('Basophil_count', 'HBP', 'Plasma_fibrinogen', 'Chlorine')
category_list = colnames(data)[-c(1:2)]
category_var_plots = list()

for (category_var in category_list){
  category_var_plots[[category_var]] = customed_boxplot(data, target = 'type',category_var_name = category_var, log2_transform = T, violin_plot = T)
  print(category_var_plots[[category_var]])
  # ggsave(category_var_plots[[category_var]], file=paste0("plot_", category_var,".png"), width = 44.45, height = 27.78, units = "cm", dpi=300)
}

## 生成多个图关键的一步
plot_grid(plotlist = category_var_plots)


## ----------------------- FGFR2 fusions and mutations stacked plot
source("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/R/stacked_barplot.R")

data<-read.xlsx("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/FGFR2_data/FGFR2_mutation.xlsx",sheet=1,colNames = T)
# data = read.csv('/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/FGFR2_data/FGFR_LS_BD.csv')
colnames(data)
## ----------------- 指定 category_var_colname 进行分析 --------------------
# data = data[data$iCC_eCC != 'eCCA',]
# data = data[data$iCCA_Pathology != 'NA_1',]
stacked_barplot(data, target = 'iCC_eCC',category_var_colname = 'FGFR2_Mut')

## ----------------------- FGFR2 fusions and mutations boxplot for RNA and Protein
source("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/R/boxplot_for_loop.R")
library('readxl')
## ------------------------ input data -----------------------------
# data <-read.csv("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/Fig1/fig1f/AA_signature_mutation.csv",header=T)
data<-read.xlsx("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/FGFR2_data/FGFR2_fusion_mutation_expression.xlsx",sheet=1,colNames = T)
colnames(data)
data = data[data$FGFR2_Mut != "NA_1",]
# data$FGFR2_protein = ifelse(data$FGFR2_protein > 1, 1, data$FGFR2_protein)
## ----------------- 指定 category_var_colname 进行分析 --------------------
customed_boxplot(data, target = "type",category_var_name = "FGFR2_protein",
                 log2_transform = T, capping_outliers = F, violin_plot = T)

## ---------------------- pathway barplot --------------------------------
# df = read.csv("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/FGFR2_Chormosome_count.csv", row.names = 1)
df = read_excel("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/FGFR2_data/pathway_enrichment.xlsx",sheet = 2)
colnames(df)
# df$count = ifelse(df$count >= 800, 800, df$count)
# Create the barplot
# factor = as.data.frame(sort(table(FGFR2_Fusions$Chormosome), decreasing = T))
# write.csv(factor,"/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/FGFR2_Chormosome_count.csv")

df$pathway = factor(df$pathway, levels = df$pathway)
# ggplot(data=df, aes(x=Chr, y=num, fill=cohort)) +
ggplot(data=df, aes(x=pathway, y=-log10(`p-value`))) +
  geom_bar(stat="identity")+
  # geom_text(aes(y=label_ypos, label=len), vjust=1.6, 
  #           color="white", size=3.5)+
  scale_fill_brewer(palette="Paired")+
  # xlim(0,25) +
  # ylim(0,22) +
  labs(title="Pathway enrichment in FGFR2 alteration vs. WT", x="Chromosome position", y = "-Log10(p-value)")+
  # theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  #       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(angle=0,hjust = 0.6,colour="black",family="Arial",size=14), 
        axis.title.x=element_text(angle=0,hjust = 0.5,colour="black",family="Arial",size=20),
        axis.text.y=element_text(family="Arial",size=16,face="plain"),
        axis.title.y=element_text(family="Arial",size = 20,face="plain"), 
        # panel.border = element_blank(),
        axis.line = element_line(colour = "black",size=0.8),
        legend.text=element_text(face="plain", family="Arial", colour="black", size=12),
        legend.title=element_text(face="plain", family="Arial", colour="black", size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + coord_flip()

#=========================================================================
# response R2-C4 - LS BD SMAD4 gene mutation
#=========================================================================
## ----------------------- SMAD4 gene mutation correlation --------------------
source("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/R/fisher_exact_test.R")
data = read.csv('/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/Fig6/LS_duct_mutation_events.csv')
colnames(data)
features = colnames(data)[4:ncol(data)]
significant_res = fisher_exact_test(data, target='iCCA_Pathology', features = features)
# write.csv(significant_res,"/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/Fig6/LS_duct_mutation_events_fishertest.csv")

## ----------------------- SMAD4 boxplot for RNA and Protein --------------------
source("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/R/boxplot_for_loop.R")
library('readxl')
## ------------------------ input data -----------------------------
data <-read.csv("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/SMAD4/SMAD4_CA199.csv",header=T)
# data<-read.xlsx("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/response_R2C4/GaoQ_SMAD4_expression.xlsx",sheet=1,colNames = T)
colnames(data)
# data$FGFR2_protein = ifelse(data$FGFR2_protein > 1, 1, data$FGFR2_protein)

## ----------------- 指定 category_var_colname 进行分析 --------------------
customed_boxplot(data, target = "SMAD4",category_var_name = "Preoperative.CA19.9_U.ml_",
                 log2_transform = T, capping_outliers = F, violin_plot = T)


## ---------------------- SAMD4 pathway barplot --------------------------------
# df = read.csv("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/FGFR2_Chormosome_count.csv", row.names = 1)
if (TRUE){
df = read_excel("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/response_R2C4/Pathway_enrichment.xlsx",sheet = 10)
colnames(df)
# df$count = ifelse(df$count >= 800, 800, df$count)
# Create the barplot
# factor = as.data.frame(sort(table(FGFR2_Fusions$Chormosome), decreasing = T))
# write.csv(factor,"/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/FGFR2_Chormosome_count.csv")

df$pathway = factor(df$pathway, levels = df$pathway)
# ggplot(data=df, aes(x=Chr, y=num, fill=cohort)) +
ggplot(data=df, aes(x=pathway, y=-log10(`p-value`))) +
  geom_bar(stat="identity")+
  # geom_text(aes(y=label_ypos, label=len), vjust=1.6, 
  #           color="white", size=3.5)+
  scale_fill_brewer(palette="Paired")+
  # xlim(0,25) +
  # ylim(0,22) +
  labs(title="Pathway enrichment in FSAMD4 alteration vs. WT", x="Chromosome position", y = "-Log10(p-value)")+
  # theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  #       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(angle=0,hjust = 0.6,colour="black",family="Arial",size=14), 
        axis.title.x=element_text(angle=0,hjust = 0.5,colour="black",family="Arial",size=20),
        axis.text.y=element_text(family="Arial",size=16,face="plain"),
        axis.title.y=element_text(family="Arial",size = 20,face="plain"), 
        # panel.border = element_blank(),
        axis.line = element_line(colour = "black",size=0.8),
        legend.text=element_text(face="plain", family="Arial", colour="black", size=12),
        legend.title=element_text(face="plain", family="Arial", colour="black", size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + coord_flip()
}

## ---------------------- 6q14.3 deletion SH3BGRL2 xpression --------------------------------
source("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/R/boxplot_for_loop.R")
library('readxl')
## ------------------------ input data -----------------------------
# data = read.csv('/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/Fig1/fig1f/AA_signature_mutation_2.csv')
data<-read.xlsx("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/SMAD4/large_small.xlsx",sheet=1,colNames = T)
colnames(data)
# data$Sig4 = ifelse(data$Sig4 > 160, 160, data$Sig4)
## ----------------- 指定 category_var_colname 进行分析 --------------------
# data = data[data$Tbil!='NA_1',]
# data = data[data$iCC_eCC =='eCCA',]
customed_boxplot(data, target = "iCCA_Pathology",category_var_name = "SH3BGRL2_Pro",
                 log2_transform = T, capping_outliers = T, violin_plot = T)



#=========================================================================
# response R1-C3 - ARID1A gene mutation
#=========================================================================

## ------------------------ BAP1 and ARID1A gene mutation -----------------------------
data = read_excel('/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/response_R1C3_ARID1A/ARID1A_index.xlsx',sheet = 1)
colnames(data)
## ----------------- 指定 category_var_colname 进行分析 --------------------
# data = data[data$iCC_eCC != 'eCCA',]
# data = data[data$iCCA_Pathology != 'NA_1',]
stacked_barplot(data, target = 'ARID1A',category_var_colname = 'BAP1')

fisher.test(table(data$BAP1,data$ARID1A))

## ----------------------- BAP1 and ARID1A boxplot for RNA and Protein
source("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/boxplot_for_loop.R")
library('readxl')
## ------------------------ input data -----------------------------
data <-read.csv("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/response_R1C3_ARID1A/BAP1/BAP1_expression.csv",header=T)
# data<-read.xlsx("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/FGFR2_data/FGFR2_expression.xlsx",sheet=1,colNames = T)
colnames(data)
# data$FGFR2_protein = ifelse(data$FGFR2_protein > 1, 1, data$FGFR2_protein)
## ----------------- 指定 category_var_colname 进行分析 --------------------
customed_boxplot(data, target = "type",category_var_name = "BAP1_Pro",
                 log2_transform = T, capping_outliers = F, violin_plot = T)

## ------------------------ BAP1 and ARID1A distribulation -----------------------------
source("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/R/stacked_barplot.R")
# data = read.csv('/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/FGFR2_data/FGFR_LS_BD.csv')
data<-read.xlsx("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/response_R1C3_ARID1A/ARID1A_index.xlsx",sheet=1,colNames = T)
colnames(data)
## ----------------- 指定 category_var_colname 进行分析 --------------------
# data = data[data$iCC_eCC != 'eCCA',]
# data = data[data$iCCA_Pathology != 'NA_1',]
stacked_barplot(data, target = 'Pathology',category_var_colname = 'ARID1A')

## ---------------------- ARID1A pathway barplot --------------------------------
# df = read.csv("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/FGFR2_Chormosome_count.csv", row.names = 1)
df = read_excel("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/response_R1C3_ARID1A/pathway_enrichment.xlsx",sheet = 10)
colnames(df)
# df$count = ifelse(df$count >= 800, 800, df$count)
# Create the barplot
# factor = as.data.frame(sort(table(FGFR2_Fusions$Chormosome), decreasing = T))
# write.csv(factor,"/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/FGFR2_Chormosome_count.csv")

df$pathway = factor(df$pathway, levels = df$pathway)
# ggplot(data=df, aes(x=Chr, y=num, fill=cohort)) +
ggplot(data=df, aes(x=pathway, y=-log10(`p-value`))) +
  geom_bar(stat="identity")+
  # geom_text(aes(y=label_ypos, label=len), vjust=1.6, 
  #           color="white", size=3.5)+
  scale_fill_brewer(palette="Paired")+
  # xlim(-8,5) +
  ylim(0,3) +
  labs(title="Pathway enrichment in FGFR2 alteration vs. WT", x="Chromosome position", y = "-Log10(p-value)")+
  # theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  #       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(angle=0,hjust = 0.6,colour="black",family="Arial",size=14), 
        axis.title.x=element_text(angle=0,hjust = 0.5,colour="black",family="Arial",size=20),
        axis.text.y=element_text(family="Arial",size=16,face="plain"),
        axis.title.y=element_text(family="Arial",size = 20,face="plain"), 
        # panel.border = element_blank(),
        axis.line = element_line(colour = "black",size=0.8),
        legend.text=element_text(face="plain", family="Arial", colour="black", size=12),
        legend.title=element_text(face="plain", family="Arial", colour="black", size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + coord_flip()

## ---------------------- ARID1A pathway moleculor vocanoplot --------------------------------
library(ggplot2)
library(hrbrthemes)
#install.packages("hrbrthemes")
library(ggpubr)
# library(ggpmisc)
library(ggrepel)

matrix <- read_excel('/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/response_R1C3_ARID1A/pathway_enrichment_volcano_protein.xlsx', col_names = T)
colnames(matrix)
matrix$Significant <- ifelse(matrix$display == 'displar', "display", "not_display")
matrix$size <- ifelse(matrix$size == 2, 2, 1)
# matrix$size <- ifelse(matrix$display == 'display', 2, 1)
ggplot(data = matrix, aes(x = logFC, y = -log10(P.Value), color=pathway, size = size)) +
  geom_point() + ## aes(color = Significant)
  scale_size_continuous(range = c(1,2))+
  scale_color_brewer(palette="Paired") +
  # scale_color_manual(values = c("#515151", "red","#7AA9CE")) + # ("gray", "#7AA9CE") c( "#EA686B","#7AA9CE","gray"))
  theme_bw(base_size = 12) + theme(legend.position = "right") +
  geom_text_repel(
    data = subset(matrix, Significant == "display"),
    aes(label = gene),
    size = 5,
    max.overlaps = Inf,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines",
    )
  )+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle("FGFR2 alteration on RNA and Protein")+
  theme(axis.text.x=element_text(angle=0,hjust = 0.6,colour="black",family="Arial",size=16), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.title.x=element_text(angle=0,hjust = 0.5,colour="black",family="Arial",size=16),
        axis.text.y=element_text(family="Arial",size=16,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Arial",size = 16,face="plain"), #设置y轴标题的字体属性
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.8), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="plain", family="Arial", colour="black",  #设置图例的子标题的字体属性
                                 size=12),
        legend.title=element_text(face="plain", family="Arial", colour="black", #设置图例的总标题的字体属性
                                  size=12)) +
  # panel.grid.major = element_blank(),   #不显示网格线
  # panel.grid.minor = element_blank())+  #不显示网格线
  geom_hline(yintercept = 2,color = 'gray', linetype="dashed") +
  # geom_hline(yintercept = 0.35,color = 'gray', linetype="dashed") +
  geom_vline(xintercept = -0.65,color = 'gray', linetype="dashed") +
  geom_vline(xintercept = 0.65,color = 'gray', linetype="dashed") +
  # scale_x_continuous(limits = c(-1.8, 1.8))+
  # scale_y_continuous(limits = c(-1.8, 1.8))
  # scale_x_continuous(breaks=seq(-2, 2, 0.5)) +
  # scale_y_continuous(breaks=seq(-2.5, 2.5, 1)) +
  ylab("Log10 [p-value]")+xlab("Log2FC (ARID1A/others)-protein")  #设置x轴和y轴的标题
  # ylab("Log2FC (FGFR2/others)-RNA")+xlab("Log2FC (FGFR2/others)-protein")  #设置x轴和y轴的标题
#=========================================================================
# response R2-C4 - gene mutation correlated with large and small BD
#=========================================================================
source("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/R/fisher_exact_test.R")
data = read.csv('/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/Fig6/LS_duct_mutation_events.csv')
colnames(data)
features = colnames(data)[4:ncol(data)]
significant_res = fisher_exact_test(data, target='iCCA_Pathology', features = features)

# write.csv(significant_res,"/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/Fig6/LS_duct_mutation_events_fishertest.csv")

#=========================================================================
# response R1-C1 - Clinical Baseline Tables cohort
#=========================================================================
library(readxl)
data = read_excel("/Volumes/Samsung_T5/downloads/CCA/Revised_data/revised_manuscript/Clinical Baseline Tables.xlsx",sheet=10, col_names = T)
colnames(data)
# data = as.matrix(data)
data = data[,-1]
data = t(data)

fisher.test(data) # simulate.p.value=TRUE

chisq.test(data, rescale.p = T)

#=========================================================================
# response R1-C1 - AA clinical relevance
#=========================================================================
source("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/R/stacked_barplot.R")

## ------------------------ iCCA and eCCA FGFR2 fusion -----------------------------
data = read.csv('/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/Fig1/fig1f/AA_signature_mutation.csv')
colnames(data)
## ----------------- 指定 category_var_colname 进行分析 --------------------
# data = data[data$iCC_eCC != 'eCCA',]
# data = data[data$iCCA_Pathology != 'NA_1',]
# stacked_barplot(data, target = 'FGFR2_Alter',category_var_colname = 'iCCA_Pathology')

## --------- for loop 对指定的 category_var_colname list 批量进行分析 -----------
# category_list = c("Heart_failure","Gender","Ethnicity","Comorbidities_modified","Smoking","Family_History","Lateral_branch_blood_supply","Processing","ACEIorARB","CCB","ARNI","Nitrates")
# category_list = significant_res$test_features
category_list = c('gender', 'iCC_eCC',
                  'PT', 'HBsAg', 'HBcAb', 'HBV', 'GGT', 'ALT', 'Tbil', 'CA199',
                  'Cholelithiasis_1', 'Type_2_diabetes', 'Hypertension',
                  'Presence_of_fluke_infection', 'Perineural_invasion',
                  'Lymphovascular_invasion')
category_var_plots = list()
for(category_var in category_list){
  category_var_plots[[category_var]] = stacked_barplot(data, target = 'TP53',category_var_colname = category_var, binary = F)
  print(category_var_plots[[category_var]])
  # ggsave(category_var_plots[[category_var]], file=paste0("plot_", category_var,".png"), width = 44.45, height = 27.78, units = "cm", dpi=300)
}

## 生成多个图关键的一步
plot_grid(plotlist = category_var_plots)


## ----------------------- AA expotrue boxplot for RNA and Protein
source("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/R/boxplot_for_loop.R")
library('readxl')
## ------------------------ input data -----------------------------
data = read.csv('/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/Fig1/fig1f/AA_signature_mutation_2.csv')
# data<-read.xlsx("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/AA_signature/AA_signature.xlsx",sheet=1,colNames = T)
colnames(data)
data$Sig4 = ifelse(data$Sig4 > 160, 160, data$Sig4)
## ----------------- 指定 category_var_colname 进行分析 --------------------
data = data[data$Tbil!='NA_1',]
data = data[data$iCC_eCC =='eCCA',]
customed_boxplot(data, target = "Tbil",category_var_name = "Sig4",
                 log2_transform = T, capping_outliers = F, violin_plot = T)

# iCC_eCC_sig4_Exposure iCC_eCC_sig4_Exposure_violin
## --------- for loop 对指定的 category_var_colname list 批量进行分析 -----------
colnames(data)
# category_list = c("Heart_failure","Gender","Ethnicity","Comorbidities_modified","Smoking","Family_History","Lateral_branch_blood_supply","Processing","ACEIorARB","CCB","ARNI","Nitrates")
category_list = colnames(data)[-c(1:13)]
# category_list = c('Basophil_count', 'HBP', 'Plasma_fibrinogen', 'Chlorine')
category_var_plots = list()

for (category_var in category_list){
  category_var_plots[[category_var]] = customed_boxplot(data, target = category_var ,category_var_name = 'Sig4', log2_transform = F)
  print(category_var_plots[[category_var]])
  # ggsave(category_var_plots[[category_var]], file=paste0("plot_", category_var,".png"), width = 44.45, height = 27.78, units = "cm", dpi=300)
}

## 生成多个图关键的一步
plot_grid(plotlist = category_var_plots)

## ----------------------- DOCR1 boxplot for RNA and Protein
source("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/R/boxplot_for_loop.R")
library('readxl')
## ------------------------ input data -----------------------------
# data = read.csv('/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/Fig1/fig1f/AA_signature_mutation_2.csv')
data<-read.xlsx("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/response_R2C2/DPCR1_mut_expre.xlsx",sheet=1,colNames = T)
colnames(data)
## ----------------- 指定 category_var_colname 进行分析 --------------------
customed_boxplot(data, target = "DPCR1_mut",category_var_name = "XIC_FOT",
                 log2_transform = F, capping_outliers = F, violin_plot = F)

