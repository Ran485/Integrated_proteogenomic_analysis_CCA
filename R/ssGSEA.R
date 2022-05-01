#source("https://bioconductor.org/biocLite.R")
#biocLite("GSVA")
library(GSVA)
library(piano)
library(limma)
library(BBmisc)
library(ComplexHeatmap)
library(circlize)
library(Hmisc)
library(GSEABase)
##------load data set ----------------------------------------------------------
# load('/Users/ranpeng/OneDrive - mail.ecust.edu.cn/Desktop-备份/科研项目/胆管癌-CCA/data/2020-08/capsule-0772343/GSEA in cancer signaling/data/Datasets/MatchedNormal-Metabric.RData')    
# load('/Users/ranpeng/OneDrive - mail.ecust.edu.cn/Desktop-备份/科研项目/胆管癌-CCA/data/2020-08/capsule-0772343/GSEA in cancer signaling/data/Datasets/MatchedTumor-Metabric.RData')    
# # path = ''

##-----create one matrix -----
# colnames(edata.normalmatched) <- paste(colnames(edata.normalmatched), "tumor", sep=".")
# colnames(edata.tumormatched) <- paste(colnames(edata.tumormatched), "tumor", sep=".")
# 
# genes <- intersect(rownames(edata.normalmatched), rownames(edata.tumormatched))
# df <- cbind(edata.normalmatched[genes,], edata.tumormatched[genes,])

##-----------------------load gmt file -----------------------------------------

KEGG <- getGmt('/Users/ranpeng/Desktop/Xcell_score/1/GSVA/geneset/c2.cp.kegg.v7.2.symbols.gmt')
Hall_mark <- getGmt('/Users/ranpeng/Desktop/Xcell_score/1/GSVA/geneset/h.all.v7.2.symbols.gmt')
Reactome <- getGmt('/Users/ranpeng/Desktop/Xcell_score/1/GSVA/geneset/c2.cp.reactome.v7.2.symbols.gmt')
GO_BP <- getGmt('/Users/ranpeng/Desktop/Xcell_score/1/GSVA/geneset/c5.go.bp.v7.2.symbols.gmt')
GO_ALL <- getGmt('/Users/ranpeng/Desktop/Xcell_score/1/GSVA/geneset/c5.all.v7.2.symbols.gmt')

ssg_res <- gsva(matrix, GO_ALL, min.sz= 10, max.sz= 300,
                method="gsva",
                verbose=TRUE, parallel.sz = 1)

write.csv(ssg_res,"/Users/ranpeng/Desktop/CCA/Data/Fig2/Arm/6q14.3/GSVA/result/GO_ALL_GSVA.csv")
## ssGSEA Fuction
ssGSEA <- function(matrix, geneset){
  
  ssg_res <- gsva(matrix, KEGG, min.sz= 10, max.sz= 300, 
                  verbose=TRUE, parallel.sz = 18)
  
  ##----- run limma to find differential pathway ---------------------------------
  design <- model.matrix(~ factor(meta$type))
  colnames(design) <- c("iCCA", "iCCAVseCCA")
  
  fit <- eBayes(lmFit(ssg_res, design))
  diffPath <- topTable(fit, coef = "iCCAVseCCA", number=Inf)
  
  ##------to plot we take top 20 pathways by logFC -------------------------------
  
  diffPath <- sortByCol(diffPath, c("logFC"))
  topPathways <- c(rev(rownames(diffPath))[1:20],rownames(diffPath)[1:20])
  mat <- ssg_res[topPathways, ]
  
  enrichment_pathway = data.frame(diffPath)
  top20_pathway = data.frame(mat)
  # Print results objects to workspace 
  
  write.csv(ssg_res,paste(getwd(), geneset, paste('ssGSEA_results',".csv",sep = ''),sep = '/'))
  write.csv(enrichment_pathway, paste(getwd(), geneset,paste('ssGSEA_results_pvalue',".csv",sep = ''),sep = '/'))
  write.csv(top20_pathway,paste(getwd(), geneset, paste('top20_pathway',".csv",sep = ''),sep = '/'))
  
  # write.csv(enrichment_pathway,"/Users/ranpeng/Desktop/CCA-data/2021-1-18/HBV/ssGSEA_results_pvalue.csv")
  # write.csv(enrichment_pathway,"/Users/ranpeng/Desktop/CCA-data/2021-1-18/HBV/top20_pathway.csv")
  
  
  ##-----------------------format pathway names ----------------------------------
  rownames(mat) <- gsub("_", " ", rownames(mat))
  rownames(mat) <- gsub("KEGG ", "", rownames(mat))
  
  ##-----------------------format pathway names 大小写转换 -----------------------
  pathway_name = rownames(mat)
  pathway_name = tolower(pathway_name) ## 将x中的字符全部转换为小写字母
  pathway_name = capitalize(pathway_name) ## 将y中的字符全部转换为大写字母
  rownames(mat) = pathway_name
  
  if (TRUE) {
    topAnn = HeatmapAnnotation(df = meta[ ,"type", drop=F], 
                               col = list(type=c("iCCA"="#2080C3","eCCA"="#89BCDF")),
                               annotation_height = unit(1, "cm"))
    
    heat.col <-  colorRamp2(c(-0.6, 0, 0.6), c('#2166ac', '#f7f7f7', 'red'))  ## 原始图例c('#2166ac', '#f7f7f7', '#b2182b'))
    
    ht <- Heatmap(mat, name="GSVA score", col = heat.col, top_annotation = topAnn, 
                  cluster_rows = F, cluster_columns = F, show_column_names = F, 
                  row_names_side = "left")
    ## pdf oupput pathway
    pdf(paste(getwd(), geneset,paste('top_20_pathway', '.pdf' ,sep = ''), sep = '/'), height=16, width = 12)
    
    maxChar <- rownames(mat)[nchar(rownames(mat))==max(nchar(rownames(mat)))]
    
    padding <- unit.c(unit(2, "mm"), 
                      grobWidth(textGrob(maxChar))-unit(50, "mm"),
                      unit(c(2, 2), "mm"))
    draw(ht, padding = padding, merge_legends = TRUE)
    dev.off()
    
  }
}

geneset = 'GO_BP'   ## geneset in c('KEGG', 'Hall_mark', 'Reactome', 'GO_BP')
try(ssGSEA(matrix, geneset))



input_path <- '/Users/ranpeng/Desktop/CCA/Data/Fig2/multiomics/multiomics_protein_delna30%.csv'
anno_path <- '/Users/ranpeng/Desktop/CCA/Data/Fig2/Arm/6q14.3_boxplot_ana.csv'
input_data <- read.csv(input_path, header=T, row.names=1, stringsAsFactors = F)
meta <- read.csv(anno_path)
rownames(meta) = colnames(input_data)
meta$id = gsub("#",".", meta$id )
matrix <- as.matrix(input_data)

## 数据存储路径
setwd('/Users/ranpeng/Desktop/CCA/Data/Fig2/Arm/6q14.3/GSVA')
outpath = getwd()

## 创建路径
KEGG_result <-paste(outpath,'KEGG',sep = '/')
Hall_mark_result <- paste(outpath,'Hall_mark',sep = '/')
Reactome_result <- paste(outpath,'Reactome',sep = '/')
GO_BP_result<- paste(outpath,'GO_BP',sep = '/')

if (file.exists(KEGG_result,Hall_mark_result,Reactome_result,GO_BP_result)){
  print('Path already exists')
} else {
  dir.create(file.path(KEGG_result))
  dir.create(file.path(Hall_mark_result))
  dir.create(file.path(Reactome_result))
  dir.create(file.path(GO_BP_result))
}

## 新建list 存储数据集
temp<-list()
temp[["KEGG"]]<-KEGG
temp[["Hall_mark"]]<-Hall_mark
temp[["Reactome"]]<-Reactome
names(temp)<-c('KEGG', 'Hall_mark', 'Reactome')

for(geneset in names(temp)){
  print(temp[geneset])
}


## 运行函数
for (geneset in c('KEGG', 'Hall_mark', 'Reactome', 'GO_BP')){
  geneset = temp[[geneset]]
  try(ssGSEA(matrix, geneset))
}

geneset = 'Reactome'
try(ssGSEA(matrix, geneset))


##----------------------- pathway 可视化 ----------------------------------
linbary(pheatmap)
chending_color = c(colorRampPalette(c("#395F90", "white"))(20),  ## #184D87  常用#1E90FF ##14417C ## 0772BC ## #395F90黑蓝
                   colorRampPalette(c("white", "red"))(20) )

breaks = unique(c(seq(-1.5, 0, length = 21), 0, seq(0,1.5, length = 21)))
pheatmap(mat,cluster_rows=0,cluster_cols=0,
         clustering_distance_cols = "correlation",fill = T, breaks=breaks,
         clustering_distance_rows = "correlation",border_color ="white", na_col = "white",
         col=chending_color,show_rownames=T,show_colnames=F,display_numbers=F)


# meta <- data.frame(id= colnames(df),
#                    type=c(rep("HBV", 153),
#                           rep("Non_HBV", 49)),
#                    stringsAsFactors = F)
# rownames(meta) <- colnames(df)