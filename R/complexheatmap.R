# ------------- 绘制选择基因注释的heatmap ---------------
#读取的示例文件
#gene.txt是基因表达值矩阵
#name.txt是带展示名称的基因列表

#加载ComplexHeatmap包绘制这种类型的热图
library(ComplexHeatmap)
library(circlize)
library(BBmisc)
library(readxl)

Z_Score <- function(x){
  res <- (x - mean(x)) / sd(x)
  return(res)}

setwd("/Volumes/Samsung_T5/downloads/CCA/Data/Fig4")
#表达矩阵，考虑到ComplexHeatmap没有scale参数，因此需事先手动做个行标准化
# mat <- read.csv('NMF_heatmap.csv', row.names = 1, check.names = FALSE)
data <- read.csv('/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/response_R2C5/protein_marker_geneheatmap.csv', row.names = 1, check.names = FALSE)
#定义样本分组，例如示例文件中共 A、B、C 3组
# samples <- rep(c('S-I', 'S-II', 'S-III'), c(93, 46, 77))
meta = read.csv("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/response_R2C5/heatmap_anatation.csv")

## 如果存在空值进行下面的处理
mat = data
mat<-mat[,-which(apply(mat,2,function(x) all(is.na(x))))]
meta <- meta[which(meta$id %in% colnames(mat)),]
groups = meta$type
# 数据转换
mat[is.na(mat)] = 0.0001
mat = log2(mat)
# for (i in 1:nrow(mat)) mat[i, ] <- scale(log(unlist(mat[i, ]) + 1, 2))
## ------  对数据求Z-score  -----------
mat<- apply(mat, 1, Z_Score)
mat = t(mat)
mat <- as.matrix(mat)

#先绘制一个不显示任何基因名称的热图
heatmap_breakNum <- 2 # 设定颜色阈值
heatmap_col <- c('#2080C3', '#f7f7f7', 'red') #(100) #定义热图由低值到高值的渐变颜色
heat.col <-  colorRamp2(c(-heatmap_breakNum, 0, heatmap_breakNum), heatmap_col)
heat <- Heatmap(mat, 
                col = heat.col,
                heatmap_legend_param = list(grid_height = unit(10,'mm')),  #图例高度设置
                show_row_names = FALSE,  #不展示基因名称
                show_column_names = FALSE,
                cluster_rows = FALSE, 
                cluster_columns = FALSE,
                top_annotation = HeatmapAnnotation(Group = groups, 
                                                   simple_anno_size = unit(2, 'mm'), 
                                                   # col = list(Group = c('S-I' = '#00DAE0', 'S-II' = '#FF9289', 'S-III' = 'blue')),  #定义样本分组的颜色
                                                   show_annotation_name = FALSE), 
                column_names_gp = gpar(fontsize = 10), row_names_gp = gpar(fontsize = 6))

heat


#读取待展示的基因名称，并添加到热图中
name <- read.csv('/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/response_R2C5/heatmap_gene_anatation.csv', header = T, check.names = FALSE)
heat + rowAnnotation(link = anno_mark(at = which(rownames(mat) %in% name$gene), 
                                      labels = name$gene, labels_gp = gpar(fontsize = 10)))






