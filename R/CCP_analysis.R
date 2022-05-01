library(ConsensusClusterPlus)
library(pheatmap)

CCP <- function(matrix, distance, clusterAlg){
  rcc = ConsensusClusterPlus(
    matrix,
    maxK=7,
    reps=1000,
    pItem=0.8,
    pFeature=1,
    distance=distance,
    clusterAlg=clusterAlg,
    title=paste(distance, clusterAlg, sep='_'),
    plot='pdf'
  )
  
  for (i in c(2, 3, 4, 5, 6, 7)){
    data1 <- data.frame(rcc[[i]]$consensusMatrix, row.names = colnames(matrix))
    colnames(data1) <- colnames(matrix)
    
    anno2 <- as.data.frame(rcc[[i]]$consensusClass)
    anno2 <- cbind(clinical_annot, anno2)
    colnames(anno2)[length(colnames(anno2))] <- c('CCP_Class') ## 临床信息mapping
    anno2 <- anno2[order(anno2$CCP_Class),]
    write.csv(anno2, paste(getwd(),paste(distance, clusterAlg,sep = '_'), paste(i, '.csv' ,sep = '') ,sep = '/'))
    pdf(paste(getwd(), paste(distance, clusterAlg,sep = '_'), paste(i, '.pdf' ,sep = '') ,sep = '/'), height=7, width = 8)
    # pdf('/Users/ranpeng/Desktop/Xcell_score/results(2021-1-25)/1.pdf', height=6, width = 8)
    heatmap <- pheatmap(
      data1[row.names(anno2),row.names(anno2)],
      color=colorRampPalette(c('white', 'blue'))(51),
      # clustering_method='average',
      show_rownames = F,
      show_colnames = F,
      annotation_col=anno2,
      border_color='white',
      cluster_cols = F,
      cluster_rows = F,
      # annotation_colors = list(Cluster=c('1'='red', '2'='blue', '3'='green', '4'='black', '5'='orange')),
      # clustering_distance_rows = 'correlation',
      # clustering_distance_cols = 'correlation',
      angle_col=0
    )
    print(heatmap)
    dev.off()
  }
}

path <- '/Users/ranpeng/Desktop/CCA/Data/Fig4/phosphosite_CCP/phospho_CCP_input_sd_top1000.csv'
anno_path <- '/Users/ranpeng/Desktop/CCA/Data/Fig4/phosphosite_CCP/ana_index.csv'
sor_ratio <- read.csv(path, header=T, row.names=1, stringsAsFactors = F)
clinical_annot <- read.csv(anno_path, header=T, row.names=1)
matrix <- as.matrix(sor_ratio)
colnames(matrix) <- rownames(clinical_annot)
matrix[is.na(matrix)] = 0.00001
setwd('/Users/ranpeng/Desktop/CCA/Data/Fig4/phosphosite_CCP/results_rep1000')

for (distance in c('pearson', 'spearman', 'euclidean', 'binary', 'canberra', 'minkowski')){
  for (clusterAlg in c('pam','km', 'hc')){
    try(CCP(matrix, distance, clusterAlg))
  }
}






