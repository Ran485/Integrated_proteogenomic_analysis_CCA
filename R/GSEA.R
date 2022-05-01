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
