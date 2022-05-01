# ---------- define fuction -----------
Acquire_sig_gene <- function(pathway_name,DEP){
  ## pathway_name格式:
  # pathway
  # KEGG_SPLICEOSOME
  # KEGG_PROTEASOME
  ## DEP 格式：
  # gene
  # HNRNPH2
  # C8orf82
##-----------------------import package -----------------------------------------
library("plyr")  #加载获取rbind.fill函数
library("dplyr")  
##-----------------------load gmt file -----------------------------------------
KEGG <- read.table('/Users/ranpeng/Desktop/CCA/database/msigdb/c2.cp.kegg.v7.2.symbols.txt',sep = "\t",header = F, fill = TRUE)
Hall_mark <- read.table('/Users/ranpeng/Desktop/CCA/database/msigdb/h.all.v7.2.symbols.txt',sep = "\t",header = F, fill = TRUE)
Reactome <- read.table('/Users/ranpeng/Desktop/CCA/database/msigdb/c2.cp.reactome.v7.2.symbols.txt',sep = "\t",header = F, fill = TRUE)
GO_BP <- read.table('/Users/ranpeng/Desktop/CCA/database/msigdb/c5.go.bp.v7.2.symbols.txt',sep = "\t",header = F, fill = TRUE)

#不等长数据合并
list1<-list()
list1[[1]]=KEGG
list1[[2]]=Hall_mark
list1[[3]]=Reactome
list1[[4]]=GO_BP
database = do.call(rbind.fill,list1)
rownames(database) = database[,1]
database = database[,-c(1,2)]

##-----------------------load dif pathway -----------------------------------------
pathway_name = pathway_name
data = database[pathway_name$pathway,]
data = t(data)
##-----------------------load dif protein -----------------------------------------
DEP = DEP

## ----------------------查找significant gene --------------------------------------
list = list()
for(i in 1:ncol(data)){
list[[i]] = data.frame(t(intersect(data[,i],DEP$gene)))
rownames(list[[i]]) = rownames(data)[i]
}

res = do.call(rbind.fill,list)
rownames(res) = colnames(data)
return(res)

}


##-----------------------load dif pathway -----------------------------------------
pathway_name = read.csv('/Users/ranpeng/Desktop/CCA/Data/Fig7/Xcell_signature_cor/Neutrophils/select_pathway.csv')
pathway_name$pathway = toupper(pathway_name$pathway )
##-----------------------load dif protein -----------------------------------------
DEP = read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig7/Xcell_signature_cor/Neutrophils/Neutrophils_GSEA_cor.csv")

##----------------------- run function -----------------------------------------
res = Acquire_sig_gene(pathway_name,DEP)
res = t(res)
write.csv(res,"/Users/ranpeng/Desktop/CCA/Data/Fig7/Xcell_signature_cor/Neutrophils/pathway_mapping_data.csv")




