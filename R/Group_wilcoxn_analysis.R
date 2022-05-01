library(dplyr)
library(greenbrown) 
library('openxlsx')
library(xlsx)
###—————————————— wilcox.test and p.val-adjust function  ——————————————————————
# handle outliers
capping_outliers_IQR <- function(df,category_var, factor = 1.5){
  qnt <- quantile(df[,category_var], probs=c(.25, .75), na.rm = T)
  caps <- quantile(df[,category_var], probs=c(.10, .90), na.rm = T)
  H <- factor * IQR(df[,category_var], na.rm = T)
  lower_bound = qnt[1] - H
  upper_bound = qnt[2] + H
  df[,category_var] = ifelse(df[,category_var] < lower_bound, caps[1], df[,category_var])
  df[,category_var] = ifelse(df[,category_var] > upper_bound, caps[2], df[,category_var])
  return(df)
}


group_DEP <- function(data,ana_index,capping_outliers = FALSE, group1,group2,method){
  if (capping_outliers) {data = apply(data, 1, capping_outliers_IQR)}
  index = ana_index[ana_index[,2] == group1|ana_index[,2] == group2, ]
  print(index)
  group_info = data.frame(table(index$type)); rownames(index) = index[,1]
  print(group_info)
  sep1 = group_info[, group_info$Var1 == group1][2]
  sep2 = sep1 + 1
  cat("sep1: ", sep1, '\n')
  cat("sep2: ", sep2, '\n' )
  # library(dplyr)
  index_col = intersect(index$id,colnames(data))
  expset = data[,index_col]
  index_1 = group1
  index_2 = group2
  index_group1 = index[index[,2] == group1,]; index_group1 = index_group1$id
  group1_expset = expset[,index_group1]
  # print(group1_expset)
  index_group2 = index[index[,2] == group2,]; index_group2 = index_group2$id
  group2_expset = expset[,index_group2]
  # print(group2_expset)
  expset = cbind(group1_expset, group2_expset)
  # print(expset)
  #  optional to choose calculate the difference method
  if (method == 't_test') {
    ### ————————————— t-test p-value  ——————————————————————————————-
    tstat.val = apply(expset, 1, function(x) { if(AllEqual(x[1:ncol(expset)])){return(1)}
      else if(AllEqual(x[1:sep1]) & AllEqual(x[sep2:ncol(expset)])){return(0)} else t.test(x[1:sep1], x[sep2:ncol(expset)],paired=F)$statistic})
    tstat.pval = apply(expset, 1, function(x) { if(AllEqual(x[1:ncol(expset)])){return(1)}
      else if(AllEqual(x[1:sep1]) & AllEqual(x[sep2:ncol(expset)])){return(0)} else t.test(x[1:sep1], x[sep2:ncol(expset)],paired=F)$p.value})
    p_val <- data.frame(tstat.val, tstat.pval)
    p_val[,"FDR_t-test"] <- NA
    p_val[,"FDR_t-test"] <- p.adjust(p_val[,"tstat.pval"],method="BH",length(p_val[,"tstat.pval"]))
  } else if (method == 'wilcoxn_test') {
    ### ————————————— wilcoxn-test p-value  ——————————————————————————————-
    wilcox_stat.val = apply(expset, 1, function(x) { wilcox.test(x[1:sep1], x[sep2:ncol(expset)],paired=F)$statistic})
    wilcox_stat_pval= apply(expset, 1, function(x) { wilcox.test(x[1:sep1],x[sep2:ncol(expset)],paired=F)$p.value})
    p_val <- data.frame(wilcox_stat.val, wilcox_stat_pval)
    p_val[,"FDR_wilcox"] <- NA
    p_val[,"FDR_wilcox"] <- p.adjust(p_val[,"wilcox_stat_pval"],method="BH",length(p_val[,"wilcox_stat.val"]))
  } else {
    print("Please choose either `t-test` or `wilcoxn-test` method")
  }
  
  ###—————————————— Caculate group average and FC  ——————————————————————
  p_val[,index_1] = apply(expset[,1:sep1], 1, mean, na.rm = T)
  p_val[,index_2] = apply(expset[,sep2:ncol(expset)], 1, mean, na.rm = T)
  p_val$count_1 = rowSums(expset[,1:sep1] >= 0, na.rm = T) ## 统计每组非 NA 值的个数
  p_val$count_2 = rowSums(expset[,sep2:ncol(expset)] >= 0, na.rm = T) ## 统计每组非 NA 值的个数
  p_val$FC = p_val[,index_2]/p_val[,index_1]
  p_val$log2_FC = log2(p_val$FC)
  cat("index_1 =",index_1,"\n");cat("index_2 =",index_2,"\n")
  cat("FC = index_2/index_1:" ,index_2,"/",index_1)
  return(p_val)
}

# data<-read.xlsx("/Users/ranpeng/Desktop/Desktop/CCA/Data/Fig4/continuous_clinical_info.xlsx", sheetIndex=1, colNames = T, rowNames = T)
data <- read.csv('/Users/ranpeng/Desktop/Desktop/CCA/Data/Fig4/data_wash/Fig_4_merge_protein_Delna30%_20210325.csv',header = T,row.names = 1,skip = 0)
ana_index = read.csv("/Users/ranpeng/Desktop/Desktop/CCA/Data/Fig4/anatation.csv",header = T)
ana_index$id = gsub("#",".",ana_index$id)
table(ana_index$type)
## 原始数据处理
# data[is.na(data)] == 0
# data = log2(data+1)
res = group_DEP(data,ana_index,group1 = "S_I",group2= "S_II", method = 'wilcoxn_test')
# Export p-value results
write.csv(res,"/Users/ranpeng/Desktop/Desktop/CCA/Data/Fig4/data_wash/results/clinical_EDA/continuous_clinical_info_DEP_t_test.csv")

