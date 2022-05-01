## define fuction
fisher_exact_test <- function(data, target, features, pval_cutoff=FALSE, cutoff_thresh=0.05){
  pval_res = c()
  test_features = c()
  for (i in features) 
  {
    if (length(unique(data[,i])) == 2) {
    res<-fisher.test(table(data[,target], data[,i]))
    pval_res <- c(pval_res, res$p.value)
    test_features <- c(test_features,i)
    }
  }
  result_merge = as.data.frame(test_features)
  result_merge["p_value"] <-  pval_res
  result_merge[,"FDR"] <- p.adjust(result_merge[,"p_value"],method="BH",length(result_merge[,"p_value"]))
  if (pval_cutoff) {
    fisher_significant = subset(result_merge, result_merge[,2] < cutoff_thresh )
    return(fisher_significant)
  } else {
    return(result_merge)
  }
  }

# data = read.csv('/Users/ranpeng/Desktop/Integrated_multiomics_analysis_CCA/data/Fig6/LS_arm_events.csv')
# colnames(data)
# features = colnames(data)[4:ncol(data)]
# significant_res = fisher_exact_test(data, target='iCCA_Pathology', features = features)

# write.csv(significant_res,"/Users/ranpeng/Desktop/Integrated_multiomics_analysis_CCA/data/Fig6/significant_LS_duct_mutation_events.csv")

