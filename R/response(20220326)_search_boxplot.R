## ------------ boxplot search gene ----------------
# ## protein
# df = read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig4/BAP1_mutation/ROS score/ROS_score_caculate.csv",row.names = 1)
# df = data.frame(t(df))
# ana = read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig3/anatation_boxplot.csv")
# # ## RNA TN_all
# df = read.csv("/Users/ranpeng/Desktop/CCA/Metadata/RNA/RNA_FPKM_codegene_deduplicate.csv",row.names = 1)
# df = data.frame(t(df))
# ana = read.csv("/Users/ranpeng/Desktop/CCA/Metadata/RNA/anatation.csv")
# ## RNA Tumor
df = read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig2/multiomics/multiomics_RNA_delna30%.csv",row.names = 1)
df = data.frame(t(df))
# ana = read.csv("/Users/ranpeng/Desktop/CCA/Metadata/RNA/anatation.csv")

# # protein
df = read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig2/multiomics/multiomics_protein_delna30%.csv",row.names = 1)
df = data.frame(t(df))
# protein
ana = read.csv("/Users/ranpeng/Desktop/Integrated_multiomics_analysis_CCA/data/Fig1/fig1b/DPCR1_anatation.csv")
# RNA
ana = read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig7/RNA/arm_ana.csv")


## ---------------- define function ----------------
Z_Score <- function(x){
  res <- (x - mean(x)) / sd(x)
  return(res)}

search_boxplot <- function(gene){
  # search = paste(gene,"$|^",gene,".y",sep = '')
  search = paste("^",gene,"$", sep="") ## 完全匹配
  index=grep(search, colnames(df))  
  # index=grep(gene, colnames(df)) 
  print(colnames(df)[index])
  result_data=df[,index];result_data = data.frame(result_data);result_data$id = rownames(df)
  
  ## merge group
  ana$id = gsub("#",".",ana$id)
  result_data = merge(ana,result_data)
  result_data$result_data = as.numeric(result_data$result_data)
  # ## filter low expression row
  # result_data = subset(result_data,result_data < 4 & result_data > -4)
  ## data transform b不需要则注释掉
  # result_data$result_data = log10(result_data$result_data)
  
  ## z-score
  # result_data$result_data = Z_Score(result_data$result_data)
  # result_data = subset(result_data,result_data < 0.5 & result_data > -0.5)
  ## boxplot
  library(ggpubr)
  result_data[result_data==0] = NA
  result_data = na.omit(result_data)
  ggboxplot(result_data, x = "type", y = "result_data",combine = F,
            # add = c("mean_se"),
            color = "type" ,
            fill = "type", alpha = 0.12,
            # order =  c("Normal","Tumor"),#,"Cluster_3","Cluster_4","Cluster_5","Cluster_6"),
            palette = c("#86AA00","#F94F21","#916CA0","#599BAD","#DBD289"), ##palette = c("#F94F21","#EFC000","#00AEBA"),
            width = 0.5 ,x.text.angle = 45)+
    # rotate_x_text(angle=15,hjust=.5,vjust=1) + #FF9E29",
    stat_compare_means(label.y = c()) +                                         # Global p-value
    # stat_compare_means(aes(group = Type),  label.y = c()) +
    labs(title= gene, x="", y = "Relative Expression")  +
    geom_jitter(alpha= 0.7,aes(colour = type),width = 0.15,height = 0) +
    theme_test()+
    theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))  ## 控制坐标轴文本的角度
  
}

search_boxplot("MUCL3")
colnames(df)