library(ggplot2)
library(magrittr)
library(ggpubr)
library(cowplot)
library(tidyr)
library('openxlsx')
# library(showtext) ## handle 中文乱码
# showtext_auto()

## ---------------- Defined functions ------------------------

# handle outliers
capping_outliers_IQR <- function(df,category_var, factor = 4){
  qnt <- quantile(df[,category_var], probs=c(.25, .75), na.rm = T)
  caps <- quantile(df[,category_var], probs=c(.10, .90), na.rm = T)
  H <- factor * IQR(df[,category_var], na.rm = T)
  lower_bound = qnt[1] - H
  upper_bound = qnt[2] + H
  df[,category_var] = ifelse(df[,category_var] < lower_bound, caps[1], df[,category_var])
  df[,category_var] = ifelse(df[,category_var] > upper_bound, caps[2], df[,category_var])
  return(df)
}

# plot boxplots
customed_boxplot<-function(data, target = "Feature",category_var_name = "SOD1", log2_transform = FALSE, capping_outliers = FALSE, violin_plot = FALSE){
  df = data[,c(target,category_var_name)]
  if (capping_outliers) {
  df1 = capping_outliers_IQR(df,category_var = category_var_name)
  } else {df1 = df}
  if (log2_transform) {
    df1[,category_var_name] = log2(df1[,category_var_name] + 1)}
  if (violin_plot) {
    ggplot(df1, aes_string(x=target, y=category_var_name,color=target)) + 
      geom_violin(trim=FALSE) + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(以下设为背景为白色，其实表示不要轮廓线)
      #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
      geom_boxplot(width=0.6,position=position_dodge(0.9),outlier.colour = "red")+ #绘制箱线图
      # geom_jitter(alpha= 0.7,aes(x=Group, y=value,colour = Attribute),width = 0.15,height = 0) +
      # scale_fill_manual(values = c("#0E9F87", "#3C5588"))+ #设置填充的颜色
      scale_color_manual(values=c("#0E9F87", "#3C5588")) +
      theme_bw()+ #背景变为白色
      stat_summary(fun.y=mean, geom="point", shape=23, size=2,position = position_dodge(width = 0.9)) + # violin plot with mean points
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
      # ggtitle("Cis efffect on Protein")+
      ggtitle(category_var_name)+
      # ylim(0,1.5)+
      theme(axis.text.x=element_text(angle=30,hjust = 0.6,colour="black",family="Arial",size=16), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
            axis.title.x=element_text(angle=0,hjust = 0.5,colour="black",family="Arial",size=16),
            axis.text.y=element_text(family="Arial",size=16,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
            axis.title.y=element_text(family="Arial",size = 16,face="plain"), #设置y轴标题的字体属性
            panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.8), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
            legend.text=element_text(face="plain", family="Arial", colour="black",  #设置图例的子标题的字体属性
                                     size=12),
            legend.title=element_text(face="plain", family="Arial", colour="black", #设置图例的总标题的字体属性
                                      size=12),
            panel.grid.major = element_blank(),   #不显示网格线
            panel.grid.minor = element_blank())+  #不显示网格线
      ylab("Relative exprssion [log2]")+xlab("") + #设置x轴和y轴的标题
      # stat_compare_means(aes(group = Attribute), label = "p.format")  # 添加统计显著型水平
      stat_compare_means(aes_string(group = target),method = "wilcox",label = "p.format", label.y = 4.5) +       # Add global annova p-value
      stat_compare_means(aes_string(group = target),label = "p.signif", method = "wilcox")  # Pairwise comparison against all t.test wilcox
   } else {
    ggboxplot(df1, x = target, y = category_var_name,combine = F,
              # add = c("mean_se"),
              color = target ,
              fill = target, 
              alpha = 0.12,
              order =  c("WT","zMutation"), #,"Cluster_3","Cluster_4","Cluster_5","Cluster_6"),
              # palette = c("#0E9F87", "#3C5588","#FF9E29","#86AA00","#F94F21","#916CA0","#599BAD","#DBD289"), ##palette = c("#F94F21","#EFC000","#00AEBA"),
              palette = c("#0E9F87", "#3C5588"),
              width = 0.5 ,x.text.angle = 45)+
      rotate_x_text(angle=45,hjust=.5,vjust=1) +
      # geom_violin(trim=FALSE) +
      ## Global p-value
      stat_compare_means(label.y = c(), paired = FALSE) + #, method = "t.test") + # "t.test" 
      # stat_compare_means(aes(group = Type),  label.y = c()) +
      ## 设置标题、标注
      labs(title="", x="", y = "Relative exprssion [log2]")  +
      ## 设置标题、标注
      geom_jitter(alpha= 0.7,aes_string(colour = target),width = 0.15,height = 0) +
      theme_test()+
      # coord_flip()
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
      ggtitle(category_var_name)+ # 设置变量名
      theme(axis.text.x=element_text(angle=18,hjust = 0.6,colour="black",family="Arial",size=16), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
            axis.title.x=element_text(angle=0,hjust = 0.5,colour="black",family="Arial",size=16),
            axis.text.y=element_text(family="Arial",size=16,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
            axis.title.y=element_text(family="Arial",size = 16,face="plain"), #设置y轴标题的字体属性
            panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.8), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
            legend.text=element_text(face="plain", family="Arial", colour="black",  #设置图例的子标题的字体属性
                                     size=12),
            legend.title=element_text(face="plain", family="Arial", colour="black", #设置图例的总标题的字体属性
                                      size=14),
            panel.grid.major = element_blank(),   #不显示网格线
            panel.grid.minor = element_blank())
  }
  ## 设置标题、标注
  # labs(
  #   title = 'category_var_colname',
  #   subtitle = 'fisher_test_pval',
  #   tag = "",
  #   caption = "",
  #   x = "",
  #   y = "Frequency [%]"
  # )
  # coord_flip()
}

# ## ------------------------ input data -----------------------------
# data <-read.csv("/Users/ranpeng/Desktop/Desktop/项目文件/对外服务/罗俊一/script/results/clinical_EDA/clinical_washed_merged_xgboost_boxplot.csv",header=T)
# # data<-read.xlsx("/Users/ranpeng/Desktop/Desktop/项目文件/对外服务/罗俊一/clinical_info/clinical_info_arranged/Clinical_continuous_variables.xlsx",sheet=1,colNames = T)
# # data = as.data.frame(data)
# # df = data[,c("Feature","SOD1")]
# colnames(data)
# 
# ## ----------------- 指定 category_var_colname 进行分析 --------------------
# customed_boxplot(data, target = "LBC_circulation",category_var_name = "Basophil_count",log2_transform = F)
# 
# ## --------- for loop 对指定的 category_var_colname list 批量进行分析 -----------
# colnames(data)
# # category_list = c("Heart_failure","Gender","Ethnicity","Comorbidities_modified","Smoking","Family_History","Lateral_branch_blood_supply","Processing","ACEIorARB","CCB","ARNI","Nitrates")
# category_list = colnames(data)[-c(1:18)]
# # category_list = c('Basophil_count', 'HBP', 'Plasma_fibrinogen', 'Chlorine')
# category_var_plots = list()
# 
# for (category_var in category_list){
#   category_var_plots[[category_var]] = customed_boxplot(data, target = 'LBC_circulation',category_var_name = category_var, log2_transform = F)
#   print(category_var_plots[[category_var]])
#   # ggsave(category_var_plots[[category_var]], file=paste0("plot_", category_var,".png"), width = 44.45, height = 27.78, units = "cm", dpi=300)
# }
# 
# ## 生成多个图关键的一步
# plot_grid(plotlist = category_var_plots)






