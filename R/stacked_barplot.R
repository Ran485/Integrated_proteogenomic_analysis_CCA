## --------------- 百分数占比 barplot ----------------
library(vcd)
library(ggplot2)
library(Hmisc)
library(cowplot)
library('openxlsx')
# library(showtext) ## handle 中文乱码
# showtext_auto()

## --------------------- define stacked_barplot function --------------------
stacked_barplot <- function(data = data, target = 'LBC_circulation',category_var_colname = 'Lateral_branch_blood_supply',binary=FALSE){
  # target = 'LBC_circulation';category_var_colname = 'Lateral_branch_blood_supply'  
  ## 统计个数
  if (binary == TRUE){
  col_1 = paste0(category_var_colname, "_delNA")
  data[,col_1] = ifelse(data[,category_var_colname] == 0, 0, 1)
  # data[,col_1] = data[,category_var_colname]
  data[,col_1] = factor(data[,col_1],levels=c('0','1'),labels=c(paste0('Non_', category_var_colname), capitalize(category_var_colname)),ordered=TRUE)
  } else {
    col_1 = paste0(category_var_colname, "_delNA")
    data[,col_1] = data[,category_var_colname] }
  # data[,col_1] = factor(data[,col_1],levels=c('0','1','2','3'),labels=c('Grade I','Grade II','Grade III','Grade IV'),ordered=TRUE)
  ggplot(data=data, mapping=aes_string(x=target,fill=col_1))+
    geom_bar(stat="count",width=0.5,position='fill')+
    scale_fill_manual(values=c("#0E9F87","#3C5588","#FF9E29","#86AA00","#F94F21","#916CA0","#599BAD","#DBD289"))+  ## '#999999','#E69F00'
    geom_text(stat='count',aes(label=..count..), color="white", size=3.5,position=position_fill(0.5))+
    theme_minimal() #+coord_flip()
  
  ## delete contain NA rows in the dataframe
  temp = data[,c(target,col_1)]
  temp = na.omit(temp)
  
  ## ------------- Fisher's Exact Test for Count Data --------------------------
  # data = read.csv("/Users/ranpeng/Desktop/CCA/Data/Fig4/gene_mutation_results/CCA_CNAs_del.csv")
  res = table(temp[,target],temp[,col_1])
  print(res)
  print(fisher.test(res))
  p_val = fisher.test(res)$p.value
  fisher_test_pval = paste0('fisher.test p_value',' = ', signif(p_val,4))
  ## 百分数占比
  ggplot(data=temp, mapping=aes_string(x=target,fill=col_1))+
    geom_bar(stat="count",width=0.7,position='fill')+
    scale_fill_manual(values=c("#3C5588","#0E9F87","#FF9E29","#86AA00","#F94F21","#916CA0","#599BAD","#DBD289"))+
    # geom_text(stat='count' ,aes(label=scales::percent(..count../sum(..count..)))
    #           , color="white", size=5,position=position_fill(0.5))+
    geom_text(stat='count' ,aes(label=..count..)
              , color="white", size=5,position=position_fill(0.5))+
    # theme_minimal() #+coord_flip()
    # geom_text(aes(label=`Quantative protein numbers [%]`), vjust=1.6, color="white", size=6)+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    ggtitle(title)+ # 设置变量名
    theme(axis.text.x=element_text(angle=18,hjust = 0.6,colour="black",family="Arial",size=16), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
          axis.title.x=element_text(angle=0,hjust = 0.5,colour="black",family="Arial",size=16),
          axis.text.y=element_text(family="Arial",size=16,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
          axis.title.y=element_text(family="Arial",size = 20,face="plain"), #设置y轴标题的字体属性
          panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.8), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
          legend.text=element_text(face="plain", family="Arial", colour="black",  #设置图例的子标题的字体属性
                                   size=12),
          legend.title=element_text(face="plain", family="Arial", colour="black", #设置图例的总标题的字体属性
                                    size=14),
          panel.grid.major = element_blank(),   #不显示网格线
          panel.grid.minor = element_blank()) +
    ## 设置标题、标注
    labs(
      title = category_var_colname,
      subtitle = fisher_test_pval,
      tag = "",
      caption = "",
      x = "",
      y = "Frequency [%]"
    )
  # geom_hline(yintercept = 25,color = 'gray', linetype="dashed")
}

## ------------------------ input data -----------------------------
# data = read.csv('/Users/ranpeng/Desktop/CCA/Data/Fig4/BAP1_mutation/BAP1_mutation_barplot.csv')
# data<-read.xlsx("/Users/ranpeng/Desktop/Desktop/项目文件/对外服务/罗俊一/clinical_info/clinical_info_arranged/clinical_info_mapped_exp_select_20220210.xlsx",sheet=1,colNames = T)
# data = read.csv('/Users/ranpeng/Desktop/Desktop/项目文件/对外服务/罗俊一/clinical_info/clinical_info_arranged/Clinical_category_variables.csv')
# colnames(data)

## ----------------- 指定 category_var_colname 进行分析 --------------------
# stacked_barplot(data, target = 'LBC_circulation',category_var_colname = 'Lateral_branch_blood_supply')
# 
# ## --------- for loop 对指定的 category_var_colname list 批量进行分析 -----------
# # category_list = c("Heart_failure","Gender","Ethnicity","Comorbidities_modified","Smoking","Family_History","Lateral_branch_blood_supply","Processing","ACEIorARB","CCB","ARNI","Nitrates")
# category_list = significant_res$test_features
# # category_list = c('SMAD4','SETBP1','FAF1','TYRO3')
# category_var_plots = list()
# for(category_var in category_list){
#   category_var_plots[[category_var]] = stacked_barplot(data, target = 'iCCA_Pathology',category_var_colname = category_var, binary = F)
#   print(category_var_plots[[category_var]])
#   # ggsave(category_var_plots[[category_var]], file=paste0("plot_", category_var,".png"), width = 44.45, height = 27.78, units = "cm", dpi=300)
# }
# 
# ## 生成多个图关键的一步
# plot_grid(plotlist = category_var_plots)


