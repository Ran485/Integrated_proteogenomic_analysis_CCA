### 清空数据
rm(list = ls())
# options(digits = 10)
### 载入所需要的包
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(readxl)
### 数据加载
data = read_excel('/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/response_R1C3_ARID1A/ARID1A_index.xlsx',sheet=1)
# data = read.csv("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/FGFR_alteration_OS.csv")

colnames(data)
surv <- function(info_data){
  sfit <- survfit(Surv(time, status) ~ iCCA_Pathology , data = data)
  splots <- list()
  splots[[1]] = ggsurvplot(sfit, conf.int=F, pval=TRUE,risk.table = TRUE, risk.table.col = "strata",
                           # xlim = c(0,60),
                           palette = c("#FF9E29", "#1E92D3", "#F94F21", "#916CA0", "#599BAD", "#DBD289", "#FB9A99", "#A4D873", "#99CAE0"))
  res<-arrange_ggsurvplots(splots, print = F,ncol = 1, nrow = 1, risk.table.height = 0.25)
  ggsave(res, file = "/Users/ranpeng/Desktop/CCP_demo/CV_protein_abundance_scatter.pdf", width = 8, height = 5.5)
}

surv(data)

data = data[data$Mutation_Alll !='Co_Mut',]
data = data[data$X6q14.3_Deletion_N!='Amp',]
exprSet = data
### 对需要做生存分析的样本分组，把连续变量变成分类变量，这里选择测试的基因是GATA3,这里使用中位数
# group = ifelse(exprSet$SH3BGRL2 >= median(exprSet$SH3BGRL2),'zHigh','Low')
# group = exprSet$BAP1
# exprSet = as.data.frame(data)
exprSet = exprSet[exprSet$time <= 50,]
colnames(exprSet)
## 首先利用survfit函数拟合得到生存对象 sfit,
sfit <- survfit(Surv(time, status) ~ BAP1, data=exprSet)
p = ggsurvplot(sfit, conf.int=F, pval=TRUE,risk.table = TRUE, risk.table.col = "strata",
               # xlim = c(0,60),
               palette = c("#3C5588","#0E9F87","#1E92D3","#FF9E29","#F94F21", "#1E92D3", "#F94F21", "#916CA0", "#599BAD", "#DBD289", "#FB9A99", "#A4D873", "#99CAE0"))
p
# "#3C5588","#0E9F87","#1E92D3","#FF9E29","#F94F21",
# "#0E9F87", "#3C5588"
# 添加COX回归hazard ratio值等相关信息**
res_cox<-coxph(Surv(time, status) ~BAP1, data=exprSet)
## 获取样本注释信息
sample_count = data.frame(table(exprSet$BAP1))
surv_median_os = data.frame(surv_median(sfit))
surv_median_os = cbind(sample_count,surv_median_os)
p$plot = p$plot + 
  ## 添加HR info
  ggplot2::annotate("text",x = 15, y = 0.12,size = 5,color='gray',family = "Arial",fontface = "plain",label = paste("HR :",round(summary(res_cox)$conf.int[1],2), "(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = "")) +
  # ggplot2::annotate("text",x = 10, y = 0.10,size = 5,color='gray',family = "Arial",fontface = "plain",label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))+
  ggplot2::annotate("text",x = 15, y = 0.05, size = 5,color='gray',family = "Arial",fontface = "plain",label = paste("logrank-test pval:",round(summary(res_cox)$coef[5],4))) +
  ## 添加中位生存 info
  ggplot2::annotate("text",x = 40, y = 0.9,size = 4,color='#AB1420', family = "Arial", fontface = "bold",label = paste("--- ",surv_median_os[1,1],", n = ",round(surv_median_os[1,2],3),", ","MST = ",round(surv_median_os[1,4],3)," months",sep = ""))+
  ggplot2::annotate("text",x = 40, y = 0.8,size = 4,color='#0070A6', family = "Arial", fontface = "bold",label = paste("--- ",surv_median_os[2,1],", n = ",round(surv_median_os[2,2],3),", ","MST = ",round(surv_median_os[2,4],1)," months",sep = ""))

p
