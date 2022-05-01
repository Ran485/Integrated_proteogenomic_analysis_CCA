#=========================================================================
#      Exploratory data analysis of overall survival
#=========================================================================

### 清空数据
rm(list = ls())
options(digits = 10)
### 载入所需要的包
library(survival)
library(survminer)
library(tidyverse)  ## 数据处理的包
library(DT)
library(dplyr)
setwd('/Users/ranpeng/Desktop/CCA/Data/Fig2/Arm/6q14.3/cis/cis_pro_os.csv')


### 数据加载
CCA = read.csv('/Users/ranpeng/Desktop/CCA/Data/Fig4/BAP1_mutation/PARP1/PARP1_protein_os.csv')
CCA <- CCA[-which(is.na(CCA$time)), ]
colnames(CCA)[1:15]
## 如果需要对整个CCA cohort算os，需要注释掉下面两行
CCA = CCA[CCA$iCC_eCC != 1,] ## 筛选icc和ecc
# CCA = subset(CCA,time <= 60)
CCA$Status_1 = ifelse(CCA$Status == ' dead', 1, 0)
CCA$time = CCA[,'Days']/30
data = CCA
dim(data); exprSet = data
deg_1 = exprSet; m = deg_1

### 对需要做生存分析的样本分组，把连续变量变成分类变量，使用中位数作为cut-off
# group = ifelse(exprSet$PARP1 >= 0.05,'AA','Non-AA')
group = ifelse(exprSet$PARP1 >= quantile(exprSet$PARP1, 0.66) ,'High',
               ifelse(exprSet$PARP1 <= quantile(exprSet$PARP1, 0.33),'low', 'median'))

# ifelse(exprSet$average == 0 ,'WT','Amp'))
# exprSet$SH3BGRL2_1 = ifelse(exprSet$SH3BGRL2_y>median(exprSet$SH3BGRL2_y) ,"High","low")
# ifelse(exprSet$SH3BGRL2< -0.2,"Del","WT"))

# exprSet<- CCA[!(CCA[,4]=='5'|CCA[,4]=='4'|CCA[,4]=='6'),]
sfit <- survfit(Surv(time, status) ~ group , data=exprSet)
p = ggsurvplot(sfit, conf.int=F, pval=TRUE,risk.table = TRUE, risk.table.col = "strata",
               xlim = c(0,60),
               palette = c("#FF9E29", "#1E92D3", "#F94F21", "#916CA0", "#599BAD", "#DBD289", "#FB9A99", "#A4D873", "#99CAE0"))
p
# ---------- 添加COX回归hazard ratio值等相关信息 ----------
res_cox<-coxph(Surv(time, status) ~group, data=exprSet)
## 获取样本注释信息
sample_count = data.frame(table(group))
surv_median_os = data.frame(surv_median(sfit))
surv_median_os = cbind(sample_count,surv_median_os)
p$plot = p$plot + 
  ## 添加HR info
  ggplot2::annotate("text",x = 15, y = 0.12,size = 5,color='gray',label = paste("HR :",round(summary(res_cox)$conf.int[1],2), "(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = "")) + 
  # ggplot2::annotate("text",x = 10, y = 0.10,size = 5,color='gray',label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))+
  ggplot2::annotate("text",x = 15, y = 0.05, size = 5,color='gray',label = paste("logrank-test pval:",round(summary(res_cox)$coef[5],4))) +
  ## 添加中位生存 info
  ggplot2::annotate("text",x = 40, y = 0.9,size = 4,color='#AB1420', family = "Arial", fontface = "bold",label = paste("--- ",surv_median_os[1,1],", n = ",round(surv_median_os[1,2],3),", ","MST = ",round(surv_median_os[1,4],3)," months",sep = ""))+
  ggplot2::annotate("text",x = 40, y = 0.8,size = 4,color='#0070A6', family = "Arial", fontface = "bold",label = paste("--- ",surv_median_os[2,1],", n = ",round(surv_median_os[2,2],3),", ","MST = ",round(surv_median_os[2,4],1)," months",sep = ""))

p


reg <- coxph(Surv(time, status) ~ group, data = exprSet)
summary(reg)

## 生存统计信息
data.survdiff = survdiff(Surv(time, status) ~ group, data=exprSet)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))

## 计算每个基因在所有样本中的均值，赋值到数据框df中
df <- apply(m[,4:26], 2, median)	#按列，求均值，1：2位genesymbol、geneid，后面的是每个基因表达数据
df = data.frame(df)

#df_1 = df
##根据均值，将样本分为高低表达组，此处29需要根据自己的数据修改，或改成length(colnames(deg))
for (i in 4:26) {
  n=colnames(deg_1)[i]
  deg_1[,n]=ifelse(deg_1[,n] >= df[n,], "high","low")
}
print(deg_1)

## ------------- 返回生存P-value的总表 ---------------------
# 首先得到得到生存对象 my.surv, 循环计算p值
exprSet$time = as.numeric(exprSet$time)
my.surv <- Surv(exprSet$time, exprSet$status)

# 使用apply循环,对数据的2到4列进行操作，实际上就是前三个基因，
# 这里面使了survdiff，用来比较差异大小，获得p值
# 如果存在one group 的情况，需要删除完全一样的列
# exprSet_mean = data.frame(apply(exprSet[,13:length(exprSet)], 2,mean));names(exprSet_mean) ="mean"
# list = rownames(subset(exprSet_mean,exprSet_mean$mean == 0.00001))
# exprSet = exprSet[,-which(names(exprSet)%in%list)] ## 删除表达值为0的矩阵
colnames(exprSet)
log_rank_p <- apply(exprSet[,4:length((exprSet))], 2, function(values1)
{ 
  # values1 %>% drop_na() ## delletion NA
  # group = ifelse(values1 > 0.3 ,'High','Low')
  group = ifelse(values1 > median(values1) ,'Del',"WT")
  # ifelse(values1< - 0.3,'Del','WT'))
  # ifelse(values1 > 0 ,'Amp','Del'))
  kmfit2 <- survfit(my.surv~group, data=exprSet)
  #plot(kmfit2)
  data.survdiff = survdiff(my.surv~group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
})


## 运行完了之后返回的找出小于0.05的P.val，返回dataframe
# log_rank_p <- log_rank_p[log_rank_p< 0.05]
names(log_rank_p[log_rank_p < 0.05])
log_rank_p = log_rank_p[log_rank_p < 1]
df = data.frame(log_rank_p)

write.csv(df,'/Users/ranpeng/Desktop/CCA/Data/Fig7/Xcell_signature_cor/Neutrophils/heatmap_O_link_glycos_os_results.csv')

## --------------- 返回生存P-value和HR的总表 ----------------
survival_dat = exprSet
my.surv <- Surv(survival_dat$time, survival_dat$status)
cox_results <- apply(survival_dat[,4:length(survival_dat)], 2, function(values1){
  group = ifelse(values1 > median(values1) ,'xAmp',"Low")
  survival_dat <- data.frame(group=group,stringsAsFactors = F)
  m=coxph(my.surv ~ group, data =  survival_dat)
  
  beta <- coef(m)
  se <- sqrt(diag(vcov(m)))
  HR <- exp(beta)
  HRse <- HR * se
  
  #summary(m)
  tmp <- round(cbind(coef = beta, se = se, z = beta/se, p = 1 - pchisq((beta/se)^2, 1),
                     HR = HR, HRse = HRse,
                     HRz = (HR - 1) / HRse, HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                     HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                     HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
  return(tmp['groupxAmp',])
  
})


## ------------- 筛选显著的基因 -----------------
(df = cox_results[,cox_results[4,]<0.05])
df_surv = data.frame(t(cox_results))
write.csv(df_surv,'/Users/ranpeng/Desktop/CCA/Data/Fig7/Xcell_signature_cor/Neutrophils/heatmap_O_link_glycos_os_HR.csv')

## ------------- HR-cox function -----------------
df = CCA
covariates <- colnames(df)[9:length(df)]
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(time, status)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = df)})
# Extract data
univ_results <- lapply(  univ_models,
                         function(x){
                           x <- summary(x)
                           p.value<-signif(x$wald["pvalue"], digits=2)
                           wald.test<-signif(x$wald["test"], digits=2)
                           beta<-signif(x$coef[1], digits=2);#coeficient beta
                           HR <-signif(x$coef[2], digits=2);#exp(beta)
                           HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                           HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                           HR <- paste0(HR, " (",
                                        HR.confint.lower, "-", HR.confint.upper, ")")
                           res<-c(beta, HR, wald.test, p.value)
                           names(res)<-c("beta", "HR (95% CI for HR)", "wald.test",
                                         "p.value")
                           return(res)
                           #return(exp(cbind(coef(x),confint(x))))
                         })
res <- t(as.data.frame(univ_results, check.names = FALSE))
HR =as.data.frame(res)
write.csv(HR,'CCA_protein_OS_HR_Continuous.csv')


## 对三组数据进行整合
CCA_os = read_csv("/Users/ranpeng/OneDrive - mail.ecust.edu.cn/项目文件/CCA-project/Data/Fig3/protein_os/results/197_tumor/CCA_protein_OS_HR_median.csv")
iCC_os = read_csv("/Users/ranpeng/OneDrive - mail.ecust.edu.cn/项目文件/CCA-project/Data/Fig3/protein_os/results/197_tumor/iCC_protein_OS_HR_median.csv")
eCC_os = read_csv("/Users/ranpeng/OneDrive - mail.ecust.edu.cn/项目文件/CCA-project/Data/Fig3/protein_os/results/197_tumor/eCC_protein_OS_HR_median.csv")
CCA_os = CCA_os %>% select(X1, p, HR) %>% rename(CCA_pval = p,CCA_HR = HR)
iCC_os = iCC_os %>% select(X1, p, HR) %>% rename(iCC_pval= p,iCC_HR= HR)
eCC_os = eCC_os %>% select(X1, p, HR) %>% rename(eCC_pval= p,eCC_HR= HR)

integra_os = merge(CCA_os,iCC_os)
integra_os = merge(integra_os,eCC_os, by = "X1")
colnames(integra_os)
# integra_os = read.csv("/Users/ranpeng/OneDrive - mail.ecust.edu.cn/项目文件/CCA-project/Data/Fig3/protein_os/results/217_tumor/protein_integra_os_1.csv")
## 筛选数据 ，与生存相关并且蛋白高表达和预后更差相关
rownames(integra_os) = integra_os[,1]
integra_os = integra_os %>% select(!Gene)
filter_os = integra_os%>% filter( CCA_pval < 0.05 & CCA_HR > 1.2 &
                                    iCC_pval < 0.05 & iCC_HR > 1.2 &
                                    eCC_pval < 0.05 & eCC_HR > 1.2) 

write.csv(filter_os,'protein_filter_os.csv')

## 筛选数据
# DT::datatable(cox_results ,
#               extensions = 'FixedColumns',
#               options = list(
#                 #dom = 't',
#                 scrollX = TRUE,
#                 fixedColumns = TRUE
#               ))




