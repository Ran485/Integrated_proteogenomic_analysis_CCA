# install.packages("table1")
# install.packages("boot") 示例数据
library(table1) 
library(boot)
library(survival)
library(car)
options(stringsAsFactors = F)

#重写函数
#改写赋给table1中render参数的函数
rd <- function(x, name, ...){
  y <- a[[name]]
  m = any(is.na(y))
  if (is.numeric(y)) {
    normt=(shapiro.test(y)$p.value>0.1)+2
    my.render.cont(x,n=normt,m=m,...)
  }else{
    render.default(x=x, name=name, ...)
  }
}

#改写赋给table1中render.continuous参数的函数
my.render.cont <- function (x,n,m, ...) {
  
  a=with(stats.apply.rounding(stats.default(x, ...), ...), c("", 
                                                             `Mean (±SD)` = sprintf("%s (±%s)", MEAN, SD),
                                                             `Median (IQR)` = sprintf("%s [%s, %s]", MEDIAN, Q1, Q3)
                                                             # `Median (Min, Max)` = sprintf("%s (%s, %s)", MEDIAN, MIN, MAX)
  ))[-n]
  if (m) {
    a <- c(a,with(stats.apply.rounding(stats.default(is.na(x), ...), ...)$Yes,
                  c(`Missing` = sprintf("%s (%s%%)", FREQ, PCT))))
  }
  a
}

#增加新列：正态性检验
normality <- function(x, ...) {
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  l <- length(x)
  
  if (is.numeric(y)) {
    #正态性检验p>0.1的水准
    zt <- c()
    for (i in 1:length(table(g))) {
      zt <- c(zt,shapiro.test(y[g==i])$p.value>0.1)
    }
    ztjy <- sum(zt)==length(zt)
    #方差齐性检验p>0.1的水准
    fcq <- leveneTest(y~g,center=median)$`Pr(>F)`[1]>0.1
    normality <- c("非正态","正态")[ztjy+1]
    
  } else {
    # For categorical variables, perform a chi-squared test of independence
    normality <- "-"
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  #sub("<", "&lt;", format.pval(p, digits=3, eps=0.001))
  normality
}

#增加新列：方差齐性检验
hom.var <- function(x, ...) {
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  l <- length(x)
  
  if (is.numeric(y)) {
    #正态性检验p>0.1的水准
    zt <- c()
    for (i in 1:length(table(g))) {
      zt <- c(zt,shapiro.test(y[g==i])$p.value>0.1)
    }
    ztjy <- sum(zt)==length(zt)
    #方差齐性检验p>0.1的水准
    fcq <- leveneTest(y~g,center=median)$`Pr(>F)`[1]>0.1
    hom_of_var <- c("不齐","齐")[fcq+1]
    
  } else {
    # For categorical variables, perform a chi-squared test of independence
    hom_of_var <- "-"
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  #sub("<", "&lt;", format.pval(p, digits=3, eps=0.001))
  hom_of_var
}

#增加新列：p值
p <- function(x,result, ...) {
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  l <- length(x)
  
  if (is.numeric(y)) {
    #正态性检验p>0.1的水准
    zt <- c()
    for (i in 1:length(table(g))) {
      zt <- c(zt,shapiro.test(y[g==i])$p.value>0.1)
    }
    ztjy <- sum(zt)==length(zt)
    #方差齐性检验p>0.1的水准
    fcq <- leveneTest(y~g,center=mean)$`Pr(>F)`[1]>0.1
    ####（正态方差齐）####
    if (ztjy&fcq) {
      if (length(table(g))>2) {
        # 多个样本均数的方差分析:其实两个样本的方差分析和t检验是等价的
        p <- summary(aov(y~g))[[1]]$`Pr(>F)`[1]
        p.method <- "test:ANOVA"
      }else{
        # 两个样本均数的t检验
        p <- t.test(y ~ g)$p.value
        p.method <- "test:t"
      }
    }else{
      #####（非正态分布或方差不齐的多个样本）###
      
      if (length(table(g))>2) {
        #多个样本均数比较
        p <- kruskal.test(y~g)$p.value 
        p.method <- "test:kruskal"
      }else{
        #两个样本均属比较
        p <- wilcox.test(y~g)$p.value
        p.method <- "test:wilcox"
      }
    }
  }else{
    # For categorical variables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g))$p.value
    p.method <- "test:chi-square"
  }
  #格式化p值
  c(sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
  
}

#增加新列：统计方法
pmethod <- function(x, ...) {
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  l <- length(x)
  
  if (is.numeric(y)) {
    #正态性检验p>0.1的水准
    zt <- c()
    for (i in 1:length(table(g))) {
      zt <- c(zt,shapiro.test(y[g==i])$p.value>0.1)
    }
    ztjy <- sum(zt)==length(zt)
    #方差齐性检验p>0.1的水准
    fcq <- leveneTest(y~g,center=mean)$`Pr(>F)`[1]>0.1
    ####（正态方差齐分布）####
    if (ztjy&fcq) {
      if (length(table(g))>2) {
        # 多个样本均数的方差分析:其实两个样本的方差分析和t检验是等价的
        #p <- summary(aov(y~g))[[1]]$`Pr(>F)`[1]
        p.method <- "ANOVA"
      }else{
        # 两个样本均数的t检验
        #p <- t.test(y ~ g)$p.value
        p.method <- "t"
      }
    }else{
      #####（非正态分布或方差不齐的多个样本）###
      
      if (length(table(g))>2) {
        #多个样本均数比较
        #p <- kruskal.test(y~g)$p.value 
        p.method <- "Kruskal"
      }else{
        #两个样本均属比较
        #p <- wilcox.test(y~g)$p.value
        p.method <- "Wilcox"
      }
      
    }
    
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g))$p.value
    p.method <- "Chi-square"
  }
  p.method
}