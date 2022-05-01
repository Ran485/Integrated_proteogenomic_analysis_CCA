#' #!/usr/bin/env Rscript
#' # -*- encoding: utf-8 -*-
#' """
#' @Manuscript  : Integrated Proteogenomic Characterization of Cholangiocarcinoma
#' @Time        : 2022/03/11 10:57:18
#' @Author      : RanPeng
#' @Version     : R version 4.1.0
#' @Contact     : 2502388440@hotmail.com
#' @License     : (C)Copyright 2021-2022, DingLab-CHINA-SHANGHAI
#' @Description : Loads the directory structure for analysis
#' """
#=======================================================================================

suppTableDir <- paste0(baseDir,"/supplementary_tables/")
outputDir    <- paste0(baseDir,"/figure_subpanels/")
outTableDir  <- paste0(baseDir,"/output_tables/")
resourcesDir <- paste0(baseDir,"/resources/")
dataDir      <- paste0(baseDir,"/data/")

source(paste0(baseDir,"/R/manuscript_theme.R"))

setwd(baseDir)

