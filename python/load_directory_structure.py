#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Manuscript  : Integrated Proteogenomic Characterization of Cholangiocarcinoma
@File     : load_directory_structure.py
@Time     : 2022/03/11 16:32:19
@Author   : RanPeng 
@Version  : python 3.8.6
@Contact  : 2502388440@hotmail.com
@License  : (C)Copyright 2021-2022, DingLab-CHINA-SHANGHAI
@Description : Loads the directory structure for analysis
'''
# here put the import lib
import os
import sys
#=======================================================================================
def get_directory_structure(base_dir = None):
    baseDir = base_dir
    os.chdir(baseDir)
    suppTableDir = baseDir + "supplementary_tables/"
    outputDir    = baseDir + "figure_subpanels/"
    outTableDir  = baseDir + "output_tables/"
    resourcesDir = baseDir + "resources/"
    dataDir      = baseDir + "data/"
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    else:
        print(f"The outputDir already exists: {outputDir} ")
    if not os.path.exists(outTableDir):
        os.makedirs(outTableDir)
    else:
        print(f"The outTableDir already exists: {outTableDir} ")
    if not os.path.exists(suppTableDir):
        os.makedirs(suppTableDir)
    else:
        print(f"The suppTableDir already exists: {suppTableDir} ")
    if not os.path.exists(resourcesDir):
        os.makedirs(resourcesDir)
    else:
        print(f"The resourcesDir already exists: {resourcesDir} ")
    if not os.path.exists(dataDir):
        os.makedirs(dataDir)
    else:
        print(f"The dataDir already exists: {dataDir} ")
    return outputDir, outTableDir, suppTableDir, resourcesDir, dataDir