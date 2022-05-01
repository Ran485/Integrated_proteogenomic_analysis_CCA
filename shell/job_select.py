#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File     : job_select.py
@Time     : 2022/03/30 23:34:52
@Author   : RanPeng 
@Version  : python 3.8.6
@Contact  : 2502388440@hotmail.com
@License  : (C)Copyright 2021-2022, DingLab-CHINA-SHANGHAI
@Function : None
'''
# here put the import lib
import pandas as pd
import numpy as np

job_success = pd.read_csv("/Users/ranpeng/Desktop/CCA/CCA_script/python_script/gene_fusion/job_success.csv", header=None, index_col=0)
job_sample = pd.read_csv("/Users/ranpeng/Desktop/CCA/CCA_script/python_script/gene_fusion/job_samples.csv", header=None, index_col=0)

ret = list(set(job_sample.index) ^ set(job_success.index))
job_selected = pd.DataFrame(ret)
print(job_selected)
job_selected.to_csv("/Users/ranpeng/Desktop/CCA/CCA_script/python_script/gene_fusion/job_samples_fail.csv", header=None, index=None)